#include <stdio.h>
#include <stdlib.h>
#include <mpreal.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <typeinfo>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <plplot/plplot.h>
#include <plplot/plstream.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <fstream>

#include "common/shader.hpp"

#include "Scene.h"
#include "particles.h"
#include <string>
#include <thread>
std::string times;


class x00 {
public:
    x00( int, const char ** );

private:
    // Class data
    static const int NSIZE;
};



const mpfr::mpreal k_e = 8.9875517873681764*pow(10,9);
const mpfr::mpreal delta_t = 1.667*pow(10,-24)/sqrt(3) * pow(10, -0);
const mpfr::mpreal G =  6.673*pow(10, -11);
const mpfr::mpreal H =  6.62606957*pow(10, -34);
const mpfr::mpreal C =  2.99792458*pow(10, 8);

bool paused = false;

extern float viewMatrix[16];

std::vector<Particle *> all_particles = std::vector<Particle* >();
std::vector<int> keys_pressed = std::vector<int>();
// new Electron(), new Proton()


GLFWwindow* createWindow(){
	if( !glfwInit() )
	{
		fprintf( stderr, "Failed to initialize GLFW\n" );
		return NULL;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint (GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint (GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint (GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint (GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    return glfwCreateWindow(1024, 640, "matrix", NULL, NULL);
}

static void onError(int error, const char* description)
{
    std::cout << "Error: " << description << std::endl;
}

Scene scene = Scene();

void onFramebufferResize(GLFWwindow* window, int width, int height)
{
    scene.reshape(width, height);
    glfwSwapBuffers(window);

}

mpfr::mpreal my_gamma(mpfr::mpreal velocity){
	return 1.0f/sqrt(1 - pow(velocity/C, 2) );
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    static bool depthClampingActive = false;

    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GL_TRUE);
    } else if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
        if(depthClampingActive) {
            glDisable(GL_DEPTH_CLAMP);
        } else {
            glEnable(GL_DEPTH_CLAMP);
        }
        depthClampingActive = !depthClampingActive;
    } else if ( action == GLFW_PRESS || action == GLFW_RELEASE) {
        scene.keyStateChanged(key, action);
    }
}

float a_var = 15.5f;
float b_var = 15.4;
float c_var = 15.2;

void compute_new_values() {
	// Compute new position
	for (uint i = 0; i < all_particles.size(); i++ ){
		Particle* current = all_particles[i];

		current->position[0] += current->velocity[0]*delta_t+0.5*current->acceleration[0]*delta_t*delta_t;
		current->position[1] += current->velocity[1]*delta_t+0.5*current->acceleration[1]*delta_t*delta_t;
		current->position[2] += current->velocity[2]*delta_t+0.5*current->acceleration[2]*delta_t*delta_t;
		current->position_back_log.push_back( current->position );
		if (current->position_back_log.size()> 2000){ current->position_back_log.pop_front(); }

		/*
		std::cout << ">>" << std::endl;
		std::cout << current->position[0] << std::endl;
		std::cout << current->position[1] << std::endl;
		std::cout << current->position[2] << std::endl;
		*/
	}

	for (uint i = 0; i < all_particles.size(); i++ ){ // calculate new values
		Particle* current = all_particles[i];
		if (current->energy != 0){
			continue;
		}

		std::vector<mpfr::mpreal> sum_of_all_forces = std::vector<mpfr::mpreal>{0,0,0};

		for (uint j = 0; j < all_particles.size(); j++ ){ // calculate new acceleration values


			if (i==j){ continue; }

			Particle* other = all_particles[j];
			if(other->energy != 0){
				continue;
			}

			// calculate
			// all_particles[i]->next_position

			mpfr::mpreal dist = sqrt(
				  pow(current->position[0]-other->position[0], 2)
				+ pow(current->position[1]-other->position[1], 2)
				+ pow(current->position[2]-other->position[2], 2) );

			// Formula to be optimized (simplified).

			mpfr::mpreal cl_scalar = - k_e * current->charge * other->charge / dist / dist / dist;
			mpfr::mpreal gravity_scalar = G * current->mass * other-> mass / dist /dist / dist;

			//mpfr::mpreal weak_scalar = -(pow(10, 15) * exp(-dist*pow(10,15) ) )/dist - (exp(-dist*pow(10,15) ))/ dist / dist;
			mpfr::mpreal weak_scalar = 	(-1 * pow(10, 15.2f)) * (exp(-dist * pow(10, 15.2f))) / dist
										- (exp(-dist * pow(10, 15.2f))) / pow(dist, 2);

			mpfr::mpreal A = (pow(10, a_var) * exp(-dist * pow(10, b_var) ) )/ dist;
			mpfr::mpreal B =  exp(-dist * pow(10, c_var) )/ pow(dist, 2);
			mpfr::mpreal A_and_B = A + B;



			std::cout << ">>>" << std::endl;
			std::cout << cl_scalar << std::endl;
			std::cout << weak_scalar << std::endl;


			std::vector<mpfr::mpreal> delta_vector = std::vector<mpfr::mpreal>{
				other->position[0]-current->position[0],
				other->position[1]-current->position[1],
				other->position[2]-current->position[2]};

			std::vector<mpfr::mpreal> new_force = std::vector<mpfr::mpreal>{
				(A_and_B + weak_scalar + gravity_scalar + cl_scalar) *delta_vector[0],
				(A_and_B + weak_scalar + gravity_scalar + cl_scalar) *delta_vector[1],
				(A_and_B + weak_scalar + gravity_scalar + cl_scalar) *delta_vector[2]};
			std::cout << new_force[0] << std::endl;

			sum_of_all_forces[0] += new_force[0];
			sum_of_all_forces[1] += new_force[1];
			sum_of_all_forces[2] += new_force[2];

			std::cout << sum_of_all_forces[0] << std::endl;


			/*
			std::cout << all_particles[i]->position[0] << " "
				<< all_particles[i]->position[1] << " "
				<< all_particles[i]->position[2] << " "
			 << std::endl;
			 */
		}
		std::vector<mpfr::mpreal> new_acceleration = std::vector<mpfr::mpreal>{
			sum_of_all_forces[0]/(current->mass*my_gamma( current->velocity[0] )),
			sum_of_all_forces[1]/(current->mass*my_gamma( current->velocity[1] )),
			sum_of_all_forces[2]/(current->mass*my_gamma( current->velocity[2] ))};


		current->last_acceleration = current->acceleration;
		current->acceleration = new_acceleration;
	}

	for (uint i = 0; i < all_particles.size(); i++ ){
		Particle* current = all_particles[i];
		if(current->energy != 0){
			continue;
		}

		std::vector<mpfr::mpreal> initial_velocity = std::vector<mpfr::mpreal> {
			current->velocity[0],
			current->velocity[1],
			current->velocity[2],
		};

		// Leap Frog
		std::vector<mpfr::mpreal> final_velocity = std::vector<mpfr::mpreal> {
			initial_velocity[0] + (current->acceleration[0]+current->last_acceleration[0])/2*delta_t,
			initial_velocity[1] + (current->acceleration[1]+current->last_acceleration[1])/2*delta_t,
			initial_velocity[2] + (current->acceleration[2]+current->last_acceleration[2])/2*delta_t
		};
		current->velocity = final_velocity;

		/*
		// equation photo: http://goo.gl/J1d2yP
		for (uint j = 0; j < 3; j++){
			int dir = 1;
			if(current->velocity[j] < 0){
				dir = -1;
			}

			mpfr::mpreal gamma_factor = my_gamma( initial_velocity[j] );

			mpfr::mpreal mass_initial = current->mass * gamma_factor;

			mpfr::mpreal relativistic_energy = current->mass * pow(C, 2) * (gamma_factor - 1);

			mpfr::mpreal delta_v = final_velocity[j] - initial_velocity[j];

			mpfr::mpreal delta_energy = 1.0f/2.0f * mass_initial * pow(delta_v, 2);

			mpfr::mpreal energy_final = delta_energy + relativistic_energy;

			mpfr::mpreal velocity_final = C * sqrt(
				1 - 1/pow(
					energy_final / current->mass / pow(C, 2) + 1
				, 2)
			);



			current->velocity[j] = velocity_final;
		}

		for (uint j = 0; j < 3; j++){
			current->kinetic_energy[j] = 1.0f/2.0f * current->mass * pow(current->velocity[j], 2) * (my_gamma(current->velocity[j])-1);
		}
		*/
	}
}

/*
float new_perspective_matrix[16];
float perspectiveMatrix[16];
*/

glm::mat4 identityMatrix = {1,0,0,0, // keep in mind that this identity matrix
	  						0,1,0,0, // has been created in row-major mode.
							0,0,1,0, // but will be used in column-major mode.
							0,0,0,1};

void passDataToArray(glm::mat4 matrix, float* my_array) {
	const float *pSource = (const float*)glm::value_ptr(matrix);
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			my_array[j+4*i] = matrix[i][j];
		}
	}
}

glm::mat4 move(glm::mat4 mat, float x, float y, float z) {

	glm::mat4 newMatrix = identityMatrix;

	newMatrix[3][0] = x; // x
	newMatrix[3][1] = y; // y
	newMatrix[3][2] = z; // z

	glm::mat4 out = newMatrix * mat;

	return out;
}

glm::mat4 scale(glm::mat4 mat, float x, float y, float z) {
	glm::mat4 newMatrix = identityMatrix;

	newMatrix[0][0] = x; // x
	newMatrix[1][1] = y; // y
	newMatrix[2][2] = z; // z

	glm::mat4 out = newMatrix * mat;

	return out;
}

glm::mat4 rotate(glm::mat4 mat, float y, float x, float z) {
	glm::mat4 out = mat;

	glm::mat4 newMatrix;

	newMatrix = identityMatrix;
	newMatrix[1][1] = cos(x);
	newMatrix[1][2] = sin(x);
	newMatrix[2][2] = cos(x);
	newMatrix[2][1] = -sin(x);

	out = newMatrix * out;

	newMatrix = identityMatrix;
	newMatrix[0][0] = cos(y);
	newMatrix[0][2] = -sin(y);
	newMatrix[2][2] = cos(y);
	newMatrix[2][0] = sin(y);

	out = newMatrix * out;

	newMatrix = identityMatrix;
	newMatrix[0][0] = cos(z);
	newMatrix[0][1] = sin(z);
	newMatrix[1][0] = -sin(z);
	newMatrix[1][1] = cos(z);

	out = newMatrix * out;

	return out;
}

x00* x;
plstream* pls;
void init_graphs(int argc, const char * argv[]){
    x = new x00(argc, argv );
}

const int x00::NSIZE = 101;

x00::x00( int argc, const char **argv ){
    PLFLT x[NSIZE], y[NSIZE];
    PLFLT xmin = 0., xmax = 1.0, ymin = 0., ymax = 10000.;
    int   i;

    // Prepare data to be plotted.
    for ( i = 0; i < NSIZE; i++ )
    {
        x[i] = (PLFLT) ( i ) / (PLFLT) ( NSIZE - 1 );
        y[i] = ymax * x[i] * x[i];
    }

    pls = new plstream();

    // Parse and process command line arguments
    pls->parseopts( &argc, argv, PL_PARSE_FULL );

    // Initialize plplot
	pls->sdev("qtwidget");
    pls->init();

    // Create a labelled box to hold the plot.
    pls->env( xmin, xmax, ymin, ymax, 0, 0 );
    pls->lab( "x", "y=100 x#u2#d", "Simple PLplot demo of a 2D line plot" );

	std::cout << pls << std::endl;
}

PLFLT * xs = new PLFLT[1];
PLFLT * ys = new PLFLT[1];

PLFLT maxt = 100;
PLFLT maxe = 100;
PLFLT currt = 0;

mpfr::mpreal total_energy;
void display_graphs() {
	currt += 0.0001f;
	xs[0] = currt;
	ys[0] = -((float)total_energy)*10e19;

	pls->poin(1, xs, ys, 1);
	pls->flush();
}

void create_photons(){
	for(uint i = 0; i < all_particles.size(); i++){
		long double num = (std::rand() / (float)RAND_MAX);
		if( num > 0.9) {
			if( all_particles[i]->charge != 0 ){
				mpfr::mpreal x = (std::rand() / RAND_MAX);
				mpfr::mpreal y = (std::rand() / RAND_MAX);
				mpfr::mpreal z = (std::rand() / RAND_MAX);
				mpfr::mpreal k = sqrt( (long double) (x*x + y*y + z*z) )/C*100000;
				std::vector<mpfr::mpreal> velocity = {x/k, y/k, z/k};

				all_particles.push_back(new Photon());
				all_particles.back()->velocity = velocity;
				all_particles.back()->position = all_particles[i]->position;
				all_particles.back()->energy = 1.0;

				std::cout << all_particles.size() << std::endl;
			}
		}
	}
}

// 1.07 * 10-17


void compute_user_input() {
	//new_perspective_matrix = perspectiveMatrix;

	glm::mat4 outputMatrix = identityMatrix;
	const float theta = 0.1f;

	for (uint i = 0; i < keys_pressed.size(); i++){
		glm::mat4 newMatrix = identityMatrix;
		int key = keys_pressed[i];
		switch(key){
			case GLFW_KEY_W:
				newMatrix = move(newMatrix, 0, 0, 0.1f);
				break;
			case GLFW_KEY_S:
				newMatrix = move(newMatrix, 0, 0, -0.1f);
				break;
			case GLFW_KEY_A:
				newMatrix = move(newMatrix, 0.1f, 0, 0.0f);
				break;
			case GLFW_KEY_D:
				newMatrix = move(newMatrix, -0.1f, 0, 0);
				break;
			case GLFW_KEY_R:
				newMatrix = move(newMatrix, 0, -0.1f, 0);
				break;
			case GLFW_KEY_F:
				newMatrix = move(newMatrix, 0, 0.1f, 0);
				break;
			case GLFW_KEY_SPACE:
				if(!paused){
					paused = true;
				}else{
					paused = false;
				}
				break;
			case GLFW_KEY_UP:
				newMatrix = rotate(newMatrix, 0, -theta, 0);
				break;
			case GLFW_KEY_DOWN:
				newMatrix = rotate(newMatrix, 0, theta, 0);
				break;
			case GLFW_KEY_LEFT:
				newMatrix = rotate(newMatrix, -theta, 0, 0);
				break;
			case GLFW_KEY_RIGHT:
				newMatrix = rotate(newMatrix, theta, 0, 0);
				break;
		}
		outputMatrix = newMatrix*outputMatrix;
	}

	glm::mat4 r = outputMatrix * glm::make_mat4(viewMatrix);

	glm::mat4 mat = glm::mat4(1.0f);
	//mat = move(r, 1.0f, 0.0f, 0.0f);

	passDataToArray(r, viewMatrix);

	scene.updateCamera();
}

mpfr::mpreal distance(Particle * first, Particle * second){
	return sqrt(
		pow(first->position[0]-second->position[0], 2) +
		pow(first->position[1]-second->position[1], 2) +
		pow(first->position[2]-second->position[2], 2)
			);
}
mpfr::mpreal norm(std::vector<mpfr::mpreal> vec){
	return sqrt(
		pow(vec[0], 2) +
		pow(vec[1], 2) +
		pow(vec[2], 2)
			);
}

void check_conditions(){
	// check if any photons are crossing any electrons
	Particle * current;
	Particle * other;
	for(uint i = 0; i < all_particles.size(); i++){
		current = all_particles[i];
		for(uint j = 0; j < all_particles.size(); j++){
			other = all_particles[j];
			if (i==j){ continue; }

			if( distance(current, other) < current->radius ) {
				// other interacts with current
			}
		}
	}
}

void save_everything(){
	std::string filename = "local/" + times;
	std::cout << filename << std::endl;
	std::ofstream ofs(filename, std::ofstream::app);

	std::string version = "5";

	{
		boost::archive::text_oarchive oa(ofs);
		ofs << "version:" << version << std::endl;
		oa << all_particles;
		ofs << std::endl;
		ofs.close();
	}
}

void load_everything(std::string timestamp, int version){
	std::string filename = "local/123";
	std::ifstream ifs(filename);
}

void check_energy(){
	Particle * current;
	Particle * other;
	mpfr::mpreal U_g = 0;
	mpfr::mpreal U_e = 0;
	mpfr::mpreal K = 0;
	mpfr::mpreal E = 0;
	for(uint i = 0; i < all_particles.size(); i++){
		current = all_particles[i];
		for(uint j = 0; j < all_particles.size(); j++){
			other = all_particles[j];
			if(other == current){ continue; }
			if(j<=i){ continue; }

			mpfr::mpreal dist = distance(current, other);

			//U_g += -G * current->mass * other->mass / dist;

			if( current->energy == 0 && other->energy == 0){
				U_e += k_e * current->charge * other->charge / dist;
			}

		}

		if( current->energy == 0 ){
			K += 1/2.0 * current->mass * pow(norm(current->velocity), 2);
		} else{
				E += current->energy;
		}
	}

	total_energy = U_g + U_e + K + E;
}

int openglmain(int argc, const char * argv[]){
	glfwSetErrorCallback(onError);

    GLFWwindow* window = createWindow();
    if (!window)
        return 0;

    glfwMakeContextCurrent(window);

    scene.init();
	init_graphs(argc, argv);

    int windowWidth = 0;
    int windowHeight = 0;

    glfwGetFramebufferSize(window, &windowWidth, &windowHeight);
    onFramebufferResize(window, windowWidth, windowHeight);
    glfwSetFramebufferSizeCallback(window, &onFramebufferResize);

    glfwSetKeyCallback(window, key_callback);

    while (!glfwWindowShouldClose(window))
    {
		std::cout << "New iteration" << std::endl;

		if( paused ){
			scene.draw();
			continue;
		}

		std::thread compute (compute_new_values);

		std::thread check_cond(check_conditions);
		std::thread check_E(check_energy);
		compute_user_input();

		std::thread show_graph(display_graphs);

		//create_photons();

		std::thread save_log(save_everything);

		//std::thread draw(scene.draw);
		scene.draw();

		compute.join();
		check_cond.join();
		check_E.join();
		show_graph.join();
		save_log.join();
		//draw.join();


        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();

    return 0;
}

int main(int argc, const char * argv[]){
	srand (time(NULL));

	time_t timev;
	time(&timev);
	times = std::to_string(timev);

	// PARTEYYYY
	// We are working with `digits` long numbers now hehe.
    const int digits = 500;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));

	std::cout.precision(digits);

	mpfr::mpreal r_not = 5.29 * pow(10,-11);
	mpfr::mpreal v_not = 2.18805743462617 * pow(10,6);
	mpfr::mpreal new_dist = 0.29 * pow(10,-13);

	mpfr::mpreal r_nuke = 0.46f * 4 * pow(10,-14);



	all_particles.push_back(new Proton());
	all_particles.push_back(new Proton());
	all_particles.push_back(new Neutron());


	all_particles[0]->position = std::vector<mpfr::mpreal>{0, 0, 0};
	all_particles[0]->velocity = std::vector<mpfr::mpreal>{0, 0, 0};

	all_particles[1]->position = std::vector<mpfr::mpreal>{r_nuke, 0, 0};
	all_particles[1]->velocity = std::vector<mpfr::mpreal>{0, v_not, 0};

	all_particles[2]->position = std::vector<mpfr::mpreal>{r_nuke/1.4f, r_nuke/1.1f, 0};
	all_particles[2]->velocity = std::vector<mpfr::mpreal>{0, 0, 0};


	for (uint i = 0; i < all_particles.size(); i++ ){ // calculate new values
		Particle* current = all_particles[i];
		if(current->energy != 0){
			continue;
		}

		for (uint j = 0; j < all_particles.size(); j++ ){ // calculate new acceleration values
			if (i==j){ continue; }

			Particle* other = all_particles[j];

			// calculate
			// all_particles[i]->next_position

			mpfr::mpreal dist = sqrt(
				  pow(current->position[0]-other->position[0], 2)
				+ pow(current->position[1]-other->position[1], 2)
				+ pow(current->position[2]-other->position[2], 2) );

			mpfr::mpreal cl_scalar = - k_e * current->charge * other->charge / dist / dist / current->mass / dist;
			cl_scalar = 0.0f;

			std::vector<mpfr::mpreal> delta_vector = std::vector<mpfr::mpreal>{
				other->position[0]-current->position[0],
				other->position[1]-current->position[1],
				other->position[2]-current->position[2]};

			std::vector<mpfr::mpreal> new_acceleration = std::vector<mpfr::mpreal>{
				cl_scalar*delta_vector[0],
				cl_scalar*delta_vector[1],
				cl_scalar*delta_vector[2]};

			current->last_acceleration = current->acceleration;
			current->acceleration = new_acceleration;
		}
	}

	openglmain(argc, argv);

	// Destroy the graphs window setup.
	//delete pls;
    //delete x;

	return 0;
}

