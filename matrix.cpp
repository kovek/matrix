#include <vector>
#include <iostream>
#include <cmath>
#include <typeinfo>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "common/shader.hpp"

#define GLM_FORCE_RADIANS
#include <glm/glm/glm.hpp>
#include <glm/glm/gtc/matrix_transform.hpp>
#include <glm/glm/gtc/type_ptr.hpp>

#include "Scene.h"
#include "particles.h"



const long double k_e = 8.9875517873681764*pow(10,9);
const long double delta_t = 1.52l*pow(10, -16)/2500l;
const long double G =  6.673*pow(10, -11);
const long double H =  6.62606957*pow(10, -34);
const long double C =  2.99792458*pow(10, 8);

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


void compute_new_values() {
	for (uint i = 0; i < all_particles.size(); i++ ){
		Particle* current = all_particles[i];

		current->position[0] += current->velocity[0]*delta_t+0.5*current->acceleration[0]*delta_t*delta_t;
		current->position[1] += current->velocity[1]*delta_t+0.5*current->acceleration[1]*delta_t*delta_t;
		current->position[2] += current->velocity[2]*delta_t+0.5*current->acceleration[2]*delta_t*delta_t;
		current->position_back_log.push_back( current->position );
		if (current->position_back_log.size()> 2000){ current->position_back_log.pop_front(); }

	}


	for (uint i = 0; i < all_particles.size(); i++ ){ // calculate new values
		Particle* current = all_particles[i];


		std::vector<long double> sum_of_all_forces = std::vector<long double>{0,0,0};

		for (uint j = 0; j < all_particles.size(); j++ ){ // calculate new acceleration values


			if (i==j){ continue; }

			Particle* other = all_particles[j];

			if(current->wavelength != 0 || other->wavelength != 0){ // that means we're dealing with a photon. don't touch its acceleration...
				continue;
			}

			// calculate
			// all_particles[i]->next_position

			long double dist = sqrt(
				  pow(current->position[0]-other->position[0], 2)
				+ pow(current->position[1]-other->position[1], 2)
				+ pow(current->position[2]-other->position[2], 2) );

			long double cl_scalar = - k_e * current->charge * other->charge / dist / dist / dist;

			std::vector<long double> delta_vector = std::vector<long double>{
				other->position[0]-current->position[0],
				other->position[1]-current->position[1],
				other->position[2]-current->position[2]};

			std::vector<long double> new_force = std::vector<long double>{
				cl_scalar*delta_vector[0],
				cl_scalar*delta_vector[1],
				cl_scalar*delta_vector[2]};

			sum_of_all_forces[0] += new_force[0];
			sum_of_all_forces[1] += new_force[1];
			sum_of_all_forces[2] += new_force[2];


			/*
			std::cout << all_particles[i]->position[0] << " "
				<< all_particles[i]->position[1] << " "
				<< all_particles[i]->position[2] << " "
			 << std::endl;
			 */
		}
		std::vector<long double> new_acceleration = std::vector<long double>{
			sum_of_all_forces[0]/current->mass,
			sum_of_all_forces[1]/current->mass,
			sum_of_all_forces[2]/current->mass};


		current->last_acceleration = current->acceleration;
		current->acceleration = new_acceleration;
	}

	for (uint i = 0; i < all_particles.size(); i++ ){
		Particle* current = all_particles[i];
		current->velocity[0] += (current->acceleration[0]+current->last_acceleration[0])/2*delta_t;
		current->velocity[1] += (current->acceleration[1]+current->last_acceleration[1])/2*delta_t;
		current->velocity[2] += (current->acceleration[2]+current->last_acceleration[2])/2*delta_t;
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

long double distance(Particle * first, Particle * second){
	return sqrt(
		pow(first->position[0]-second->position[0], 2) +
		pow(first->position[1]-second->position[1], 2) +
		pow(first->position[2]-second->position[2], 2)
			);
}
long double norm(std::vector<long double> vec){
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

void check_energy(){
	Particle * current;
	Particle * other;
	long double U_g = 0;
	long double U_e = 0;
	long double K = 0;
	long double E = 0;
	for(uint i = 0; i < all_particles.size(); i++){
		current = all_particles[i];
		for(uint j = 0; j < all_particles.size(); j++){
			other = all_particles[j];
			if(other == current){ continue; }
			if(j<=i){ continue; }

			long double dist = distance(current, other);

			//U_g += -G * current->mass * other->mass / dist;
			U_e += k_e * current->charge * other->charge / dist;

			if( typeid(current) == typeid(Photon *) ){
				E += H*C/((Photon *)current)->wavelength;
			}

		}

		K += 1/2.0 * current->mass * pow(norm(current->velocity), 2);
	}

	//std::cout << "U_g:" << U_g << std::endl;
	std::cout << "U_e:" << U_e << std::endl;
	std::cout << "K:" << K << std::endl;
	std::cout << "E:" << E << std::endl;

	long double total_energy = U_g + U_e + K + E;
	std::cout << "total_energy: " << total_energy << std::endl;
}

int openglmain(int argc, const char * argv[]){
	glfwSetErrorCallback(onError);

    GLFWwindow* window = createWindow();
    if (!window)
        return 0;

    glfwMakeContextCurrent(window);

    scene.init();

    int windowWidth = 0;
    int windowHeight = 0;

    glfwGetFramebufferSize(window, &windowWidth, &windowHeight);
    onFramebufferResize(window, windowWidth, windowHeight);
    glfwSetFramebufferSizeCallback(window, &onFramebufferResize);

    glfwSetKeyCallback(window, key_callback);

    while (!glfwWindowShouldClose(window))
    {
		if( !paused ){
			compute_new_values();
		}
		check_conditions();
		check_energy();

		compute_user_input();

		/*
		 std::cout << all_particles[0]->position[0]
			 << " " << all_particles[0]->position[1]
			 << " " << all_particles[0]->position[2]
			 << std::endl;
			 */

        scene.draw();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();

    return 0;
}

int main(){

	long double r_not = 5.29 * pow(10,-11);
	long double v_not = 2.18805743462617 * pow(10,6);


	all_particles.push_back(new Proton());
	all_particles.push_back(new Proton());
	all_particles.push_back(new Electron());
	all_particles.push_back(new Electron());

	all_particles[0]->position = std::vector<long double>{0.1*r_not, 0, 0};
	all_particles[0]->velocity = std::vector<long double>{0, 0, 0};

	all_particles[1]->position = std::vector<long double>{-0.1*r_not, 0, 0};
	all_particles[1]->velocity = std::vector<long double>{0, 0, 0};

	all_particles[2]->position = std::vector<long double>{r_not, 0, 0};
	all_particles[2]->velocity = std::vector<long double>{0, -v_not, 0};

	all_particles[3]->position = std::vector<long double>{-r_not, 0, 0};
	all_particles[3]->velocity = std::vector<long double>{0, v_not, 0};

	for (uint i = 0; i < all_particles.size(); i++ ){ // calculate new values

		Particle* current = all_particles[i];

		for (uint j = 0; j < all_particles.size(); j++ ){ // calculate new acceleration values
			if (i==j){ continue; }

			Particle* other = all_particles[j];

			// calculate
			// all_particles[i]->next_position

			long double dist = sqrt(
				  pow(current->position[0]-other->position[0], 2)
				+ pow(current->position[1]-other->position[1], 2)
				+ pow(current->position[2]-other->position[2], 2) );

			long double cl_scalar = - k_e * current->charge * other->charge / dist / dist / current->mass / dist;

			std::vector<long double> delta_vector = std::vector<long double>{
				other->position[0]-current->position[0],
				other->position[1]-current->position[1],
				other->position[2]-current->position[2]};

			std::vector<long double> new_acceleration = std::vector<long double>{
				cl_scalar*delta_vector[0],
				cl_scalar*delta_vector[1],
				cl_scalar*delta_vector[2]};

			current->last_acceleration = current->acceleration;
			current->acceleration = new_acceleration;

			/*
			std::cout << all_particles[i]->position[0] << " "
				<< all_particles[i]->position[1] << " "
				<< all_particles[i]->position[2] << " "
			 << std::endl;
			 */
		}
	}

	// std::cout << delta_t << std::endl;

	//return 0;

	openglmain(0, NULL);

	return 0;
}

