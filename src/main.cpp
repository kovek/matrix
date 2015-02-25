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
#define null_callback NULL
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
#include "initialize.h"
#include <string>
#include <thread>
#include <time.h>
#include <regex>
#include <boost/algorithm/string.hpp>
#include "md5.h"
#include <sqlite3.h>

#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(Proton, "Proton")
BOOST_CLASS_EXPORT_GUID(Electron, "Electron")
BOOST_CLASS_EXPORT_GUID(Neutron, "Neutron")

std::string times;
bool followup_exists = false;
std::string load_filename;
std::string load_iteration;
void load_everything(std::string filename, std::string iteration);
std::string current_commit;


class x00 {
public:
    x00( int, const char ** );

private:
    // Class data
    static const int NSIZE;
};



const mpfr::mpreal k_e = 8.9875517873681764*pow(10,9);
const mpfr::mpreal delta_t = 1.52l*pow(10, -15)/22500l;
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

void compute_new_values() {
	std::cout << load_iteration << std::endl;
	std::cout << "acc_x: " << all_particles[0]->last_acceleration[0] << std::endl;
	if (followup_exists){
		load_everything(load_filename, load_iteration);
		std::cout << "NOT COMPUTING SHIT!" << std::endl;
		return;
	}

	std::cout << "Number of particles: " << all_particles.size() << std::endl;

	// Compute new position
	for (uint i = 0; i < all_particles.size(); i++ ){
		Particle* current = all_particles[i];

		current->position[0] += current->velocity[0]*delta_t+0.5*current->acceleration[0]*delta_t*delta_t;
		current->position[1] += current->velocity[1]*delta_t+0.5*current->acceleration[1]*delta_t*delta_t;
		current->position[2] += current->velocity[2]*delta_t+0.5*current->acceleration[2]*delta_t*delta_t;
		current->position_back_log.push_back( current->position );
		if (current->position_back_log.size()> 50){ current->position_back_log.pop_front(); }

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

			mpfr::mpreal awe_scalar = -(pow(10, 15) * exp(-dist*pow(10,15) ) )/dist - (exp(-dist*pow(10,15) ))/ dist / dist;

			//std::cout << awe_scalar << std::endl;

			/*
			awe_scalar = 0;
			cl_scalar = 0;
			gravity_scalar = 0;
			*/

			std::vector<mpfr::mpreal> delta_vector = std::vector<mpfr::mpreal>{
				other->position[0]-current->position[0],
				other->position[1]-current->position[1],
				other->position[2]-current->position[2]};

			std::vector<mpfr::mpreal> new_force = std::vector<mpfr::mpreal>{
				(awe_scalar + gravity_scalar + cl_scalar) *delta_vector[0],
				(awe_scalar + gravity_scalar + cl_scalar) *delta_vector[1],
				(awe_scalar + gravity_scalar + cl_scalar) *delta_vector[2]};

			sum_of_all_forces[0] += new_force[0];
			sum_of_all_forces[1] += new_force[1];
			sum_of_all_forces[2] += new_force[2];


			std::cout << "SUP!" << std::endl;
			std::cout << all_particles[i]->position[0] << " "
				<< all_particles[i]->position[1] << " "
				<< all_particles[i]->position[2] << " "
			 << std::endl;
		}
		std::vector<mpfr::mpreal> new_acceleration = std::vector<mpfr::mpreal>{
			sum_of_all_forces[0]/(current->mass*my_gamma( current->velocity[0] )),
			sum_of_all_forces[1]/(current->mass*my_gamma( current->velocity[1] )),
			sum_of_all_forces[2]/(current->mass*my_gamma( current->velocity[2] ))};


		current->last_acceleration = current->acceleration;
		current->acceleration = new_acceleration;

		std::cout << "Last acceleration then acceleration then mass: " << std::endl;
		std::cout << current->last_acceleration[0] << std::endl;
		std::cout << current->last_acceleration[1] << std::endl;
		std::cout << current->last_acceleration[2] << std::endl;
		std::cout << current->acceleration[0] << std::endl;
		std::cout << current->acceleration[1] << std::endl;
		std::cout << current->acceleration[2] << std::endl;
		std::cout << "mass" << current->mass_ev << std::endl;
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

		/*
		std::cout <<
			">>>" << std::endl
			<< current->velocity[0]/C << std::endl
			<< current->velocity[1]/C << std::endl
			<< current->velocity[2]/C << std::endl;
		*/


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

			/*
			std::cout
				<< "Gamma factor: " << gamma_factor << std::endl
				<< "Initial Mass: " << mass_initial << std::endl
				<< "Relative Energy : " << relativistic_energy << std::endl
				<< "Delta v : " << delta_v << std::endl
				<< "Delta Energy : " << delta_energy << std::endl
				<< "Energy Final : " << energy_final << std::endl
				<< "Velocity Final : " << velocity_final << std::endl
					<< std::endl;
			*/

			//mpfr::mpreal energy = 1.0f/2.0f * current->mass * C * C * ( gamma_factor - 1);

			current->velocity[j] = velocity_final;
		}

		for (uint j = 0; j < 3; j++){
			current->kinetic_energy[j] = 1.0f/2.0f * current->mass * pow(current->velocity[j], 2) * (my_gamma(current->velocity[j])-1);
		}

		/*
		std::cout << "in sqrt" <<
			1.0f - 1.0f/pow(
				1.0f/sqrt(
					1.0f- pow(
						initial_velocity[0]/C
					,2)
				) + abs(final_velocity[0]-initial_velocity[0])/2/C/C
			, 2) << std::endl;

		std::cout << "in pow in sqrt" <<
				1.0f/sqrt(
					1.0f- pow(
						initial_velocity[0]/C
					,2)
				) + (final_velocity[0]-initial_velocity[0])/2/C/C << std::endl;

		std::cout << "in pow in pow in sqrt" <<
					1.0f- pow(
						initial_velocity[0]/C
					,2) << std::endl;

		std::cout << "ratio" <<
						initial_velocity[0]/C << std::endl;


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
    // pls->parseopts( &argc, argv, PL_PARSE_FULL );

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

				//std::cout << all_particles.size() << std::endl;
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

uint iteration = 0;
void save_everything(){
	iteration++;
	std::string filename = "local/" + times;
	//std::cout << filename << std::endl;
	std::ofstream ofs(filename, std::ofstream::app);

	std::string version = "5";

	{
		boost::archive::text_oarchive oa(ofs);
		ofs << "iteration:" << iteration << std::endl;
		oa << all_particles;
		ofs << std::endl;
		ofs.close();
	}
}

std::string get_blob(){
    std::ostringstream oss("");

	{
		boost::archive::text_oarchive oa{oss};

		// This tells boost that we need to use the proper class
		oa.template register_type<Proton>();
		oa.template register_type<Electron>();
		oa.template register_type<Neutron>();

		oa << all_particles;
	}

    return oss.str();
}

void load_blob(std::string initial_cond){
    std::vector<Particle *> my_particles;
	{
		std::istringstream my_blob_container(initial_cond);
		boost::archive::text_iarchive ia{my_blob_container};

		// This tells boost that we need to use the proper class
		ia.template register_type<Proton>();
		ia.template register_type<Electron>();
		ia.template register_type<Neutron>();

		ia >> my_particles;
	}

    all_particles = my_particles;
}

class gps_position
{
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version) {
			ar & degrees;
			ar & minutes;
			ar & seconds;
		}
		int minutes;
		float seconds;

	public:
		gps_position(){};
		gps_position(int d, int m, float s) :
		degrees(d), minutes(m), seconds(s)
		{}
		int degrees;
};

std::string return_match_string(std::string source, std::string expression){
	std::smatch m;
	std::regex e(expression);
	if(std::regex_search(source, m, e)){
		return (std::string)m[1];
	}else{
		return std::string();
	}
}

std::string get_current_commit(){
	FILE * output = popen("git describe --abbrev=0 --always", "r");

	if(!output){
		// couldn't get output from syscall.
        throw std::runtime_error("Could not get the commit hash");
	}

	char current_commit_hash[41];
	fgets(current_commit_hash, 41, output); // moves output to variable

    return std::string(current_commit_hash);
}

bool correct_commit_hash(std::string commit_hash){

	if (commit_hash == current_commit){
		return true;
	}

	return false;
}

// Fuck sqlite and how it is implemented in c++
std::string sqlite_initial_condition_result;
int callback(void *NotUsed, int argc, char **argv, char **azColName){
	sqlite_initial_condition_result = argv[0];
    return 0;
}

uint parent_id = 0;

std::string get_initial_condition_string(std::string filename, std::string iteration){
	sqlite3 * db;
	char *errmsg = 0;
	int rc;

	rc = sqlite3_open(filename.c_str(), &db);
	if(rc){
        throw std::runtime_error("Couldn't open db");
	}

	std::string sql = "SELECT data FROM timeline WHERE id=" + iteration + ";";

	rc = sqlite3_exec(db, sql.c_str(), callback, 0, &errmsg);
	if (rc != SQLITE_OK){
		std::cout << sql << std::endl;
		std::cout << "filename: " << filename << std::endl;
        throw std::runtime_error("Couldn't call sql");
	}
	sqlite3_close(db);

	parent_id = atoi(iteration.c_str());

	return sqlite_initial_condition_result;
}

bool is_file_exist(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

std::string initial_cond_hash(){
	FILE * pFile;
	char contents[100000];

	pFile = fopen("src/initialize.cpp", "r");
	if (pFile == NULL) { perror("Error opening file"); }
	else {
		fgets(contents, 100000, pFile);
	}

	std::string filecontents = std::string(contents);
	std::string hash = md5(filecontents);
	return std::string(hash);
}

void save_information(){
    system( ("mkdir local/" + current_commit).c_str() );

    std::string initial_conditions_hash = initial_cond_hash();

    std::string filename = "local/" + current_commit + "/" + initial_conditions_hash;

	std::cout << "checking for " << filename << std::endl;
	system( std::string("echo `ls local/93eb4ce2ea1c00fdcb9e6513a1f9a342d71e4b9d`").c_str() );

	sqlite3 * db;
	char *errmsg = 0;
	int rc;

	rc = sqlite3_open(filename.c_str(), &db);
	if(rc){
		throw std::runtime_error("Couldn't open db" );
	}


    std::string blob = get_blob();

    std::string sql = "INSERT INTO timeline (id, parent, data) VALUES (NULL, " +
      std::to_string(parent_id) + ", \"" +
      blob +
    "\");";

    rc = sqlite3_exec(db, sql.c_str(), null_callback, 0, &errmsg);
    if (rc != SQLITE_OK){
		std::cout << sql << std::endl;
        throw std::runtime_error("Coudn't call sql" );
    }

    parent_id = sqlite3_last_insert_rowid(db);

	sqlite3_close(db);
}


int followup_callback(void *NotUsed, int argc, char **argv, char **azColName){
	std::cout << "argv: " << argv[0] << std::endl;
	std::cout << " in it! " << std::endl;
	if(argv[0] != NULL && argv[0][0] != (char)0){
		std::cout << "followup exists" << std::endl;
		followup_exists = true;
	}else{
		std::cout << "followup does not exist" << std::endl;
		followup_exists = false;
	}

	return 0;
}

void check_if_followup_exists(std::string filename, std::string iteration){
	sqlite3 * db;
	char *errmsg = 0;
	int rc;

	rc = sqlite3_open(filename.c_str(), &db);
	if(rc){
		throw std::runtime_error("Couldn't open db" );
	}

	iteration = std::to_string(atoi(iteration.c_str()) + 1);
	load_iteration = iteration;
	parent_id = atoi(load_iteration.c_str());

    std::string sql = "SELECT id FROM timeline WHERE id=" + iteration + ";";

	followup_exists = false;
    rc = sqlite3_exec(db, sql.c_str(), followup_callback, 0, &errmsg);

    if (rc != SQLITE_OK){
		std::cout << sql << std::endl;
        throw std::runtime_error("Coudn't call sql" );
    }

	sqlite3_close(db);
}

void load_everything(std::string filename, std::string iteration){
	std::cout << "in load_everything. iteration: " << iteration << std::endl;

	std::string commit_hash = return_match_string(filename, "^local\\/(.*)\\/.*$");

    // Just making sure we are executing the right data.
	if(!correct_commit_hash(commit_hash)){
        // Forget about this for now...
		// return;
	}

	get_initial_condition_string(filename, iteration);

    load_blob(sqlite_initial_condition_result);

	std::cout << "Checking" << std::endl;
	check_if_followup_exists(filename, iteration);
}


// Convert initial conditions file to hash
// Save that file with "commit_hash/hash" as the filename
// create database with that hash name
// run

/*
 * build/main commit_hash/hash iteration
 * open database. read iteration
 * run while setting parent as the iteration selected.
 */

/*
 * build/main commit_hash/hash
 * open database. read last iteration
 * run.
 */

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

	//std::cout << "U_g:" << U_g << std::endl;

	/*
	std::cout << "U_e:" << U_e << std::endl;
	std::cout << "K:" << K << std::endl;
	std::cout << "E:" << E << std::endl;
	*/

	total_energy = U_g + U_e + K + E;

	//std::cout << "total_energy: " << total_energy << std::endl;
}

clock_t timestamp;
clock_t timestamp2;
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
		//std::cout << "New iteration" << std::endl;

		if( paused ){
			scene.draw();
			compute_user_input();
			continue;
		}


		timestamp = clock();

		//*
		std::thread compute (compute_new_values);

		std::thread check_cond(check_conditions);
		std::thread check_E(check_energy);

		std::thread show_graph(display_graphs);
		compute_user_input();

		//create_photons();

		std::thread save_log;
		if(!followup_exists){
			std::cout << "spawning thread" << std::endl;
			save_log = std::thread(save_information);
		}
		//*/

		/*
		compute_new_values();
		check_conditions();
		check_energy();
		compute_user_input();
		display_graphs();
		save_everything();
		//*/

		/*
		 std::cout << all_particles[0]->position[0]
			 << " " << all_particles[0]->position[1]
			 << " " << all_particles[0]->position[2]
			 << std::endl;
			 */

		//std::thread draw_thread(&Scene::draw, scene);
		scene.draw();

		//*
		compute.join();
		check_cond.join();
		check_E.join();
		show_graph.join();

		if(!followup_exists && save_log.joinable()){
			save_log.join();
		}
		//draw_thread.join();
		//*/

		timestamp2 = clock();
		//std::cout << timestamp2 - timestamp << " " << std::flush;

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();

    return 0;
}

void parse_input(int argc, const char * argv[]){
	for (int i = 1; i < argc; i++) {
        if ( strcmp(argv[i], "-h") == 0 ) {
            load_filename = argv[i + 1];
        } else if ( strcmp(argv[i], "-i") == 0 ) {
            load_iteration = argv[i + 1];
        }
    }
}

int main(int argc, const char * argv[]){
	current_commit = get_current_commit();
	parse_input(argc, argv);

	std::string filename_temp = "local/"+current_commit+"/"+initial_cond_hash();
    if(!is_file_exist(filename_temp.c_str())){
		std::cout << "file does not exist" << std::endl;
		system( ("mkdir local/" + current_commit).c_str() );
		system( ("cp local/template " + filename_temp).c_str() );
    }

	srand (time(NULL));

	time_t timev;
	time(&timev);
	times = std::to_string(timev);

	// PARTEYYYY
	// We are working with `digits` long numbers now hehe.
    const int digits = 500;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));

	std::cout.precision(digits);

    if(!load_filename.empty() && !load_iteration.empty()){
      load_everything(load_filename, load_iteration);
	  std::cout << "LOADING " << std::endl;
    }else{
		try{
			std::cout << "auto load " << std::endl;
			load_filename = "local/"+current_commit+"/"+initial_cond_hash();
			load_iteration = "1";
			load_everything(load_filename, load_iteration);
		} catch(std::exception e){
			initialize(all_particles);
			std::cout << "NOT LOADING " << std::endl;
			std::cout << load_filename << std::endl;
			std::cout << load_iteration << std::endl;
		}
	}



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

	openglmain(argc, argv);

	// Destroy the graphs window setup.
	//delete pls;
    //delete x;

	return 0;
}

