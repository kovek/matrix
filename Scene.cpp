/* changes:
 * disabled culling
 * added but did not finish perspective stuff
 * added object loader and trying to draw object...
 */

#include "Scene.h"
#include <iostream>
#include "particles.h"
#include "debug.h"
#include "glhelpers.h"
#include "GLFW/glfw3.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// these needed for object loader...---------------
#define GLM_FORCE_RADIANS
#define GLM_GTX_norm
#include <glm/glm/glm.hpp>
#include <glm/glm/gtx/vec1.hpp>
#include <glm/glm/gtc/matrix_transform.hpp>
#include <glm/glm/gtc/type_ptr.hpp>
// end of the imports needed for object loader...

extern std::vector<Particle *> all_particles;
extern std::vector<int> keys_pressed;
extern glm::mat4 move(glm::mat4 mat, float x, float y, float z);
extern glm::mat4 rotate(glm::mat4 mat, float x, float y, float z);
extern glm::mat4 scale(glm::mat4 mat, float x, float y, float z);
extern void passDataToArray(glm::mat4 matrix, float* my_array);
extern float identityMatrix[16];


const float vertexData[] = {
    -0.05f,  0.0f, -0.0f, 1.0f,
    0.0f, 0.05f, -0.0f, 1.0f,
	0.05f,  0.0f, -0.0f, 1.0f,

    -0.05f,  0.0f, -0.0f, 1.0f,
    0.05f, 0.0f, -0.0f, 1.0f,
	0.0f,  -0.050f, -0.0f, 1.0f,

    1.0f, 1.0f, 0.0f, 1.0f,
	1.0f, 0.0f, 1.0f, 1.0f,
	0.0f, 1.0f, 0.0f, 1.0f,

    1.0f, 1.0f, 0.0f, 1.0f,
	1.0f, 0.0f, 1.0f, 1.0f,
	0.0f, 1.0f, 0.0f, 1.0f
};

const float particle[] = {
    -0.25f,  0.0f, 0.0f, 1.0f,
    0.0f, 0.25f, 0.0f, 1.0f,
	0.25f,  0.0f, 0.00f, 1.0f,

    1.0f, 0.0f, 0.0f, 1.0f,
    1.0f, 0.0f, 0.0f, 1.0f,
    1.0f, 0.0f, 0.0f, 1.0f,
};

const float trail[] = {
    -0.25f,  0.0f, -0.0f, 1.0f,
    0.0f, 0.25f, -0.0f, 1.0f,
	0.25f,  0.0f, -0.0f, 1.0f,

    1.0f, 0.0f, 0.0f, 1.0f,
    0.0f, 1.0f, 0.0f, 1.0f,
    0.0f, 0.0f, 1.0f, 1.0f,
};

const float ivector[] = {
	0.0f, 0.0f, 0.0f, 1.0f,
	0.0f, 1.0f, 0.0f, 1.0f,
	1.0f, 0.0f, 0.0f, 1.0f,

	1.0f, 0.0f, 0.0f, 1.0f,
	0.0f, 1.0f, 0.0f, 1.0f,
	0.0f, 0.0f, 1.0f, 1.0f,
};

const GLushort traili[] = {
	0, 1, 2,

	5, 6, 7
};

GLuint offsetUniform;
GLuint objMatrixUniform;
GLuint perspectiveMatrixUniform;
GLuint viewMatrixUniform;

float perspectiveMatrix[16];
float viewMatrix[16];
const float frustumScale = 1.0f;
const float pi = atan(1)*4;



/*
 * For MOVE:
 *
 * x: right
 * y: up
 * z: out of screen
 *
 * For ROTATE:
 *
 * x: up
 * y: left
 * z: out of screen
 */

glm::mat4 getVectorFromTwoPoints(glm::vec3 a, glm::vec3 b){
	glm::mat4 out = glm::mat4(1.0f);

	glm::vec3 dir = glm::vec3(b[0]-a[0], b[1]-a[1], b[2]-a[2]);
	glm::vec3 ivec = glm::vec3(1.0f, 0, 0);

	float norm = sqrt( pow(dir[0], 2) + pow(dir[1], 2) + pow(dir[2], 2) );

	out = scale(out, 1, 1, 1);

	//out = rotate(out, 0, 0, pi/2);


	/*
	float thetay = atan(dir[2]/dir[0]); // x/z
	float thetaz = atan2(dir[1], dir[2]); // y/z
	*/

	float thetay = atan2(dir[2], dir[1]); // x/z
	float thetaz = acos(dir[0]/norm); // y/z


	out = rotate(out, 0, 0, thetaz);
	out = rotate(out, 0, thetay, 0);

	//out = rotate(out, pi/2, 0, -pi/2);
	//out = rotate(out, pi/2, 0, 0);
	//out = rotate(out, pi/2, -pi/2, 0);

	return out;
}

void load_obj(const char* filename, std::vector<glm::vec4> &vertices, std::vector<glm::vec3> &normals, std::vector<GLushort> &elements) {
  std::ifstream in(filename, std::ios::in);
  if (!in) { std::cerr << "Cannot open " << filename << std::endl; exit(1); }

  std::string line;
  while (getline(in, line)) {
    if (line.substr(0,2) == "v ") {
		std::istringstream s(line.substr(2));
      glm::vec4 v; s >> v.x; s >> v.y; s >> v.z; v.w = 1.0f;
      vertices.push_back(v);
    }  else if (line.substr(0,2) == "f ") {
		std::istringstream s(line.substr(2));
      GLushort a,b,c;
      s >> a; s >> b; s >> c;
      a--; b--; c--;
      elements.push_back(a); elements.push_back(b); elements.push_back(c);
    }
    else if (line[0] == '#') { /* ignoring this line */ }
    else { /* ignoring this line */ }
  }

  normals.resize(vertices.size(), glm::vec3(0.0, 0.0, 0.0));
  for (int i = 0; i < elements.size(); i+=3) {
    GLushort ia = elements[i];
    GLushort ib = elements[i+1];
    GLushort ic = elements[i+2];
    glm::vec3 normal = glm::normalize(glm::cross(
      glm::vec3(vertices[ib]) - glm::vec3(vertices[ia]),
      glm::vec3(vertices[ic]) - glm::vec3(vertices[ia])));
    normals[ia] = normals[ib] = normals[ic] = normal;
  }
}

////////////////////////////////////////
GLuint attribute_v_coord;
GLuint attribute_v_normal;

GLuint vbo_mesh_vertices;
GLuint vbo_mesh_normals;

GLuint line_vertices_buffer;

GLuint ibo_mesh_elements;

std::vector<glm::vec4> sphere_vertices;
std::vector<glm::vec3> sphere_normals;
std::vector<GLushort> sphere_elements;

std::vector<glm::vec4> cone_vertices;
std::vector<glm::vec3> cone_normals;
std::vector<GLushort> cone_elements;

GLuint _shaderProgram;
GLuint spherevbo;
GLuint conevbo;
GLuint linevbo;
int sphereSize;
int coneSize;
int lineSize;

void cone() {
    glGenVertexArrays(1, &conevbo);
	glBindVertexArray(conevbo);

	//attribute_v_coord = glGetAttribLocation(_shaderProgram, "position");
	//attribute_v_normal = glGetAttribLocation(_shaderProgram, "normal");

	glEnableVertexAttribArray(0); // position
	glEnableVertexAttribArray(1);
	//glEnableVertexAttribArray(attribute_v_normal);
	// Describe our vertices array to OpenGL (it can't guess its format automatically)

	glGenBuffers(1, &vbo_mesh_vertices);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_mesh_vertices);
	glBufferData(GL_ARRAY_BUFFER, cone_vertices.size() * sizeof(glm::vec4) * 2, &cone_vertices[0], GL_STATIC_DRAW);
	glVertexAttribPointer(
		0,
		4,                  // number of elements per vertex, here (x,y,z,w)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		0                   // offset of first element
	);

	/*
	glGenBuffers(1, &vbo_mesh_normals);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_mesh_normals);
	glBufferData(GL_ARRAY_BUFFER, sizeof(cone_normals), &cone_normals[0], GL_STATIC_DRAW);
	glVertexAttribPointer(
		attribute_v_normal, // attribute
		3,                  // number of elements per vertex, here (x,y,z)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		0                   // offset of first element
	);
	*/

	glVertexAttribPointer(
		1,
		4,                  // number of elements per vertex, here (x,y,z,w)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		0                   // offset of first element
	);

	glGenBuffers(1, &ibo_mesh_elements);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_mesh_elements);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, cone_elements.size() * sizeof(GLushort), &cone_elements[0], GL_STATIC_DRAW);

   	glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &coneSize);

	glBindVertexArray(0);

	//glDrawArrays(GL_TRIANGLES, 0, sizeof(cone_vertices)*sizeof(glm::vec4) );
}

void sphere() {

    glGenVertexArrays(1, &spherevbo);
	glBindVertexArray(spherevbo);

	//attribute_v_coord = glGetAttribLocation(_shaderProgram, "position");
	//attribute_v_normal = glGetAttribLocation(_shaderProgram, "normal");

	glEnableVertexAttribArray(0); // position
	glEnableVertexAttribArray(1);
	//glEnableVertexAttribArray(attribute_v_normal);
	// Describe our vertices array to OpenGL (it can't guess its format automatically)

	glGenBuffers(1, &vbo_mesh_vertices);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_mesh_vertices);
	glBufferData(GL_ARRAY_BUFFER, sphere_vertices.size() * sizeof(glm::vec4) * 2, &sphere_vertices[0], GL_STATIC_DRAW);
	glVertexAttribPointer(
		0,
		4,                  // number of elements per vertex, here (x,y,z,w)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		0                   // offset of first element
	);

	/*
	glGenBuffers(1, &vbo_mesh_normals);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_mesh_normals);
	glBufferData(GL_ARRAY_BUFFER, sizeof(sphere_normals), &sphere_normals[0], GL_STATIC_DRAW);
	glVertexAttribPointer(
		attribute_v_normal, // attribute
		3,                  // number of elements per vertex, here (x,y,z)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		0                   // offset of first element
	);
	*/

	glVertexAttribPointer(
		1,
		4,                  // number of elements per vertex, here (x,y,z,w)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		0                   // offset of first element
	);

	glGenBuffers(1, &ibo_mesh_elements);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_mesh_elements);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sphere_elements.size() * sizeof(GLushort), &sphere_elements[0], GL_STATIC_DRAW);

   	glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &sphereSize);

	glBindVertexArray(0);

	//glDrawArrays(GL_TRIANGLES, 0, sizeof(sphere_vertices)*sizeof(glm::vec4) );
}
///////////////////////////////////////////

void line() {

    glGenVertexArrays(1, &linevbo);
	glBindVertexArray(linevbo);

	glEnableVertexAttribArray(0); // position
	glEnableVertexAttribArray(1); // color
	//glEnableVertexAttribArray(attribute_v_normal);
	// Describe our vertices array to OpenGL (it can't guess its format automatically)

	glGenBuffers(1, &line_vertices_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, line_vertices_buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(ivector)*sizeof(float), ivector, GL_STATIC_DRAW);
	glVertexAttribPointer(
		0,
		4,                  // number of elements per vertex, here (x,y,z,w)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		0                   // offset of first element
	);

	size_t size = sizeof(ivector)/2;

	glVertexAttribPointer(
		1,
		4,                  // number of elements per vertex, here (x,y,z,w)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		(void*) size                    // offset of first element
	);

	glBindVertexArray(0);
}


Scene::Scene()
{
}

void Scene::init()
{
	{ // Load sphere object
		load_obj("sphere.obj", sphere_vertices, sphere_normals, sphere_elements);
		load_obj("cone.obj", cone_vertices, cone_normals, cone_elements);
	}

	sphere();
	cone();
	line();

    _shaderProgram = createShaderProgramWithFilenames("MatrixPerspective.vert", "StandardColors.frag");
    glUseProgram(_shaderProgram);

    // Initialize Vertex Buffer
    glGenBuffers(1, &_vertexBufferObject);

	glBindBuffer(GL_ARRAY_BUFFER, _vertexBufferObject);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Vertex array object
    glGenVertexArrays(1, &_vertexArrayObject);
	glBindVertexArray(_vertexArrayObject);

	/*
    // Enable cull facing
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CW);
	*/

	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	glDepthFunc(GL_LEQUAL);
	glDepthRange(0.0f, 1.0f);

    // Uniforms
	offsetUniform = glGetUniformLocation(_shaderProgram, "offset");
	objMatrixUniform = glGetUniformLocation(_shaderProgram, "objMatrix");

	perspectiveMatrixUniform = glGetUniformLocation(_shaderProgram, "perspectiveMatrix");
	viewMatrixUniform = glGetUniformLocation(_shaderProgram, "viewMatrix");

	float zNear = 0.5f; float zFar = 30.0f;

	memset(perspectiveMatrix, 0, sizeof(float) * 16);
	memset(viewMatrix, 0, sizeof(float) * 16);

	viewMatrix[0] = 1;
	viewMatrix[5] = 1;
	viewMatrix[10] = 1;
	viewMatrix[15] = 1;
	viewMatrix[14] = 1;

	perspectiveMatrix[0] = frustumScale;
	perspectiveMatrix[5] = frustumScale;
	perspectiveMatrix[11] = 1.0f;
	perspectiveMatrix[10] = (zFar + zNear) / (zNear - zFar);
	perspectiveMatrix[14] = (2 * zFar * zNear) / (zNear - zFar);
	perspectiveMatrix[11] = -1.0f;

	glUniformMatrix4fv(perspectiveMatrixUniform, 1, GL_FALSE, perspectiveMatrix);
    glUniformMatrix4fv(viewMatrixUniform, 1, GL_FALSE, viewMatrix);
}

Scene::~Scene()
{
    glDeleteProgram(_shaderProgram);
    glDeleteBuffers(1, &_vertexBufferObject);
}

void Scene::updateCamera() {
    glUniformMatrix4fv(viewMatrixUniform, 1, GL_FALSE, viewMatrix);

	this->draw();
}

void Scene::reshape(int width, int height)
{
    perspectiveMatrix[0] = frustumScale / (width / (float)height);
    perspectiveMatrix[5] = frustumScale;

    glUniformMatrix4fv(perspectiveMatrixUniform, 1, GL_FALSE, perspectiveMatrix);

    this->draw();
}


GLuint trailBufferObject;
GLuint tvaoObject;
GLuint trailibuf;

uint trail_counter;

float t;

float norm(glm::vec3 vec){
	return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
}

int z = 0;

void Scene::draw()
{
    glClearColor(0.2f, 0.0f, 0.0f, 0.0f);
	glClearDepth(1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	t += 0.1f;


	glm::mat4 transform;
	float length_of_vector;
	{
		float floats_array[16] = {0};
		glBindVertexArray(spherevbo);
		glUniform3f(offsetUniform, 0, 0, 0);

		glm::mat4 mod = glm::mat4(1.0f);
		mod = scale(mod, 1/7.0f, 1/7.0f, 1/7.0f);
		passDataToArray(mod, floats_array);

		glUniformMatrix4fv(objMatrixUniform, 1, GL_FALSE, floats_array);

		glDrawElements(GL_TRIANGLES, sphereSize, GL_UNSIGNED_SHORT, 0);
	}

	{
		float floats_array[16] = {0};
		glm::mat4 mod = glm::mat4(1.0f);
		mod = scale(mod, 1/7.0f, 1/7.0f, 1/7.0f);


		std::vector<long double> acc = all_particles[0]->acceleration;
		std::vector<long double> pos = all_particles[0]->position;

		length_of_vector = norm(glm::vec3(acc[0], acc[1], acc[2]));

		mod = move(mod, 0, length_of_vector/9/pow(10,8), 0);

		mod = rotate(mod, pi/2, 0, 0);




		mod = rotate(mod, 0, 0, -pi/2);

		//mod = move(mod, 0, 0, 3);

		transform = getVectorFromTwoPoints( glm::vec3{0,0,0}, glm::vec3{acc[0], acc[1], acc[2]} );
		mod = transform*mod;

		mod = move(mod, pos[0], pos[1], pos[2]);

		passDataToArray(mod, floats_array);

		glUniformMatrix4fv(objMatrixUniform, 1, GL_FALSE, floats_array);




		glBindVertexArray(conevbo);
		//glUniform3f(offsetUniform, all_particles[0]->position[0], all_particles[0]->position[1], all_particles[0]->position[2]);
		glUniform3f(offsetUniform, 0, 0, 0);
		//glDrawElements(GL_TRIANGLES, coneSize, GL_UNSIGNED_SHORT, 0);

		glUniformMatrix4fv(objMatrixUniform, 1, GL_FALSE, identityMatrix);
	}

	{
		float floats_array[16] = {0};
		glm::mat4 new_mod = glm::mat4(1);
		new_mod = scale(new_mod, 0, length_of_vector/9/pow(10,8), 0);
		new_mod = rotate(new_mod, pi/2, 0, -pi/2);
		std::vector<long double> pos = all_particles[0]->position;
		new_mod = transform*new_mod;
		new_mod = move(new_mod, pos[0], pos[1], pos[2]);

		passDataToArray(new_mod, floats_array);
		glUniformMatrix4fv(objMatrixUniform, 1, GL_FALSE, floats_array);

		glBindVertexArray(linevbo);
		glUniform3f(offsetUniform, 0.0f, 0.0f, 0.0f);
		glDrawArrays(GL_LINES, 0, 2);

		glUniformMatrix4fv(objMatrixUniform, 1, GL_FALSE, identityMatrix);
	}


	for(uint t = 0; t < all_particles.size(); t++) {
		Particle * current = all_particles[t];

		glBindVertexArray(_vertexArrayObject);

		glGenBuffers(1, &trailBufferObject);
		glBindBuffer(GL_ARRAY_BUFFER, trailBufferObject);
		glBufferData(GL_ARRAY_BUFFER, sizeof(trail), trail, GL_STATIC_DRAW);

		glGenBuffers(1, &trailibuf);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, trailibuf);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(traili), traili, GL_STATIC_DRAW);

		/*
		glGenVertexArrays(1, &tvaoObject);
		glBindVertexArray(tvaoObject);
		*/

		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);

		size_t colorpos = sizeof(trail)/sizeof(float)/2;

		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, (void*) colorpos);

		std::list< std::vector<long double> >* log = &(current->position_back_log);
		std::list< std::vector<long double> >::iterator cur = log->begin();
		std::list< std::vector<long double> >::iterator end = log->end();

		uint skip_every_nth = 2;
		uint i;
		float j = log->size();
		bool stop_iterating = false;
		trail_counter++;
		trail_counter = trail_counter % skip_every_nth;

		for(uint k = 0; k < 18; k++){
			cur++;
		}


		while (cur != end) { // draw particles
			j -= (float)skip_every_nth;

			//glUniform3f(offsetUniform, (*cur)[0]*pow(10,10), (*cur)[1]*pow(10,10), (*cur)[2]*pow(10,10));
			glUniform3f(offsetUniform, 0, 0, 0);
			//glDrawArrays(GL_TRIANGLES, 0, sizeof(trail)/sizeof(float)/2 );

			float floats_array[16] = {0};
			glm::mat4 mod = glm::mat4(1.0f);
			mod = scale(mod, 5/j, 5/j, 5/j);
			mod = move(mod, (*cur)[0]*pow(10,10), (*cur)[1]*pow(10,10), (*cur)[2]*pow(10,10));
			passDataToArray(mod, floats_array);

			glUniformMatrix4fv(objMatrixUniform, 1, GL_FALSE, floats_array);

			int size;
			glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);

			glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_SHORT, 0);


			i = 0;
			while(i!=skip_every_nth)
			{
				cur++;
				i++;
				if(cur == end){ stop_iterating = true; }
			}
			if (stop_iterating){ break; }
		}
	}



	for (uint k = 0; k < all_particles.size(); k++){
		Particle* current = all_particles[k];

		glUniform3f(offsetUniform, (float)current->position[0], (float)current->position[1], (float)current->position[2]);

		if(z == 70){
			//throw std::invalid_argument("test");
		}
		z++;

		size_t colorData = sizeof(vertexData) / 2;
		glBindBuffer(GL_ARRAY_BUFFER, _vertexBufferObject);
		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, (void*)colorData);


		glDrawArrays(GL_TRIANGLES, 0, sizeof(vertexData)/sizeof(float)/2/4 );

		glBindVertexArray(spherevbo);


		float floats_array[16] = {0};
		glm::mat4 mod = glm::mat4(1.0f);
		mod = scale(mod, 1/7.0f, 1/7.0f, 1/7.0f);
		passDataToArray(mod, floats_array);

		glUniformMatrix4fv(objMatrixUniform, 1, GL_FALSE, floats_array);

		glUniform3f(offsetUniform, (float)current->position[0]*7*pow(10,10), (float)current->position[1]*7*pow(10,10), (float)current->position[2]*7*pow(10,10));
		glDrawElements(GL_TRIANGLES, sphereSize, GL_UNSIGNED_SHORT, 0);
		glBindVertexArray(0);
	}
}

void Scene::keyStateChanged(int key, int action)
{
	if (action == GLFW_PRESS) {
		std::vector<int>::iterator p = std::find(keys_pressed.begin(), keys_pressed.end(), key);
		if (p == keys_pressed.end() ){
			keys_pressed.push_back(key);
		}
	}else if(action == GLFW_RELEASE) {
		std::vector<int>::iterator p = std::find(keys_pressed.begin(), keys_pressed.end(), key);
		keys_pressed.erase(p);
	}
}




