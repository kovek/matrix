#include <mpreal.h>
#include <cmath>
#include <vector>

#include "particles.h"

mpfr::mpreal r_not = 5.29 * pow(10,-11);
mpfr::mpreal v_not = 2.18805743462617 * pow(10,6);
mpfr::mpreal new_dist = 0.29 * pow(10,-13);

void initialize(std::vector<Particle *> & all_particles){
	all_particles.push_back(new Proton());
	all_particles.push_back(new Proton());
	all_particles.push_back(new Electron());
	all_particles.push_back(new Electron());

	all_particles[0]->position = std::vector<mpfr::mpreal>{0, 0.05*r_not, 0};
	all_particles[0]->velocity = std::vector<mpfr::mpreal>{0, 0, v_not*0.00131};
	//all_particles[0]->charge = 1.602176565*pow(10,-19)*1;

	all_particles[1]->position = std::vector<mpfr::mpreal>{0, 0, 3*r_not};
	all_particles[1]->velocity = std::vector<mpfr::mpreal>{0, 0, 0.00131*v_not};

	all_particles[2]->position = std::vector<mpfr::mpreal>{-r_not, 0, 0};
	all_particles[2]->velocity = std::vector<mpfr::mpreal>{0, -0.6*v_not, 0};

	all_particles[3]->position = std::vector<mpfr::mpreal>{r_not, 0, 3*r_not};
	all_particles[3]->velocity = std::vector<mpfr::mpreal>{0, 0.6*v_not, 0};
}
