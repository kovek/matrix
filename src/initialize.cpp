#include <mpreal.h>
#include <cmath>
#include <vector>

#include "particles.h"

mpfr::mpreal r_not = 2.29 * pow(10,-12);
mpfr::mpreal v_not = 2.18805743462617 * pow(10,4);
mpfr::mpreal new_dist = 0.29 * pow(10,-13);

void initialize(std::vector<Particle *> & all_particles){
	all_particles.push_back(new Proton());
	all_particles.push_back(new Electron());

	all_particles[0]->position = std::vector<mpfr::mpreal>{0, 0, 0};
	all_particles[0]->velocity = std::vector<mpfr::mpreal>{0, 0, 0};

	all_particles[1]->position = std::vector<mpfr::mpreal>{r_not, r_not, r_not};
	all_particles[1]->velocity = std::vector<mpfr::mpreal>{0, 0, 0};
}
