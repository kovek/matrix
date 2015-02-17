#include <stdio.h>
#include <mpreal.h>
#include <list>
#include <cmath>
#include <string>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_free.hpp>


namespace boost {
namespace serialization {
	template<class Archive>
	void save(Archive & ar, mpfr::mpreal const & m, const unsigned int version) {
		std::string foo = m.toString();
		ar & foo;
		//std::cout << ">>>" << version << std::endl;
	}

	template<class Archive>
	void load(Archive & ar, mpfr::mpreal & m, const unsigned int version){
		std::string s;
		ar & s;
		m = mpfr::mpreal(s.c_str());
	}

	template<class Archive>
	void serialize(Archive & ar, mpfr::mpreal & m, const unsigned int version){
		split_free(ar, m, version);
	}
}
}

struct Particle {
    private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version) {
			ar & position;
			ar & velocity;
			ar & acceleration;
		}


	public:
		mpfr::mpreal mass;
		mpfr::mpreal mass_ev;
		mpfr::mpreal charge;
		double radius;
		double magnetic_moment;
		bool spin;
		mpfr::mpreal energy = 0;
		std::vector<mpfr::mpreal> kinetic_energy = std::vector<mpfr::mpreal>{0, 0, 0};

		mutable std::vector<mpfr::mpreal> position = std::vector<mpfr::mpreal>{0.0, 0.0, 0.0};
		mutable std::vector<mpfr::mpreal> velocity = std::vector<mpfr::mpreal>{0.0, 0.0, 0.0};
		mutable std::vector<mpfr::mpreal> acceleration = std::vector<mpfr::mpreal>{0.0, 0.0, 0.0};

		std::list< std::vector<mpfr::mpreal> > position_back_log;
		std::vector<mpfr::mpreal> last_velocity;
		std::vector<mpfr::mpreal> last_acceleration = std::vector<mpfr::mpreal>{0,0,0};
};

struct Electron: public Particle {
	public:
		Electron(){
			mass_ev = 0.510998928l; // in MeV/c^2
			mass = 9.11l*pow(10,-31); // in kg
			charge = -1.602176565l*pow(10,-19); // in Coulombs
			radius = 0; // Doesn't matter
			magnetic_moment = -1.00115965218076; // in /mu/B
			spin = true;
		}
};

struct Proton: public Particle {
	public:
		Proton(){
			mass_ev = 938.272046; // in MeV/c^2
			mass = 1.672621777 * pow(10,-27); // in kg
			charge = 1.602176565*pow(10,-19); // in Coulombs
			radius = 0; // Doesn't matter
			magnetic_moment = 0.001521032210; // in /mu/B
			spin = true;
		}
};

struct Neutron: public Particle {
	public:
		Neutron(){
			mass_ev = 938.272046; // in MeV/c^2
			mass = 1.672621777 * pow(10,-27); // in kg
			charge = 0; // in Coulombs
			radius = 0; // Doesn't matter
			magnetic_moment = 0.001521032210; // in /mu/B
			spin = true;
		}
};

struct Photon: public Particle {
	public:
		Photon(){
			mass_ev = 0;
			mass = 0;
			charge = 0;
			radius = 0;
			magnetic_moment = 0;
			spin = 0;
			energy = 1;
		}
};
