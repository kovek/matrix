#include <list>
#include <cmath>

struct Particle {
	public:
		long double mass;
		long double mass_ev;
		long double charge;
		double radius;
		double magnetic_moment;
		bool spin;
		long double wavelength = 0;

		mutable std::vector<long double> position = std::vector<long double>{0.0, 0.0, 0.0};
		mutable std::vector<long double> velocity = std::vector<long double>{0.0, 0.0, 0.0};
		mutable std::vector<long double> acceleration = std::vector<long double>{0.0, 0.0, 0.0};

		std::list< std::vector<long double> > position_back_log;
		std::vector<long double> last_velocity;
		std::vector<long double> last_acceleration = std::vector<long double>{0,0,0};
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

struct Photon: public Particle {
	public:
		Photon(){
			mass_ev = 0;
			mass = 0;
			charge = 0;
			radius = 0;
			magnetic_moment = 0;
			spin = 0;
			wavelength = 1;
		}
};
