//simulates atoms in a 2D box that interact via Lennard-Jones potential, using Verlet integration algorithm
//compilation: g++ -O3 -std=c++11 -o verlet verlet.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class Particle {
    public:
        double x;
        double y;
        double vx;
        double vy;
        double ax = 0.0;
        double ay = 0.0;
        double prev_ax = 0.0;
        double prev_ay = 0.0;
        Particle(double x_, double y_, double vx_, double vy_) {
            x = x_; y = y_; vx = vx_; vy = vy_;
        }
        //updates particle position, with timestep h
        void update(double h) {
            x += vx*h + 0.5*ax*pow(h, 2);
            y += vy*h + 0.5*ay*pow(h, 2);
        }
        //returns magnitude of velocity squared
        double velSquared() {
            return (pow(vx, 2) + pow(vy, 2));
        }
        ~Particle() {}
};

class Grid {
    public:
        double density;
        int N;    //number of particles
        double L; //side length
        double h;
        vector<Particle> particles;
        //constructor with specifications as requested in "guia 6"
        Grid(int N_c, double h_, double density_) {
            density = density_;
            N = N_c*N_c;
            h = h_;
            L = sqrt(N/density);
            double a = (double)L/(N_c + 1);
            double v_0[2] {-1.1, 1.1};
            for (int n = 1; n <= N_c; n++) {
                for (int m = 1; m <= N_c; m++) {
                    particles.push_back( Particle(n*a, m*a, v_0[rand() % 2], 0.0) );
                }
            }
            //initial accelerations
            for (int i = 0; i < N; i++) {
                double ax = 0;
                double ay = 0;
                for (int j = 0; j < N; j++) {
                    if (i != j) F(particles[i], particles[j], ax, ay);
                }
                particles[i].ax = ax;
                particles[i].ay = ay;
            }
            cout << "\nGrid has been initialized. Info:\n";
            cout << " - Number of particles N = " << N_c << "x" << N_c << " = " << N << endl;
            cout << " - Density           rho = " << density << endl;
            cout << " - Side of box length  L = " << L << endl;
            cout << " - Timestep            h = " << h << endl;
            cout << "\nStarting simulation...\n\n";
        }

        void updateParticles() {
            //update particles positions
            for (int i = 0; i < N; i++) {
                //update particle[i]'s position, with timestep h
                particles[i].update(h);
                //check if out of bounds
                if ( particles[i].x < 0 ) {
                    particles[i].vx *= -1;
                    particles[i].x = -particles[i].x;
                }
                if ( particles[i].x > L ) {
                    particles[i].vx *= -1;
                    particles[i].x = 2*L - particles[i].x;
                }
                if ( particles[i].y < 0 ) {
                    particles[i].vy *= -1;
                    particles[i].y = -particles[i].y;
                }
                if ( particles[i].y > L ) {
                    particles[i].vy *= -1;
                    particles[i].y = 2*L - particles[i].y;
                }
            }
            //update accelerations and velocities
            for (int i = 0; i < N; i++) {
                //set prev_ax, prev_ay of particle to newer, updated values
                particles[i].prev_ax = particles[i].ax;
                particles[i].prev_ay = particles[i].ay;
                //temporary values to store the cumulative acceleration inflicted
                //on particle i by all other particles j != i
                double ax = 0;
                double ay = 0;
                for (int j = 0; j < N; j++) {
                    //as we pass ax and ay to F(), the latter adds the acceleration inflicted
                    //on particle i by particle j to ax and ay
                    if (i != j) F(particles[i], particles[j], ax, ay);
                }
                //we now set the particles acceleration components
                particles[i].ax = ax;
                particles[i].ay = ay;
                //finally, we update velocities
                particles[i].vx += 0.5*(particles[i].ax + particles[i].prev_ax)*h;
                particles[i].vy += 0.5*(particles[i].ay + particles[i].prev_ay)*h;
            }
        }

        //returns cinetic energy
        double E_cin() {
            double E_k = 0.0;
            for (int i = 0; i < N; i++) E_k += particles[i].velSquared();
            return 0.5*E_k;
        }

        //return potential energy
        double E_pot() {
            double E_p = 0.0;
            for (int i = 0; i < N; i++) {
                for (int j = i+1; j < N; j++) {
                    E_p += U(particles[i], particles[j]);
                }
            }
            return E_p;
        }

        //returns total energy
        double E_tot() {
            return E_cin()+E_pot();
        }

        ~Grid() {}

    private:
        //calculate potential between two particles
        double U(Particle a, Particle b) {
            double r = pow(a.x - b.x, 2) + pow(a.y - b.y, 2);
            double u = 4 * (pow(r, -6) - pow(r, -3));
            return u;
        }
        //calculate acceleration inflicted on particle a by particle b
        //arguments are the two particles, and references to two doubles
        //that represent acceleration_x and acceleration_y
        //the function adds the components calculated to ax and ay
        void F(Particle a, Particle b, double &ax, double &ay) {
            double r = pow(a.x - b.x, 2) + pow(a.y - b.y, 2);
            double f = 24 * (2 * pow(r, -7) - pow(r, -4));
            ax += (a.x - b.x) * f;
            ay += (a.y - b.y) * f;
        }
};

int main(int argc, char* argv[]) {
    
    //some example code that will simulate 900 particles, placed with density 0.3, and timestep h = 0.005
    //the 3 energies will be saved on each iteration to file "energies_out.csv"
    //a total of 2000 iterations are saved

    char const* filename = "energies_out.csv";
    fstream file;
    file.open(filename, ios::out | ios::trunc);
    if (file.is_open()) {
        cout << "Will save to file: " << filename << endl;
    }
    else {
        cout << "Couldn't open file: " << filename << ". Exiting...\n";
        return -1;
    }

    int N_c = 30;
    double h = 0.005;
    double density = 0.3;

    //create a Grid object where simulation will take place
    Grid g = Grid(N_c, h, density);
    
    //save initial energy values. fomart (delimiter=" "): E_TOTAL E_CINETIC E_POTENTIAL
    file << g.E_tot() << " " << g.E_cin() << " " << g.E_pot() << endl;
    for (int i = 1; i <= 2000; i++) {
        g.updateParticles();
        //print progress of simulation
        cout << "\rIteration: " << i << flush;
        //save energies to file
        file << g.E_tot() << " " << g.E_cin() << " " << g.E_pot() << endl;
    }
    cout << "\nDone...\n";
    file.close();
    return 0;
}