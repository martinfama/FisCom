#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <iterator>
#include <random>
#include <cmath>
#include <ctime>

using namespace std;

//random devices
random_device rd;
mt19937 gen(rd());

//modulus function 
//different from simply remainder, denoted by "%"
//we want, for example, mod(-1, 10) to be 9 and not -1
//so that we can access indices of grid periodically
int mod(int x, int N) {
    return (x % N + N) % N;
}

class Particle {
    public:
        int y;
        int x;
        Particle(int Y, int X) {
            y = Y;
            x = X;
        }
        ~Particle() {}
};

class Grid {

    public:
        int W;
        int H;
        int maxParticles;
        int stp = 0;
        list<Particle> particles;
        vector<vector<int>> aggregate; //points given as (y, x)
        
        Grid(int h, int w, int N) {
            //height and width of grid
            W = w;
            H = h;
            maxParticles = N;
            cout << "Grid size: " << W << " x " << H << " = " << W*H << " cells\n";
            //populate the grid
            uniform_int_distribution<> distribY(0, H-1);
            uniform_int_distribution<> distribX(0, W-1);
            for (int i = 0; i < N; i++) {
                particles.push_back(Particle(distribY(gen), distribX(gen)));
            }
            cout << "Created: " << particles.size() << " particles\n";
            
            //make agreggate
            for (int i = 0; i < H; i++) {
                vector<int> a;
                aggregate.push_back(a);
                for (int j = 0; j < W; j++) {
                    aggregate[i].push_back(0);
                }
            }
            //place seed
            aggregate[H/2][W/2] = 1;
            cout << "Placed seed at: " << H/2 << ", " << W/2 << "\n";
        }

        bool updateParticles() {
            //if there are less than 5% of particles left, we exit (takes too long if not)
            if ( (double)(particles.size())/maxParticles < 0.05) {
                cout << "\n5% of particles left. Updating done...\n";
                return false;
            }
            for (auto i = particles.begin(); i != particles.end(); ) {
                if (fixedNear(i->y, i->x)) {
                    aggregate[i->y][i->x] = 1+stp;
                    i = particles.erase(i);
                }
                else {
                    switch(uniform_int_distribution<int>{0, 3}(gen)) {
                        case 0:
                            i->y -= 1;
                            break;
                        case 1:
                            i->x += 1;
                            break;
                        case 2:
                            i->y += 1;
                            break;
                        case 3:
                            i->x -= 1;
                            break;
                    }
                    i->y = mod(i->y, H);
                    i->x = mod(i->x, W);
                    i++;
                }
            }
            stp++;
            return true;
        }

        //save the agreggate to file
        void saveToFile(string filename = "out_default.csv") {
            fstream file;
            file.open(filename, fstream::out);
            if (!file.is_open()) {
                cout << "\nCouldn't open file: " << filename << "\n";
                return;
            }
            cout << "\nSaving current aggregate to: " << filename << endl;
            for (int y = 0; y < H; y++) {
                for (int x = 0; x < W; x++) {
                    if (aggregate[y][x] >= 1 ) file << x << " " << y << " " << aggregate[y][x] << "\n";
                }
            }
            file.close();
            cout << "Saving done\n";
        }

        ~Grid() {}

    private:

        //checks 8 neighbouring cells
        bool fixedNear(int y, int x) {
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    if (aggregate[mod(y+i, H)][mod(x+j, W)] >= 1) return true;
                }
            }
            return false;
        }
};

int main(int argc, char** argv) {

    //seed random with current time
    gen.seed(time(NULL));

    //check if user has passed a filename to save agreggate to
    //this is passed by passing "-f <filename>"
    //by default, the agreggate is saved to "out_default.csv"
    //be careful not to overwrite a previous agreggate like this
    bool customSaveFile = false;
    string filename;
    for (int i = 1; i < argc; i++) {
        if (string(argv[i]) == "-f") {
            cout << "Accepting command line argument as filename to save to... ";
            filename = string(argv[i+1]);
            customSaveFile = true;
            cout << "will save to: " << filename << endl;
        }
    }
    if (customSaveFile == false) cout << "No filename given. Will save to default: out_default.csv" << endl; 
    
    //create a grid
    Grid g(1024, 1024, 20000);
    
    //maxSteps guarantees exiting after this many iterations
    int maxSteps = 3000000;
    bool update = true;
    
    //in this loop we update and print some information on the current state of the simulation
    cout << "Starting simulation...\n\n";
    while (g.stp < maxSteps) {
        if (update == true) {
            update = g.updateParticles();
            if (g.stp % 50 == 0) {    
                cout << "\rIteration: " << g.stp << " (max: " << maxSteps << ") | Particles: " << g.particles.size() << flush;
            }
        }
        //if updating is done, save to file and exit
        else {
            if (customSaveFile) g.saveToFile(filename);
            else g.saveToFile();
            cout << "Exiting...\n\n";
            break;
        }
    }
    return 0;
}
