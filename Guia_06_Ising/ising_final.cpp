// Ising Model simulator
// possible compilation command: g++ -O3 -std=c++11 -o ising ising_final.cpp

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
#include <chrono>

using namespace std;

//random generator
random_device rd;
mt19937 gen(rd());

//mod function used to index grid as periodic
int mod(int x, int N) {
    return (x % N + N) % N;
}

class Grid {
    public:
        //N will be the grid side length in number of cells
        int N;
        double T;
        double E_tot = 0.0;
        double M_tot = 0.0;
        vector<vector<int>> grid;

        Grid(int N_, double T_) {       
            //seed random generator
            gen.seed(chrono::high_resolution_clock::now().time_since_epoch().count());
            N = N_;
            T = T_;
            //initialilze grid by filling it with spins
            for(int y = 0; y < N; y++) {
                vector<int> tmp_vec;
                grid.push_back(tmp_vec);
                for (int x = 0; x < N; x++) {
                    grid[y].push_back(2*uniform_int_distribution<int>{0, 1}(gen)-1); //this writes a -1 or a 1
                }
            }
            //initialize total Energy and Magentization
            for(int y = 0; y < N; y++) {
                for (int x = 0; x < N; x++) {
                    E_tot += spinE(y, x);
                    M_tot += grid[y][x];
                }
            }
        }

        //update the grid, and the total energy and magnetization
        void updateGrid() {
            //choose random point on grid
            int y = uniform_int_distribution<int>{0, N-1}(gen);
            int x = uniform_int_distribution<int>{0, N-1}(gen);
            //calculate cost of flipping
            double dE = flipEnergy(y, x);
            //calculate Boltzmann factor
            double p = min(1.0, exp(-dE/T));
            //generate a uniform double in range (0, 1)
            double w = uniform_real_distribution<>{0, 1}(gen);
            //check if flip is accepted, and if so, update grid, E and M
            if (w <= p) {
                grid[y][x] *= -1;
                E_tot += 2*spinE(y, x);
                M_tot += 2*grid[y][x];
            }
        }

        double E_per_spin() {
            return E_tot/N/N;
        }
        double M_per_spin() {
            return M_tot/N/N;
        }

        ~Grid() {}
    
    private:
        //returns the energy of spin at (y, x)
        double spinE(int y, int x) {
            return -grid[y][x] * ( grid[mod(y-1, N)][x] + grid[mod(y+1, N)][x] + grid[y][mod(x-1, N)] + grid[y][mod(x+1, N)] );
        }
        //returns the energy required to flip spin at (y, x)
        double flipEnergy(int y, int x) {
            return -2*spinE(y, x); //why? E_yx' = -E_yx -> dE = E' - E = -2*E_yx
        }

};

int main(int argc, char* argv[]) {
    
    //the following code will simulate a 32x32 grid at T=2.35
    //the energy per spin will be saved on each time unit passage
    //to file "E_espines.csv"

    //open file to save data to
    ofstream file;
    file.open("E_espines.csv", ios::out | ios::trunc);

    //create 32x32 grid, at T=2.35
    Grid g = Grid(32, 2.35);
    //run 1000 time units (1000*N^2 iterations)
    for (int i = 0; i < 1000 * g.N * g.N; i++) {
        g.updateGrid();
        //save E/spin to file when a whole time unit has passed
        if (i % (g.N*g.N) == 0) {
            file << g.E_per_spin() << " ";
        }
    }
    cout << "\nDone...\n";
    file.close();
    return 0;
}