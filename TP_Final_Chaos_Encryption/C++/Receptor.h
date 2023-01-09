//The receptor object can simulate a Lorenz system, coupled to an encrypted signal
//s(t), and recover (to an extent) the original message m(t) sent by the transmissor

#ifndef RECEPTOR_H
#define RECEPTOR_H

#include <vector>
#include <tuple>
#include <random>
#include <ctime>

class receptor {
    
    public:
        double r = 60.0, sigma = 10.0, b = 8.0/3.0; //system coefficients 
        double h; //RK4 stepsize
        double epsilon; //original amplitude modulation of m(t)
        double u;
        double v;
        double w;
        std::vector<double> ur; //u_r(t) trajectory
        std::vector<double> s;  //s(t)
        
        receptor(double freq, std::vector<double> s_, double epsilon_) {
            h = 1/freq;
            s = s_;
            epsilon = epsilon_;

            //random initial values of trajectory
            srand(time(NULL));
            u = doubleRand(-1, 1);
            v = doubleRand(-1, 1);
            w = doubleRand(-1, 1);
        }

        std::vector<double> run() {
            double ut;
            for (int i = 0; i < s.size(); i++) {
                std::tuple <double, double, double> uvw = RK4(h, u, v, w, s[i]);

                u = std::get<0>(uvw);
                v = std::get<1>(uvw);
                w = std::get<2>(uvw);

                ur.push_back((s[i]-u)/epsilon);
            }
            return ur;
        }

        ~receptor() {
        }

    private:
        //returns the three differential equations evaluated at u,v,w
        std::tuple <double, double, double> f(double u, double v, double w, double ut) {
            return std::make_tuple(sigma*(v-u), r*ut-v-20*ut*w, 5*ut*v-b*w);
        }
        //gives the next RK4 step from (u,v,w), using stepsize h_
        std::tuple <double, double, double> RK4(double h_, double u, double v, double w, double ut) {
            std::tuple <double, double, double> k1 = f(u,                 v,                 w,                 ut);
            std::tuple <double, double, double> k2 = f(u+h_*std::get<0>(k1)/2, v+h_*std::get<1>(k1)/2, w+h_*std::get<2>(k1)/2, ut);
            std::tuple <double, double, double> k3 = f(u+h_*std::get<0>(k2)/2, v+h_*std::get<1>(k2)/2, w+h_*std::get<2>(k2)/2, ut);
            std::tuple <double, double, double> k4 = f(u+h_*std::get<0>(k3),   v+h_*std::get<1>(k3),   w+h_*std::get<2>(k3),   ut);
            
            double u_ = u + h_*(std::get<0>(k1)+2*std::get<0>(k2)+2*std::get<0>(k3)+std::get<0>(k4))/6.0; // + o(h^5)
            double v_ = v + h_*(std::get<1>(k1)+2*std::get<1>(k2)+2*std::get<1>(k3)+std::get<1>(k4))/6.0; // + o(h^5)
            double w_ = w + h_*(std::get<2>(k1)+2*std::get<2>(k2)+2*std::get<2>(k3)+std::get<2>(k4))/6.0; // + o(h^5)
            
            return std::make_tuple(u_, v_, w_);
        }
        //returns a random double in range [min, max]
        double doubleRand(double min, double max) {
            double f = (double)rand() / RAND_MAX;
            return min + f*(max-min);
        }
};

#endif