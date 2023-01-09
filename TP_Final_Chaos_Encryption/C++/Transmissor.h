//The transmissor object can simulate a Lorenz system, and add a 
//message m(t) to the signal to be sent as an encrypted message

#ifndef TRANSMISSOR_H
#define TRANSMISSOR_H

#include <vector>
#include <tuple>
#include <random>
#include <ctime>

class transmissor {
    
    public:
        double r = 60.0, sigma = 10.0, b = 8.0/3.0; //system coefficients
        double h; //RK4 step size
        double t_f; //final time
        double t_grace; //grace time: amount of time given to synchronize systems
        int sampleRatio; //atractor sampling distance
        double epsilon; //amplitude modulation of m(t)
        double u; 
        double v;
        double w;
        std::vector<double> ut; //u_t(t) trajectory
        std::vector<double> m;  //m(t)

        transmissor(double t_grace_, double freq, std::vector<double> m_, int sampleRatio_, double epsilon_) {
            h = 1/freq;
            m = m_;
            t_grace = t_grace_;
            sampleRatio = sampleRatio_;
            epsilon = epsilon_;

            //random initial values of trajectory
            srand(time(NULL));
            u = doubleRand(-1, 1);
            v = doubleRand(-1, 1);
            w = doubleRand(-1, 1);
        }

        std::vector<double> run() {
            //simulate trajectory for t_grace, without adding m(t)
            for (int i = 0; i < t_grace/h; i++) {
                std::tuple <double, double, double> uvw = RK4(h, u, v, w);

                u = std::get<0>(uvw);
                v = std::get<1>(uvw);
                w = std::get<2>(uvw);

                ut.push_back(u);
            }
            //continue simulating trajectory, and now adding epsilon*m(t)
            for (int i = 0; i < m.size()*sampleRatio; i++) {
                std::tuple <double, double, double> uvw = RK4(h, u, v, w);

                u = std::get<0>(uvw);
                v = std::get<1>(uvw);
                w = std::get<2>(uvw);

                if (i % sampleRatio == 0) ut.push_back(u + epsilon*m[i/sampleRatio]);
                else ut.push_back(u);
            }
            return ut;
        }

        ~transmissor() {}

    private:
        //returns the three differential equations evaluated at u,v,w
        std::tuple <double, double, double> f(double u, double v, double w) {
            return std::make_tuple(sigma*(v-u), r*u-v-20*u*w, 5*u*v-b*w);
        }
        //gives the next RK4 step from (u,v,w), using stepsize h_
        std::tuple <double, double, double> RK4(double h_, double u, double v, double w) {
            std::tuple <double, double, double> k1 = f(u,           v,           w);
            std::tuple <double, double, double> k2 = f(u+h_*std::get<0>(k1)/2, v+h_*std::get<1>(k1)/2, w+h_*std::get<2>(k1)/2);
            std::tuple <double, double, double> k3 = f(u+h_*std::get<0>(k2)/2, v+h_*std::get<1>(k2)/2, w+h_*std::get<2>(k2)/2);
            std::tuple <double, double, double> k4 = f(u+h_*std::get<0>(k3),   v+h_*std::get<1>(k3),   w+h_*std::get<2>(k3));
            
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