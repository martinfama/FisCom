 import numpy as np
 
#transmissor solves the Lorenz system numerically using RK4
#Arguments: initial values (u_0, v_0, w_0)
#           initial stepsize (h)
#           final time (t_f)
def transmissor(u_0, v_0, w_0, t_f, freq):
    
    #Lorenz system coefficients
    r, sigma, b = 60, 10, 8/3
        
    #stepsize
    h = 1/freq

    #returns the three differential equations evaluated at u,v,w
    def f(u, v, w):
        return sigma*(v-u), r*u-v-20*u*w, 5*u*v-b*w
    
    #gives the next RK4 step from (u,v,w), using stepsize h_
    def RK4(h_, u, v, w):
        k_1_u, k_1_v, k_1_w = f(u,           v,           w)
        k_2_u, k_2_v, k_2_w = f(u+h_*k_1_u/2, v+h_*k_1_v/2, w+h_*k_1_w/2)
        k_3_u, k_3_v, k_3_w = f(u+h_*k_2_u/2, v+h_*k_2_v/2, w+h_*k_2_w/2)
        k_4_u, k_4_v, k_4_w = f(u+h_*k_3_u,   v+h_*k_3_v,   w+h_*k_3_w)
        
        u_ = u + h_*(k_1_u+2*k_2_u+2*k_3_u+k_4_u)/6 # + o(h^5)
        v_ = v + h_*(k_1_v+2*k_2_v+2*k_3_v+k_4_v)/6 # + o(h^5)
        w_ = w + h_*(k_1_w+2*k_2_w+2*k_3_w+k_4_w)/6 # + o(h^5)
        
        return u_, v_, w_

    #time points 
    t = np.arange(0, t_f, h)
    
    #trajectory of system (we will return these arrays (and t) at the end)
    u = [u_0]
    v = [v_0]
    w = [w_0]
    
    for j in range(1, len(t)):
        u_, v_, w_ = RK4(h, u[j-1], v[j-1], w[j-1])
            
        u.append(u_)
        v.append(v_)
        w.append(w_)
    
    return t, u, v, w
