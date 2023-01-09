#receptor solves the Lorenz system numerically using RK4 (adapative step)
#and coupled to the transmitter signal u_t(t)
#Arguments: initial values (u_0, v_0, w_0)
#           list of time points and trajectory of u of transmissor
def receptor(u_0, v_0, w_0, t, u_T):
    
    #Lorenz system coefficients
    r, sigma, b = 60, 10, 8/3

    #returns the three differential equations evaluated at u,v,w
    def f(u, v, w, uT):
        return sigma*(v-u), r*uT-v-20*uT*w, 5*uT*v-b*w
    
    #gives the next RK4 step from (u,v,w), using stepsize h_
    def RK4(h_, u, v, w, uT):
        k_1_u, k_1_v, k_1_w = f(u,           v,           w,              uT)
        k_2_u, k_2_v, k_2_w = f(u+h_*k_1_u/2, v+h_*k_1_v/2, w+h_*k_1_w/2, uT)
        k_3_u, k_3_v, k_3_w = f(u+h_*k_2_u/2, v+h_*k_2_v/2, w+h_*k_2_w/2, uT)
        k_4_u, k_4_v, k_4_w = f(u+h_*k_3_u,   v+h_*k_3_v,   w+h_*k_3_w,   uT)
        
        u_ = u + h_*(k_1_u+2*k_2_u+2*k_3_u+k_4_u)/6 # + o(h^5)
        v_ = v + h_*(k_1_v+2*k_2_v+2*k_3_v+k_4_v)/6 # + o(h^5)
        w_ = w + h_*(k_1_w+2*k_2_w+2*k_3_w+k_4_w)/6 # + o(h^5)
        
        return u_, v_, w_

    #trajectory of system (we will return these arrays at the end)
    u = [u_0]
    v = [v_0]
    w = [w_0]
    
    for j in range(0, len(t)-1):
        #next values with step h
        h = t[j+1]-t[j]
        u_1, v_1, w_1 = RK4(h, u[j], v[j], w[j], u_T[j])
            
        u.append(u_1)
        v.append(v_1)
        w.append(w_1)
    
    return u, v, w
