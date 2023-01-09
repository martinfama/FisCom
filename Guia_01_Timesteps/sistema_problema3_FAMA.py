import numpy as np #importamos numpy para manejar arrays

#sistema nos resuelve el sistema de ecuaciones numéricamente
#se le pasan como argumento los parámetros a,b
#                           los valores iniciales x_0, y__0
#                           el paso h
#                         y hasta que tiempo t_f resolver
def sistema(a, b, x_0, y_0, h, t_f):
        
    #devuelve las dos ecuaciones diferenciales evaluadas
    def f(x, y):
        return a - (b+1)*x + x**2*y, b*x - x**2*y
    
    #array que contiene los intervalos de tiempo
    t = np.arange(0, t_f, h)
    
    #aca guardaremos la trayectoria del sistema
    x = np.array([x_0])
    y = np.array([y_0])
    
    #loopeamos sobre todos los intervalos de tiempo
    for j in range(0, len(t)):
        
        #calculamos los k para el metodo RK4
        k_1_x, k_1_y = f(x[j], y[j])
        k_2_x, k_2_y = f(x[j]+h*k_1_x/2, y[j]+h*k_1_y/2)
        k_3_x, k_3_y = f(x[j]+h*k_2_x/2, y[j]+h*k_2_y/2)
        k_4_x, k_4_y = f(x[j]+h*k_3_x, y[j]+h*k_3_y)
        
        #agregamos los nuevos valores calculados de x e y a los
        #arrays que contienen la trayectoria
        x = np.append(x, x[j] + h*(k_1_x+2*k_2_x+2*k_3_x+k_4_x)/6) # + o(h^5)
        y = np.append(y, y[j] + h*(k_1_y+2*k_2_y+2*k_3_y+k_4_y)/6) # + o(h^5)
        
    return x, y, t

#un ejemplo que corre el metodo y guarda a un .csv:
x, y, t = sistema(a=1, b=1, x_0=3, y_0=3, h=0.1, t_f=10, err=1e-05)
np.savetxt(f"fijo_a_1_b_1.csv", np.c_[t, x, y])