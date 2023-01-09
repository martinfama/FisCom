import numpy as np #importamos numpy para manejar arrays

#sistemaAdapt ahora tiene pasos adaptativos 
#se le pasan los mismos argumentos, mas la
#tolerancia (err)
def sistemaAdapt(a, b, x_0, y_0, h, t_f, err):
    
    #devuelve las dos ecuaciones diferenciales evaluadas
    def f(x, y):
        return a - (b+1)*x + x**2*y, b*x - x**2*y

    def RK4(h_, x, y):
        k_1_x, k_1_y = f(x             , y)
        k_2_x, k_2_y = f(x + h_*k_1_x/2, y + h_*k_1_y/2)
        k_3_x, k_3_y = f(x + h_*k_2_x/2, y + h_*k_2_y/2)
        k_4_x, k_4_y = f(x + h_*k_3_x  , y + h_*k_3_y)
        
        #calculamos los (posibles) nuevos valores de x, y
        #digo posibles porque mas adelante vamos a chequear
        #si cumplen con la cota del error para ser aceptados
        x_ = x + h_*(k_1_x+2*k_2_x+2*k_3_x+k_4_x)/6 # + o(h^5)
        y_ = y + h_*(k_1_y+2*k_2_y+2*k_3_y+k_4_y)/6 # + o(h^5)
        
        return x_, y_
    
    #array que contiene los intervalos de tiempo
    t = np.array([0.0])
    
    #aca guardaremos la trayectoria del sistema
    x = np.array([x_0])
    y = np.array([y_0])
    
    #en este while loopeamos hasta que el tiempo llegue a t_f
    #puede tener una cantidad arbitraria de pasos, ya que el 
    #largo del paso "h" se va modificando
    j = 0 #indice
    while t[-1] < t_f:
        while True:
            
            #proximos valores con paso h
            x_1, y_1 = RK4(h, x[j], y[j])
            #proximos valores con paso h/2
            #hay que hacer dos pasos aca
            x_2, y_2 = RK4(h/2, x[j], y[j])
            x_2, y_2 = RK4(h/2, x_2, y_2)
            dx = np.abs(x_2 - x_1)
            dy = np.abs(y_2 - y_1)
            delta = max(dx, dy)
            
            #aca chequeamos si se cumplen las cotas de errores
            if delta <= err/2:
                t = np.append(t, t[-1]+h/2)
                h *= 1.5
                break #salimos del while
            elif delta > err/2 and delta <= err:
                t = np.append(t, t[-1]+h/2)
                break #salimos del while
            
            #si ninguna condicion se cumplio, rechazamos x_2, y_2 y 
            #reducimos el paso
            h /= 1.5
        
        #si llegamos hasta aca, nos aseguramos que queremos aceptar x_2, y_2
        x = np.append(x, x_2) # + o(h^5)
        y = np.append(y, y_2) # + o(h^5)
        j += 1
        
    return x, y, t

#un ejemplo que corre el metodo y guarda a un .csv:
x, y, t = sistema(a=1, b=3, x_0=3, y_0=3, h=0.01, t_f=10, err=1e-05)
np.savetxt(f"adapt_a_1_b_3.csv", np.c_[t, x, y])