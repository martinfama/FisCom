import numpy as np

#funcion del mapeo logistico
def f(x, r):
    return r*x*(1- x)

#genera una serie de puntos del mapeo logistico
#devuelve dos arrays: los puntos n, y los correspondientes x
def mL(r = 0.0, x0 = 0.0, N=0):
    x = np.array([x0])
    n = np.arange(0, N, 1)
    
    for i in range(1, N):
        x = np.append(x,  f(x[-1], r))
        
    return n, x

#genera un conjunto de puntos que forman el diagrama de bifurcaciones
#throwN le dice cuantos valores iniciales de N descartar, asi nos sacamos de encima el transitorio
#encontre que funciona bien para N = 1000, throwN = 500
def bifurcaciones(dr=0.01, r_min=0, r_max=4, N=100, throwN=10):
    r = []
    atractores = []
    
    for r_ in np.arange(r_min, r_max+dr, dr):
        atractores.append([])                      #esto de aca incrementa un poco N al aumentar r
        n, x = mL(r=r_, x0=np.random.uniform(0,1), N=int(N+(r_-r_min)*40))
        atractores[-1] = x[throwN:] #tiramos los primeros throwN puntos
        #para cada atractor encontrado, adjuntamos r_ de nuevo a r, asi nos quedan la misma cantidad de puntos en ambos
        for j in range(len(atractores[-1])):
            r.append(r_)
    
    return r, np.concatenate(atractores) #lo ultimo aplasta la lista de atractores (que es un conjunto de listas), en una sola dimension

#calcula el exponente de lyapunov en el intervalo r_min, r_max
#funciona mas o menos similar al de bifurcaciones (es decir throwN descarta el transitorio)
def lyapunov(dr, r_min, r_max, N, throwN):
    r = []
    lya = []
    
    for r_ in np.arange(r_min, r_max+dr, dr):
        r.append(r_)
        n, x = mL(r=r_, x0=np.random.uniform(0,1), N=int(N+(r_-r_min)*40))
        x = x[throwN:]
        #esto que sigue es la sumatoria con la que esta definido el exponente de lyapunov
        lya.append(0.0)
        for x_ in x:                   #esto seria f'(x)
            lya[-1] += np.log( np.abs( r_*(1 - 2*x_) ) )
        lya[-1] /= (len(x)) #dividimos por la cantidad de puntos en la serie 
    
    #devuelve dos arrays, los puntos considerados en r, y el correspondiente exponente calculado
    return r, lya