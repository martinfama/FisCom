#parametros:
# a da la amplitud de los valores iniciales aleatorios de u y v.
# D es el coeficiente de difusion (se uso mu = D, nu = D/2).
# t_f el tiempo final medido en intervalos (es decir numero entero).
# dt y dx el tamaño de los intervalos de la grilla.
# seed es para iniciar el generador de numeros aleatorios (por default, seed = None, 
#                     pero si se pasa un valor fijo se pueden reproducir resultados).
# writeTimes es un array de numeros enteros que le dicen a la funcion en que tiempos guardar los datos
#    (por default, siempre guardamos los valores iniciales de t=0).
def morfogenesis(a, D, t_f, dt, dx, seed = None, writeTimes=[]):
    
    #seedeamos el generador de numeros aleatorios
    np.random.seed(seed)
    
    def g_1(u, v):
        return (-7*u*u - 50*u*v + 57)/32
    def g_2(u, v):
        return (7*u*u + 50*u*v - 2*v - 55)/32
    
    mu = D
    nu = D/2
    
    #cantidad de puntos en el intervalo de x (si dx = 0.01, hay 101 puntos entre 0 y 1 inclusive)
    xn = int( (1+dx)/dx )
    
    #estos dos arrays van a guardar los datos de los tiempos que pidamos en writeTimes
    U = []
    V = []
    
    #cada funcion tiene asociada dos estados: el anterior y el proximo. es decir, un u[0] y v[0]
    #se calcula el proximo paso, el cual se guarda en u[1] y v[1]. 
    u = [np.zeros(xn), np.zeros(xn)]
    v = [np.zeros(xn), np.zeros(xn)]
    #lista de puntos de x
    x = np.arange(0, 1+dx, dx)
    
    #valores inciales aleatorios de u y v
    for j in range(0, xn):
        u[0][j] = 1 + a*(2*np.random.uniform(0, 1) - 1)
        v[0][j] = 1 + a*(2*np.random.uniform(0, 1) - 1)
    #por default, guardamos los datos en t=0
    U.append(u[0][:])
    V.append(v[0][:])
    
    #la simulacion empieza aca
    for n in range(0, t_f+1):
        for j in range(0, xn):
            u[1][j] = u[0][j] + dt/(dx*dx) * mu * (u[0][(j+1)%101] - 2*u[0][j] + u[0][(j-1)%101]) + dt*g_1(u[0][j], v[0][j])
            v[1][j] = v[0][j] + dt/(dx*dx) * nu * (v[0][(j+1)%101] - 2*v[0][j] + v[0][(j-1)%101]) + dt*g_2(u[0][j], v[0][j])
        
        #aca chequeamos si n es uno de los tiempos que pedimos para guardar, 
        #y si es asi, anexamos los valores actuales a U y V
        if n in writeTimes:
            u_ = []
            v_ = []
            for j in range(0, xn):
                u_.append(u[1][j])
                v_.append(v[1][j])
            U.append(u_[:])
            V.append(v_[:])
        
        #actualizamos u[0] y v[0] con los nuevos valores calculados 
        u[0] = u[1][:]
        v[0] = v[1][:]
        
        #un contador que sirve para ver el progreso de la simulacion
        #hace que vaya mas lento, por eso esta comentado y queda a disposicion
        #del usuario si usarlo o no
        #print(f"\r{100*(n+1)/(t_f-1):.2f}% | {n}",end="")
        
    #le devolvemos al pedido de morfogenesis() el array con los puntos de x, y los valores de U y V
    return x, U, V


#un ejemplo de como llamar a la funcion:
x, U, V = morfogensis( a=0.1, D=0.00075, t_f=200000, dt=0.001, dx=0.01, seed=100, writeTimes=[2500, 25000, 100000, 200000] )
