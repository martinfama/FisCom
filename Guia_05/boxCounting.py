import numpy as np
                                      #d = box side length in number of cells
def boxCount(grid_Height, grid_Width, d, filename=""):
    #data should be a list given as:
    # x_1 y_1 t_1
    # x_2 y_2 t_2
    #     ...
    # x_n y_n t_n
    #where each row represents a point of the agreggate, and the time at which it stuck
    data = np.transpose(np.genfromtxt(filename))
    x, y, t = data[0], data[1], data[2]    
    B = 0.0

    #generate a grid to represent the boxes. each box has a value set to 0
    box_grid = []
    for j in range(int(grid_Height/d)):
        box_grid.append([])
        for i in range(int(grid_Width/d)):
            box_grid[-1].append(0)
    
    #here we go over the entire agreggate. for each cell, we check which box it is at
    #and set that box equal to 1
    for i in range(len(x)):
        X = int( x[i] // d)
        Y = int( y[i] // d)
        box_grid[Y][X] = 1
    
    #sum all instances of boxes that have value 1, therefore giving the boxcount
    B = sum([b.count(1) for b in box_grid])
    #and return that value
    return B

print(boxCount(1024, 1024, 4, "out.csv"))