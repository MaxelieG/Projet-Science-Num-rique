import numpy as np 

### Calculation of the infinite sum ###

#Calculate the term n of the sum
def term_of_the_sum (n, x, y, Lx, Ly):

    n_term = (1/(2*n + 1))*np.sin((2*n + 1)*np.pi*x/Lx)*np.sinh((Ly - y)*(2*n + 1)*np.pi/Lx)/np.sinh((2*n + 1)*np.pi*Ly/Lx)

    return n_term

#Calculate the sum
def infinite_sum (x, y, Lx, Ly):

    index = 0
    sum = 0

    while abs(term_of_the_sum(index, x, y, Lx, Ly) - term_of_the_sum(index + 1, x, y, Lx, Ly)) > 10**(-10):

        sum += term_of_the_sum(index, x, y, Lx, Ly)
        index += 1
    
    return sum

### Calculation of the analytical solution ###

def solution(x_list, y_list, Lx, Ly, temperature_T1, temperature_T2, temperature_list):

    for i, x in enumerate(x_list):

        for j, y in enumerate(y_list):

            temperature = temperature_T2 + 4*(temperature_T1 - temperature_T2)*infinite_sum(x, y, Lx, Ly)/np.pi 

            temperature_list[i, j] = temperature
    
    return temperature_list