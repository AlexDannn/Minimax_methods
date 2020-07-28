import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import random
from scipy.signal import argrelextrema 
from numpy import array
import time 

'''
Fractional rational Remez method:

1) Accepts:
a) Approximate function, given discretely:
• x - a list containing a set of points on the approximation interval [-1, 1];
• y - a list containing the values ​​of the function at the points of the column vector x;
The number of points is formed based on the entered powers of the numerator and denominator of the approximating function m and n (m + n + 2);

b) m and n are integers denoting the degrees of the numerator and denominator of the approximating function.

2) Returns:
a) Error_func is the maximum value of the error function on the approximation interval;
b) The graph of the error function and the graph of the approaching and approximating function;
c) a_coeff, b_coeff are the coefficients of the numerator and denominator of the approximating function;
d) iter is the number of iterations performed by the method.
'''


start_time = time.time()

functions = {
   '1': lambda x: math.sin(25*x),
    '2': lambda x: math.sin(2.5*math.cos(x)),
    '3': lambda x: math.sin(x)*math.cos(x),
    '4': lambda x: math.cos(4*math.sin(x)),
    '5': lambda x: math.exp(math.sin(x)),
    '6': lambda x: math.exp(x),
    '7': lambda x: math.sinh(x),
    '8': lambda x: math.atan(x),
    '9': lambda x: x**6,
    '11': lambda x: math.log(x+3),

    '12': lambda x: 1/(x**2+1),
    '13': lambda x: (x**5+x**3+x)/(x**6+x**2+3),
    '14': lambda x: 1/(1+x+x**2),

    '16': lambda x: math.sin(x)/(x**2+1),
    '17': lambda x: (x**2+2)/(x**2+1),
    '18': lambda x: math.exp(-(x+1)**2)/(4*(x)**2 +  1.1),
    '19': lambda x: max(math.sin(20*x), math.exp(x-1)),

    '20': lambda x: math.tanh(50*x),
    '21': lambda x: math.tanh(x+0.5) - math.tanh(x-0.5),
    '22': lambda x: 1 - math.sin(5*abs(x - 0.5)),
    '23': lambda x: math.log(1.001 + x),
    '25': lambda x: x**3 + (x**(1/3.0)*math.exp(-x**2))/8, 
    '26': lambda x: (100*math.pi*(x**2-0.36))/(math.sinh(100*math.pi*(x**2-0.36))), 
    '27': lambda x: abs(x)*math.sqrt(abs(x)),
    '28': lambda x: (x-(x**3)/math.factorial(3)+ (x**5)/math.factorial(5))/(x**2+1)
}



def Rat_Remez(m, n, formula, x_equiv=None):
    mmErrors = []
    errvec = []
    N = n + m + 2
    a=-1
    b=1
    x, y = [], []
    E_0 = 0
    k = 0
    a_b_E_coeff = []
    a_coeff = np.zeros(m+1)
    b_coeff = np.zeros(n+1)
    delta = 1e-15

    # Узлы Чебышева 
    def func(x, formula):
        return formula(x)

    if x_equiv is None:
        for i in range(1,N+1):
            x.append((a+b)/2+(b-a)/2*math.cos((math.pi*(2*i-1)/(2*N))))
            y.append(
                func((a+b)/2+(b-a)/2*math.cos((math.pi*(2*i-1)/(2*N))), formula))
        x.reverse()
        y.reverse()
    else:
        x = np.linspace(a, b, N)
        y = [func(i, formula) for i in x]
    a_matr = np.zeros((N, N))
    y_matr = np.zeros(N)
    min_E = -1
    max_E = 5
    iter = 1

    while abs(max_E - min_E) > delta:
        for i in range(N):
            deg = 1
            y_matr[i] = y[i]
            for j in range(0, m+1):
                a_matr[i][j] = x[i]**j
            for j in range(m+1,N):
                a_matr[i][j] = (((-1)**i)*E_0-y[i])*x[i]**deg
                deg += 1
            a_matr[i][-1] = (-1)**i
        k += 1
        a_b_E_coeff = np.linalg.solve(a_matr, y_matr)
        print("Начальные точки: ")
        print(x)

        for i in range(N-1):
            if i < m+1:
                a_coeff[i] = a_b_E_coeff[i]
            else:
                b_coeff[i-m] = a_b_E_coeff[i]

        b_coeff[0] = 1 

        def rat_res(a_coeff, b_coeff, x):
            a = 0
            b = 1
            for i in range(len(a_coeff)):
                a+=a_coeff[i]*x**i
            for i in range(1,len(b_coeff)):
                b+=b_coeff[i]*x**i
            return a/b

        xnew = np.linspace(a, b, 1000)
        y_func = [func(i, formula) for i in xnew]

        ynew = [rat_res(a_coeff, b_coeff, i) for i in xnew]
        y_diff = [rat_res(a_coeff, b_coeff, i)-func(i, formula)
                    for i in xnew]
        
        extremums_x = []
        extremums_y = []
        s = pd.Series(y_diff)
        grp = s.groupby((np.sign(s).diff().fillna(0).ne(0)).cumsum())
        extremums_y_n = grp.apply(
            lambda x: x.abs().max()*np.sign(x[x.abs().idxmax()]))
        extremums_y = []
        print(len(extremums_y_n))
        for i in range(len(extremums_y_n)):
            # print(i)
            if extremums_y_n.iloc[i] != 0:
                extremums_y.append(extremums_y_n[i])

        print(extremums_y, 'Экстремумы y')
        print('----------------')
        extremums_x = []
        it = 0
        for i in xnew:
            if rat_res(a_coeff, b_coeff,i)-func(i, formula) == extremums_y[it]:
                extremums_x.append(i)
                if it != (len(extremums_y)-1):
                    it += 1
                else:
                    break

        print("Экстремумы")
        print(extremums_x)
        print(extremums_y)

        max_E = 0
        local_max = 0
        x_new = x

        x_old = x.copy()
        def closest(list, Number):
                aux = []
                for valor in list:
                    aux.append(abs(Number-valor))
                return aux.index(min(aux))

        if len(extremums_x) == len(x):
            for i in range(len(extremums_x)):
                if abs(extremums_y[i]) > abs(a_b_E_coeff[-1]):
                    x[i] = extremums_x[i]
                    y[i] = func(extremums_x[i], formula)
        else:
            for i in range(len(extremums_y)):

                if math.fabs(extremums_y[i]) > abs(a_b_E_coeff[-1]):

                    if extremums_x[i] not in x:
                        closest_x_index = closest(x, extremums_x[i])
                        x[closest_x_index] = extremums_x[i]

            for i in range(len(x)):
                y[i] = func(x[i], formula)

        if a_b_E_coeff[-1] > 0:
            Delta_plot = a_b_E_coeff[-1]
        else:
            Delta_plot = (-1)*a_b_E_coeff[-1]
        E_0 = a_b_E_coeff[-1]
        print('Iter number: ', iter)
        print('Error E: ', a_b_E_coeff[-1])
        extremums_y_l = []
        for i in range(len(extremums_y)):
            if extremums_y[i] < 0:
                extremums_y_l.append(abs(extremums_y[i]))
            else:
                extremums_y_l.append(extremums_y[i])
        
        max_E = max(extremums_y_l)
        min_E = min(extremums_y_l)
        print('Error func: ', max(extremums_y, key=abs))
        errvec.append(max(extremums_y))
        print('max_E: ', max_E)
        print('min_E: ', min_E)
        mmErrors.append(abs(max_E-min_E))
        print('Cur_delta: ', abs(max_E - min_E))
        y_test = np.zeros(len(x))

        iter += 1
        result = list(set(x) & set(x_old))
        Delta_plot_p = [Delta_plot for i in xnew]
        Delta_plot_n = [(-1)*Delta_plot for i in xnew]
        if len(result) == len(x):
            print("Совпадение с предыдущими значениями")

            break
        if iter == 51:

            break
        
        plt.figure(iter)
        plt.plot(xnew, y_diff, label='Error function')
        plt.plot(xnew, Delta_plot_p, '--', color='red', label='E')
        plt.plot(xnew, Delta_plot_n, '--', color='red')
        plt.plot(extremums_x, extremums_y, 'o', label='Extremums')
        plt.plot(x_new, y_test, 'o', label='X grid')
        plt.legend(fontsize='x-small',frameon=True,framealpha=1)
        plt.grid()    
        
    plt.show()
    print(errvec)
    print("--- %s seconds ---" % (time.time() - start_time))
    print(mmErrors)
    return 0
for i in range(11,12):
    #for j in range(0,2):
    Rat_Remez(5, 2, functions['6'])
