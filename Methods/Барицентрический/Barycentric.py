import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import random
from scipy.signal import argrelextrema 
from numpy import array
import time 

start_time = time.time()

'''
Barycentric method:
1) Accepts:
    a) Approximate function, given discretely:
        • x - a list containing a set of points on the approximation interval [-1, 1];
        • y - a list containing the values ​​of the function at the points of the column vector x;
        The number of points is formed on the basis of the entered powers of the numerator and denominator of the approximating barycentric function n and n (2n + 2);
    b) n is an integer denoting the degree of the numerator and denominator of the approximating function.
2) Returns:
    a) Error_func is the maximum value of the error function on the approximation interval;
    b) The graph of the error function and the graph of the approaching and approximating function;
    c) iter is the number of iterations performed by the method.
'''


functions = {
    '1': lambda x: math.sin(2*x),
    '2': lambda x: math.sin(2.5*math.cos(x)),
    '3': lambda x: math.sin(x)*math.cos(x),
    '4': lambda x: math.cos(4*math.sin(x)),
    '5': lambda x: math.exp(math.sin(x)),
    '6': lambda x: math.exp(x),
    '7': lambda x: math.sinh(x),
    '8': lambda x: math.atan(x),
 '9': lambda x: x**6,
    # '10': lambda x: x**7+3*x**5+2*x-1,
    '11': lambda x: math.log(x+3),
    '12': lambda x: 1/(x**2+1),
    '13': lambda x: (x**5+x**3+x)/(x**6+x**2+3),
    '14': lambda x: 1/(1+x+x**2),
    '16': lambda x: math.sin(x)/(x**2+1),
    '17': lambda x: (x**2+2)/(x**2+1),
    '18': lambda x: math.exp(-(x-1)**2/(2*.5**2)),
    '19': lambda x: math.sin(25*x),
    '27': lambda x: abs(x)*math.sqrt(abs(x)),
    '20': lambda x: math.sin(math.exp(30*x)),
    '21': lambda x: (x-(x**3)/math.factorial(3)+ (x**5)/math.factorial(5))/(x**2+1)
}

def Barycentric(n, formula, x_equiv=None):
    errvec = []
    mmErrors = []
    N = n + 2
    a=-1
    b=1
    x, y = [], []

    Delta = 1e-16

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
    iter = 0

    while abs(max_E - min_E) > Delta:
        w = []
        Tau = (b-a)/4
        for i in range(len(x)):
            res_top = 1
            res_bottom = 0 
            for j in range(len(x)):
                if i != j:
                    res_top *= np.sign(x[i]-x[j])
            for j in range(len(x)):
                if i != j:
                    res_bottom += math.log(abs(x[i]-x[j]))
            w.append(res_top/math.exp(len(x)*math.log(1/Tau) + res_bottom))
        
        
        if iter == 0:
            h_t_start = 0
            h_b_start = 0

            for i in range (len(x)):
                h_t_start += w[i]*func(x[i], formula)
                h_b_start += (-1)**i*w[i]
            delta_start = h_t_start/h_b_start
            delta = delta_start
        if iter > 0:
            h_t = 0
            h_b = 0

            for i in range (len(x)):
                h_t += w[i]*func(x[i], formula)
                h_b += (-1)**i*w[i]
            delta = h_t/h_b

        if iter > 0  and abs(delta) > abs(delta_start):
            delta_start = delta 
        
        p = []
        for i in range(len(x)):
            p.append((-1)**i*delta_start + func(x[i], formula))

        def rat_res(X):
            top = 0
            bottom = 0
            for i in range(len(x)):
                if X != x[i]:
                    top += w[i]*p[i]/(X-x[i])
                    bottom += w[i]/(X-x[i])
                else:
                    return (rat_res(X-1e-15)+rat_res(X+1e-15))/2
            return top/bottom

        xnew = np.linspace(a, b, 1000)
        y_func = [func(i, formula) for i in xnew]
        ynew = [rat_res(i) for i in xnew]
        y_diff = [rat_res(i)-func(i, formula)
                    for i in xnew]
        
        extremums_x = []
        extremums_y = []
        y_diff = array(y_diff)
        
        def check_for_alter(y_for_check):
            for i in range(len(y_for_check)-1):
                if y_for_check[i]*y_for_check[i+1] >= 0:
                    print("Not alternant")
                    return False
            print("Alternant")
            return True

        s = pd.Series(y_diff)
        grp = s.groupby((np.sign(s).diff().fillna(0).ne(0)).cumsum())
        extremums_y_n = grp.apply(
            lambda x: x.abs().max()*np.sign(x[x.abs().idxmax()]))
        extremums_y = []
        for i in range(len(extremums_y_n)):
            if extremums_y_n.iloc[i] != 0:
                extremums_y.append(extremums_y_n[i])

        extremums_x = []
        it = 0
        for i in xnew:
            if rat_res(i)-func(i, formula) == extremums_y[it]:
                extremums_x.append(i)
                if it != (len(extremums_y)-1):
                    it += 1
                else:
                    break
        errvec.append(max(abs(y_diff)))

        check_for_alter(extremums_y) 


        extremums_y_abs = extremums_y.copy()
        for i in range(len(extremums_y_abs)):
            extremums_y_abs[i] = abs(extremums_y_abs[i])
        
        def del_min_index(extremums_y_abs):

            Maxum = 10000    
            for i in range(1,len(extremums_y_abs)-2):
                if abs(extremums_y_abs[i]) + abs(extremums_y_abs[i+1]) < Maxum:
                    Maxum = abs(extremums_y_abs[i]) + abs(extremums_y_abs[i+1])
                    min_index_1 = i
                    min_index_2 = i+1
            del extremums_y_abs[min_index_2]
            del extremums_y_abs[min_index_1]
            del extremums_y[min_index_2]
            del extremums_y[min_index_1]
            del extremums_x[min_index_2]
            del extremums_x[min_index_1]
        
        

        N1 = len(extremums_x)
        N2 = len(x)

        if len(x) != len(extremums_x):
            for i in range( (N1 - N2) // 2):
                del_min_index(extremums_y_abs)
                a = extremums_x[0]
                b = extremums_x[-1]

        if check_for_alter(extremums_y) is False:
            print("Not an alternant")
            #break

        index_y = []
        for i in range(len(extremums_y)-1):
            if extremums_y[i] > 0 and extremums_y[i+1] > 0:
                if extremums_y[i] >= extremums_y[i+1]:
                    index_y.append(i+1)
                else:
                    index_y.append(i)
            if extremums_y[i] < 0 and extremums_y[i+1] < 0:
                if extremums_y[i] <= extremums_y[i+1]:
                    index_y.append(i+1)
                else:
                    index_y.append(i)
        for i in reversed(index_y):
            del extremums_y[i]
            del extremums_x[i] 
        
    

        max_E = 0
        local_max = 0

        x_old = x.copy()

        def closest(list, Number):
                aux = []
                for valor in list:
                    aux.append(abs(Number-valor))
                return aux.index(min(aux))

        if len(extremums_x) == len(x):
            for i in range(len(extremums_x)):
                if abs(extremums_y[i]) > abs(delta_start):
                    x[i] = extremums_x[i]
                    y[i] = func(extremums_x[i], formula)
        else:
            for i in range(len(extremums_y)):
                if math.fabs(extremums_y[i]) > abs(delta_start):
                    if extremums_x[i] not in x:
                        closest_x_index = closest(x, extremums_x[i])
                        x[closest_x_index] = extremums_x[i]

        for i in range(len(x)):
            y[i] = func(x[i], formula)

        x_new = x.copy()
        result = list(set(x) & set(x_old))
        if len(result) == len(x):
            break

        E_0 = delta_start
        extremums_y_l = []
        for i in range(len(extremums_y)):
            if extremums_y[i] < 0:
                extremums_y_l.append(abs(extremums_y[i]))
            else:
                extremums_y_l.append(extremums_y[i])
        
        max_E = max(extremums_y_l)
        min_E = min(extremums_y_l)
        print('Error func: ', max(extremums_y, key=abs))
        mmErrors.append(abs(abs(max_E) - abs(min_E)))
        y_test = np.zeros(len(x))
        iter += 1
        if iter == 10:
            break
        if delta_start < 0:
            delta_n = delta_start
            delta_p = -1*delta_start
        else:
            delta_p = delta_start
            delta_n = -1*delta_start
        y_delta_n = [delta_n for i in xnew]
        y_delta_p = [delta_p for i in xnew]
        
    plt.figure(iter)
    plt.plot(xnew, y_diff, label='Error function')
    plt.plot(xnew, y_delta_p, '--', color='red', label='E')
    plt.plot(xnew, y_delta_n, '--', color='red')
    plt.plot(extremums_x, extremums_y, 'o', label='Extremums')
    plt.plot(x_new, y_test, 'o', label='X grid')
    plt.legend(fontsize='x-small',frameon=True,framealpha=0.5, loc='upper right')
    plt.grid()
        
    
    plt.figure(iter+1)
    plt.plot(xnew, y_func, label='f(x)')
    plt.plot(xnew, ynew, label='g(x)')
    plt.legend(loc='center left', fontsize='x-small')
    plt.grid()       
    plt.show()
    
    print(errvec)
    print("--- %s seconds ---" % (time.time() - start_time))
    print(mmErrors)
    return 0
for i in range(81,82, 2):
    Barycentric(18, functions['5'])
