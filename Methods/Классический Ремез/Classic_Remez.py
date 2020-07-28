import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import random
from fpdf import FPDF
import time 

start_time = time.time()
'''
Классический алгоритм Ремеза.
1)	Принимает:
    а)	 Приближаемую функцию, заданную дискретно:
        •	x – список, содержащий набор точек на интервале аппроксимации [-1, 1];
        •	y – список, содержащий значения функции в точках списка x;
    Количество точек формируется исходя из степени приближающего полинома n (n + 2); 

    б)	n – целое число, обозначающее степень приближающего полинома.

2)	Возвращает:
    а)	best_result_max – максимальное значение функции ошибки на интервале аппроксимации;
    б)	График функции ошибки и график приближающей и приближаемой функции;
    в)	b_coeff – коэффициенты полинома наилучшего приближения;
    г)	iter_n – количество итераций, совершенных методом.

'''


functions = {
    # Периодические функции
    '1': lambda x: math.sin(2*x),
    '2': lambda x: math.sin(2.5*math.cos(x)),
    '3': lambda x: math.sin(x)*math.cos(x),
    '4': lambda x: math.cos(4*math.sin(x)),
    '5': lambda x: math.exp(math.sin(x)),
    # Монотонные функции
    '6': lambda x: math.exp(x),
    '7': lambda x: math.sinh(x),
    '8': lambda x: math.atan(x),
    '9': lambda x: x**6,
    # '10': lambda x: x**7+3*x**5+2*x-1,
    '11': lambda x: math.log(x+3),
    # Дробно-рациональные функции
    '12': lambda x: 1/(x**2+1),
    '13': lambda x: (x**5+x**3+x)/(x**6+x**2+3),
    '14': lambda x: 1/(1+x+x**2),
    # '15': lambda x: math.sin(x)/x, '''Избегать нуля! '''
    '16': lambda x: math.sin(x)/(x**2+1),
    '17': lambda x: (x**2+2)/(x**2+1),
    '18': lambda x: math.exp(-(x-1)**2/(2*.5**2)),
    '19': lambda x: math.sin(25*x),
    '27': lambda x: abs(x)*math.sqrt(abs(x)),
    '20': lambda x: math.sin(math.exp(30*x)),
     '28': lambda x: (x-(x**3)/math.factorial(3)+ (x**5)/math.factorial(5))/(x**2+1)
}


def Remez(n, formula, x_equiv=None):
    errvec = []
    mmErrors = []
    print('Степень многочлена n: ', n)
    x = []
    y = []
    a = -1
    b = 1
    N = n+1
    if x_equiv is None:
        a_matr = np.zeros((n+2, n+2))
        m = np.zeros(n+2)
    iter_n = 0
    Error = 1
    delta = 1e-12

    def func(x, formula):
        return formula(x)

    if x_equiv is None:
        # Узлы Чебышева
        for i in range(n+2):
            x.append((a+b)/2+(b-a)/2*math.cos((math.pi*(n+1-i)/(n+1))))
            y.append(
                func((a+b)/2+(b-a)/2*math.cos((math.pi*(n+1-i)/(n+1))),
                     formula))
    else:
        # Равноудаленные узлы
        x = np.linspace(a, b, x_equiv)
        y = [func(i, formula) for i in x]
        n = x_equiv-2
        a_matr = np.zeros((n+2, n+2))
        m = np.zeros(n+2)
        print(n, 'n')

    max_E_prev = 0
    local_max_prev = 0
    local_min_prev = 0
    max_E_1_prev = 0
    min_E_1_prev = 0

    while Error > delta:
        print(x, 'Начальные точки x')
        print(y, 'Начальные точки y')
        for i in range(len(x)):
            m[i] = func(x[i], formula)
            for j in range(len(x)):
                a_matr[i][j] = x[i]**j
            a_matr[i][len(x)-1] = (-1)**i

        # print(a_matr, 'матрица')
        b_coeff_E = np.linalg.solve(a_matr, m)  # Находим коэффициенты b  и E
        print(np.linalg.solve(a_matr, m))
        b_coeff = np.zeros(len(x))
        for i in range(len(x)):
            b_coeff[i] = b_coeff_E[i]

        # print(b_coeff)
        print('NNNN ', n, 'l_x', len(x))

        def b_res(X, b_coeff):
            sum = 0
            for i in range(n+1):
                sum += X**i*b_coeff[i]
            return sum

        xnew = np.linspace(a, b, 1000)
        y_func = [func(i, formula) for i in xnew]  # Значение функции
        # Построение интерполянта Лагрнжа по n-1 точкам
        ynew = [b_res(i, b_coeff) for i in xnew]
        y_diff = [b_res(i, b_coeff)-func(i, formula)
                  for i in xnew]  # Разница между Лагранжом и функцией

        # Поиск экстремумов
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
            if b_res(i, b_coeff)-func(i, formula) == extremums_y[it]:
                extremums_x.append(i)
                # print(it,' ',len(extremums_y),' ',extremums_y[it])
                if it != (len(extremums_y)-1):
                    it += 1
                else:
                    break

        print(extremums_x, 'Экстремумы x')
        print('----------------')
        # '''

        extremums_y_p = []
        extremums_y_n = []
        extremums_x_p = []
        extremums_x_n = []

        for i in range(len(extremums_y)):
            if extremums_y[i] > 0:
                extremums_y_p.append(extremums_y[i])
                extremums_x_p.append(extremums_x[i])
            else:
                extremums_y_n.append(extremums_y[i])
                extremums_x_n.append(extremums_x[i])

        min_E = 1000000
        max_E = 0
        local_max = 0
        local_min = 0
        max_E_1 = 0
        min_E_1 = 0

        if(len(extremums_y) >= n+2):
            print("@")
            for i in range(len(extremums_y)):
                #  print(extremums_y[i],' ',max_E)
                if math.fabs(extremums_y[i]) > max_E:
                    max_E = math.fabs(extremums_y[i])
                    local_max = i
                    if extremums_y[i] < 0:
                        max_E_1 = max_E*(-1)
                    else:
                        max_E_1 = max_E

            for i in range(len(extremums_y)):
                # print(extremums_y[i],' ',min_E)
                if math.fabs(extremums_y[i]) < min_E:
                    min_E = math.fabs(extremums_y[i])
                    local_min = i
                    if extremums_y[i] < 0:
                        min_E_1 = min_E*(-1)
                    else:
                        min_E_1 = min_E

        print('E максимальное', local_max, '  ', max_E_1,' ', extremums_x[local_max])
        print('E минимальное', local_min, '  ', min_E_1)
        Error = math.fabs(math.fabs(max_E_1) - math.fabs(min_E_1))
        mmErrors.append(Error)
        print(Error, ' Error')
        print('l_m', local_max, 'len_x', len(x), ' extre_len ',len(extremums_x))

        
        def closest(list, Number):
                aux = []
                for valor in list:
                    aux.append(abs(Number-valor))
                return aux.index(min(aux))
        if len(extremums_x) > len(x):
            closest_x_index = closest(x, extremums_x[local_max])
            x[closest_x_index] = extremums_x[local_max]
        else:
            if local_max == len(x):
                x[-1] = extremums_x[local_max]
                b = x[-1]
            elif local_max == 0:
                x[0] = extremums_x[local_max]
                a = x[0]
            else:
                x[local_max] = extremums_x[local_max]
        
        print("==========================")
        print("Конец итерации номер ", iter_n+1)
        print("==========================") 
        y_o = np.zeros(len(x))
        iter_n += 1
        
        if(iter_n > 30):
            print('Лучший результат за 50 итераций:', Error)
            break

        if(max_E_prev == max_E and local_max_prev == local_max and
           local_min_prev == local_min and max_E_1_prev == max_E_1 and
           min_E_1_prev == min_E_1):
            print("Совпадение с предидущими значениями")
            break

        max_E_prev = max_E
        local_max_prev = local_max
        local_min_prev = local_min
        max_E_1_prev = max_E_1
        min_E_1_prev = min_E_1
        if(Error < delta):
            print('Лучший результат за ', iter_n, ' итераций:', Error)

    '''
    print('Коэффициенты многочлена наилучшего приближения')
    for i in range(n+1):
        print(b_coeff[i], 'x^', i, end=' ')
    '''
    print('Максимальное отклонение:', max(extremums_y))
    errvec.append(extremums_y[0])
    print("--- %s seconds ---" % (time.time() - start_time))
    print(mmErrors)
    plt.figure(iter_n)
    plt.plot(xnew, y_diff, extremums_x, extremums_y, 'o',x,y_o,'o', extremums_x_p,
                 extremums_y_p, '--', extremums_x_n, extremums_y_n, '--')
    plt.grid()
    return max(extremums_y)


# -----------------------------------------------------------------------------#
# Конец алгоритма Ремеза
# -----------------------------------------------------------------------------#
best_result_max = 1000

for i in range(11, 12):
    best_result = Remez(9, functions['6'])
    if best_result < best_result_max:
        best_result_max = best_result
        n_degree = i
    print()
    plt.show()
        



