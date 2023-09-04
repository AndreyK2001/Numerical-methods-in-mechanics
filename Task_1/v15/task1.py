import math
import numpy as np

# Значение alpha
ALPHA = 10
# Допустимая погрешность
eps = 10e-4
# Концы отрезка, на котором локализован корень
a = -6
b = 5
# Начальное приближение
x0 = -6

# Функция и её производные - вариант 15


def f(x):
    """ Функция f(x) с параметром alpha
    """
    return x + pow(math.e, -1 / (1 + x ** 2)) - ALPHA + 5


def f1(x):
    """ f'(x) - Первая производная f(x)
    """
    return 1 + (2 * x * pow(math.e, -1 / (1 + x ** 2))) / ((1 + x ** 2) ** 2)


def f2(x):
    """ f"(x) - Вторая производная f(x)
    """
    return (pow(math.e, -1 / (1 + x ** 2)) * (2 - 6 * (x ** 4))) / ((1 + x ** 2) ** 4)


def M():
    h = (b-a)*eps
    x = [a+i*h for i in range(int(1/eps+1))]
    f_x = list(map(f1, x))
    return np.linalg.norm(np.array(f_x), ord=np.inf)


def m():
    h = (b-a)*eps
    x = [a+i*h for i in range(int(1/eps+1))]
    f_x = list(map(f1, x))
    return np.linalg.norm(np.array(f_x), ord=-np.inf)


def sign(x):
    if (x == 0):
        return 0
    elif (x < 0):
        return -1
    else:
        return 1


def phi(x):
    """ Функция phi(x) представляет собой эквивалент уравнению: f(x) = 0 <=> phi(x) = x
    """
    return x - sign(f1(x))*f(x)/M()


def phi1(x):
    """ Функция phi(x) представляет собой эквивалент уравнению: f(x) = 0 <=> phi(x) = x
    """
    return 1 - abs(f1(x))/M()


def q():
    h = (b-a)*eps
    x = [a+i*h for i in range(int(1/eps+1))]
    f_x = list(map(phi1, x))
    return np.linalg.norm(np.array(f_x), ord=np.inf)

# Численные методы


def bisection(f, a, b, eps):
    # print('Расчет методом половинного деления:')
    n = 0

    if (b < a):  # Проверяем, что а - начало отрезка, b - конец
        a, b = b, a

    while abs(a - b) > eps:
        n += 1
        xi = (a + b) / 2
        if f(a) * f(xi) > 0:
            a = xi
        else:
            b = xi
        # print(f'n: {n}, x = {xi}')

    return [n, xi]


def hords(f, a, b, eps):
    # print('Расчет методом хорд:')
    n = 0
    x_prev = a
    x = b
    while abs(x - x_prev) > m()/(M()-m())*eps:
        n += 1
        temp = x
        x = x - f(x) * (x - x_prev) / (f(x) - f(x_prev))
        x_prev = temp
        # print(f'n: {n}, x = {x}')

    # print()
    return [n, x]


def iterations(phi, x0, eps):
    # print('Расчет методом простых итераций:')
    n = 0
    x = phi(x0)
    x_prev = x0
    while abs(x - x_prev) > (1-q())/q()*eps:
        n += 1
        x_prev = x
        x = phi(x)
        # print(f'n: {n}, x = {x}')

    # print()
    return [n, x]


def aitken(phi, a, b, eps):
    # print('Расчет методом Эткена:')
    n = 0
    x_prev = a
    x = phi(a)
    while abs(x - x_prev) > (1-q())/q()*eps:
        n += 1
        temp = x
        x = phi(x) - ((phi(x) - phi(x_prev)) * (phi(x) - x)) / \
            ((phi(x) - phi(x_prev)) - (x - x_prev))
        x_prev = temp
        # print(f'n: {n}, x = {x}')

    # print()
    return [n, x]


def newton(f, f1, x0, eps):
    # print('Расчет методом Ньютона:')
    n = 0
    x_prev = x0
    x = x0 - f(x0) / f1(x0)
    while abs(x - x_prev) > m()/(M()-m())*eps:
        n += 1
        x_prev = x
        x = x - f(x) / f1(x)
        # print(f'n: {n}, x = {x}')

    # print()
    return [n, x]


# Расчеты

with open("table.csv", "w") as file:
    file.write("alpha={}, [a b] = [{} {}] \n".format(ALPHA, a, b))
    file.write("Метод, Корень, Невязка, x0, Количество итераций, M, m, q \n")
    [n, x] = bisection(f, a, b, eps)
    file.write("Половинное деление, {}, {}, {}, {}, , , \n".format(x, f(x), x0, n))

    [n, x] = hords(f, a, b, eps)
    file.write("Метод хорд, {}, {}, {}, {}, {}, {}, \n".format(
        x, f(x), x0, n, M(), m()))

    [n, x] = iterations(phi, x0, eps)
    file.write("Метод простой итерации, {}, {}, {}, {}, , , {} \n".format(
        x, f(x), x0, n, q()))

    [n, x] = aitken(phi, a, b, eps)
    file.write("Метод Эткина, {}, {}, {}, {}, , , {} \n".format(
        x, f(x), x0, n, q()))

    [n, x] = newton(f, f1, x0, eps)
    file.write("Метод Ньютона, {}, {}, {}, {}, {}, {}, \n".format(
        x, f(x), x0, n, M(), m()))

    file.write("\n")

    file.close()
