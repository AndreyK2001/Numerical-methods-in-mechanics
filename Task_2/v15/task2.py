import math
import numpy as np

# Допустимая погрешность
eps = 1e-4
# Концы отрезка, на котором локализован корень

# Начальное приближение
x0 = 1
y0 = 1

x_min = 0.3
y_min = 0.7
x_max = 0.5
y_max = 0.9

# Функции и их производные - вариант 15


def f(x, y):
    return x * x + y * y - 2 * x


def g(x, y):
    return x - math.e ** (-y)


def der_f_x(x, y):
    """ d/dx f(x, y) = 2x - 2
    """
    return 2 * x - 2


def der_f_y(x, y):
    """ d/dy f(x, y) = 2y
    """
    return 2 * y


def der_g_x(x, y):
    """ d/dx g(x, y) = 1
    """
    return 1


def der_g_y(x, y):
    """ d/dy g(x, y) = e^(-y)
    """
    return math.e ** (-y)


def phi(x, y):
    """ phi(x,y) = x <=> f(x,y) = 0
    """
    return math.sqrt(2 * x - x*x)  # (x*x + y*y) / 2


def dphidx(x, y):
    return (1-x)/(math.sqrt(2*x-x**2))  # (x*x + y*y) / 2


def gamma(x, y):
    """ gamma(x,y) = y <=> g(x,y) = 0
    """
    return pow(math.e, -y)  # math.log(1 / x)


def dgammady(x, y):
    return -pow(math.e, -y)  # math.log(1 / x)
# Числовые методы


def q():
    D = np.zeros((10, 10))
    step_x = (x_max-x_min)/10
    step_y = (y_max-y_min)/10

    for i in range(10):
        for j in range(10):
            D[i][j] = np.linalg.norm(np.array([[0, abs(dgammady(x_min+i*step_x, y_min+j*step_y))],
                                               [abs(dphidx(x_min+i*step_x, y_min+j*step_y)), 0]]), ord=np.inf)
    # print(D)
    return np.max(D)


def mu():
    D = np.zeros((10, 10))
    step_x = (x_max-x_min)/10
    step_y = (y_max-y_min)/10

    for i in range(10):
        for j in range(10):
            x = x_min+i*step_x
            y = y_min+j*step_y
            A = np.array([[der_f_x(x, y), der_f_y(x, y)],
                          [der_g_x(x, y), der_g_y(x, y)]])
            D[i][j] = np.linalg.norm(A)*np.linalg.norm(np.linalg.inv(A))

    return np.max(D)


def iterations(phi, gamma, x0, y0, eps):
    print('Расчет методом простых итераций:')
    n = 0
    x_prev = x0
    y_prev = y0
    y = phi(x0, y0)
    x = gamma(x0, y0)
    vector_prev = np.array([x_prev, y_prev])
    vector = np.array([x, y])
    norm = np.linalg.norm(vector - vector_prev, ord=np.inf)

    #print(f'n: {n}, x0 = {x0}, y0 = {y0}')
    #print(f'x1 = {x}, y1 = {y}')
    #print(f'diff = {vector - vector_prev}, norm = {norm}')
    # print()

    while norm > (1-q())/q()*eps:
        n += 1
        #print(f'n: {n}')
        x_now, y_now = x, y
        y = phi(x_now, y_now)
        x = gamma(x_now, y_now)
        x_prev, y_prev = x_now, y_now

        vector_prev = np.array([x_prev, y_prev])
        vector = np.array([x, y])
        norm = np.linalg.norm(vector - vector_prev, ord=np.inf)

        #print(f'x = {x_now}, y = {y_now}')
        #print(f'x+1 = {x}, y+1 = {y}')
        #print(f'diff = {vector - vector_prev}, norm = {norm}')
        # print()

    print(f'Окончательное решение на итерации {n}: (x, y) = ({x}, {y})')
    print()
    return(x, y, n)


def seidel(phi, gamma, x0, y0, eps):
    print('Расчет методом Зейделя:')
    n = 0
    x_prev = x0
    y_prev = y0
    y = phi(x0, y0)
    x = gamma(x0, y0)
    vector_prev = np.array([x_prev, y_prev])
    vector = np.array([x, y])
    norm = np.linalg.norm(vector - vector_prev, ord=np.inf)

    #print(f'n: {n}, x0 = {x0}, y0 = {y0}')
    #print(f'x1 = {x}, y1 = {y}')
    #print(f'diff = {vector - vector_prev}, norm = {norm}')
    # print()

    while norm > (1-q())/q()*eps:
        n += 1
        #print(f'n: {n}')
        x_now, y_now = x, y
        y = phi(x, y)
        x = gamma(x, y)
        x_prev, y_prev = x_now, y_now

        vector_prev = np.array([x_prev, y_prev])
        vector = np.array([x, y])
        norm = np.linalg.norm(vector - vector_prev, ord=np.inf)

        #print(f'x = {x_now}, y = {y_now}')
        #print(f'x+1 = {x}, y+1 = {y}')
        #print(f'diff = {vector - vector_prev}, norm = {norm}')
        # print()

    print(f'Окончательное решение на итерации {n}: (x, y) = ({x}, {y})')
    print()

    return(x, y, n)


def newton(f, g, x0, y0, eps):
    print('Расчет методом Ньютона:')
    n = 0
    # Вектор переменных
    vector_prev = np.array([x0, y0]).reshape(-1, 1)
    # Вектор функций
    F = np.array([f(x0, y0), g(x0, y0)]).reshape(-1, 1)
    # Матрица Якоби
    W = np.array([
        [der_f_x(x0, y0), der_g_x(x0, y0)],
        [der_f_y(x0, y0), der_g_y(x0, y0)]
    ])

    # обобщение метода Ньютона на системы уравнений: X+1 = X - W^-1 * F
    vector = vector_prev - np.linalg.inv(W) @ F

    norm = np.linalg.norm(vector - vector_prev, ord=np.inf)

    #print(f'n: {n})')
    # print(f'vector0:\n{vector_prev}')
    # print(f'vector:\n{vector}')
    #print(f'diff:\n{vector - vector_prev},\nnorm = {norm}')
    # print()

    while norm > mu()*eps:
        n += 1
        #print(f'n: {n}')
        vector_now = vector

        # Разложим вектор на отдельные переменные для вычисления функций
        x, y = vector.reshape((1, 2))[0]

        # Вычисляем вектор функций
        F = np.array([f(x, y), g(x, y)]).reshape(-1, 1)
        #print(f'Значение вектора функций F = {F}')

        # Вычисляем якобиан
        W = np.array([
            [der_f_x(x, y), der_g_x(x, y)],
            [der_f_y(x, y), der_g_y(x, y)]
        ])

        #print('Матр. произведение обратной матрицы Якоби и вектора функций')
        #print(f'W^-1 * F:\n{np.linalg.inv(W) @ F}')

        vector = vector - np.linalg.inv(W) @ F

        vector_prev = vector_now
        norm = np.linalg.norm(vector - vector_prev, ord=np.inf)

        # print(f'\nvector_prev:\n{vector_prev}')
        # print(f'vector:\n{vector}')
        #print(f'diff:\n{vector - vector_prev},\nnorm = {norm}')
        # print()

    x, y = vector.reshape((1, 2))[0]
    print(f'Окончательное решение на итерации {n}: (x, y) = ({x}, {y})')
    print()
    return(x, y, n)


# Вычисления

with open("otchet.csv", "w") as file:
    file.write(", x, y, ||F(x y)||, N+1, x0, y0, q, m \n")
    [x, y, n] = iterations(phi, gamma, x0, y0, eps)
    file.write("Простая итерация, {}, {}, {}, {}, {}, {}, {},  \n".format(
        x, y, max(f(x, y), g(x, y)), n, x0, y0, q()))

    [x, y, n] = seidel(phi, gamma, x0, y0, eps)
    file.write("Метод Зейдаля, {}, {}, {}, {}, {}, {}, {},  \n".format(
        x, y, max(f(x, y), g(x, y)), n, x0, y0, q()))

    [x, y, n] = newton(f, g, x0, y0, eps)
    file.write("Метод Ньютона, {}, {}, {}, {}, {}, {}, , {} \n".format(
        x, y, max(f(x, y), g(x, y)), n, x0, y0, mu()))
file.close()
