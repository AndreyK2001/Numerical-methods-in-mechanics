from math import *
import numpy as np

# интервал
a_0 = 0
b_0 = 1

# u(0)-2u'(0)=0
mu_a = 1
lambda_a = -2
psi_a = 0

# u(1)=-1/sqrt(2)

mu_b = 1
lambda_b = 0
psi_b = -1/sqrt(2)

n_vals = [25, 50, 100, 200, 500, 1000, 2000, 4000, 8000]

# u''(x) + p(x)u'+q(x)u=g(x)


def p(x):
    return -sin(x)


def q(x):
    return 2./(x+1)**2


def g(x):
    return 9./(2*(x+1)**(2./3))


def grid(a, b, n):
    h = (b - a) / n
    return np.array([a+h*i for i in range(n+1)])


def solver(a, b, c, d, n):
    alpha = np.zeros(n)
    beta = np.zeros(n)
    y = np.zeros(n)

    y[0] = b[0]
    alpha[0] = -c[0] / y[0]
    beta[0] = d[0] / y[0]

    for i in range(1, n-1):
        y[i] = b[i] + a[i] * alpha[i - 1]
        alpha[i] = -c[i] / y[i]
        beta[i] = (d[i] - a[i] * beta[i - 1]) / y[i]

    y[n - 1] = b[n - 1] + a[n - 1] * alpha[n - 2]
    beta[n - 1] = (d[n - 1] - a[n - 1] * beta[n - 2]) / y[n - 1]

    x = np.zeros(n)
    x[n - 1] = beta[n - 1]

    for i in range(n-2, -1, -1):
        x[i] = alpha[i] * x[i + 1] + beta[i]

    return x


def u(a_, b_, n):
    A = np.zeros(n+1)
    B = np.zeros(n+1)
    C = np.zeros(n+1)

    h = (b_ - a_) / n
    x = grid(a_, b_, n)

    D = np.array([g(x[i]) for i in range(n+1)])

    B[0] = -2 * lambda_a + 2 * mu_a * h - 2 * \
        p(x[1]) * h * lambda_a + mu_a * p(x[1]) * h**2
    C[0] = 2 * lambda_a + 2 * p(x[1]) * h * lambda_a + \
        q(x[1]) * h ** 2 * lambda_a
    D[0] = 2 * psi_a * h + psi_a * p(x[1]) * h**2 + g(x[1]) * h ** 2 * lambda_a

    for i in range(1, n):
        A[i] = -(1 / (h * h) - p(x[i]) / (2 * h))
        B[i] = 2 / (h * h) - q(x[i])
        C[i] = -(1 / (h * h) + p(x[i]) / (2 * h))

    B[n] = -2 * lambda_b - 2 * h * mu_b - \
        p(x[n - 1]) * h - h**2 * p(x[n - 1]) * mu_b
    A[n] = 2 * lambda_b - 2 * p(x[n - 1]) * h + lambda_b * h ** 2 * q(x[n - 1])
    D[n] = -2 * psi_b * h - p(x[n - 1]) * h ** 2 * \
        psi_b + lambda_b * h * h * g(x[n - 1])

    return solver(A, B, C, D, n)


def v(a_, b_, n):
    A = np.zeros(n+1)
    B = np.zeros(n+1)
    C = np.zeros(n+1)

    h = (b_ - a_) / n
    x = grid(a_, b_, n)

    D = np.array([g(x[i]) for i in range(n+1)])

    B[0] = (2 / (h * h) - q(x[0])) * lambda_a / 2 - (1 / h - p(x[0]) / 2)
    C[0] = -lambda_a / (h * h)
    D[0] = -(psi_a * (1 / h - p(x[0]) / 2) + g(x[0]) * lambda_a / 2)

    for i in range(1, n):
        A[i] = -(1 / (h * h) - p(x[i]) / (2 * h))
        B[i] = 2 / (h * h) - q(x[i])
        C[i] = -(1 / (h * h) + p(x[i]) / (2 * h))

    A[n] = lambda_b / (h * h)
    B[n] = -(lambda_b * (1 / (h * h) - q(x[n]) / 2) +
             mu_b * (1 / h + p(x[n]) / (2 * h)))
    D[n] = -((1 / h + p(x[n]) / (2 * h)) * psi_b - g(x[n]) * lambda_b / 2)

    return solver(A, B, C, D, n)


with open("table.csv", "w") as f:

    f.write("n1/n2")
    for i in range(len(n_vals)-1):
        f.write(", {}/{}".format(n_vals[i], n_vals[i+1]))
    f.write("\n")

    f.write("||•||")
    for i in range(len(n_vals)-1):
        u_ = u(a_0, b_0, n_vals[i])
        v_ = v(a_0, b_0, n_vals[i])

        res1 = np.linalg.norm(u_-v_, ord=np.inf)

        u_ = u(a_0, b_0, n_vals[i+1])
        v_ = v(a_0, b_0, n_vals[i+1])

        res2 = np.linalg.norm(u_-v_, ord=np.inf)

        f.write(", {:.4f}".format(res1 / res2))
    f.write("\n")

    f.write("||•||1")
    for i in range(len(n_vals)-1):
        u_ = u(a_0, b_0, n_vals[i])
        v_ = v(a_0, b_0, n_vals[i])

        res1 = np.linalg.norm(u_-v_, ord=1)

        u_ = u(a_0, b_0, n_vals[i+1])
        v_ = v(a_0, b_0, n_vals[i+1])

        res2 = np.linalg.norm(u_-v_, ord=1)

        f.write(", {:.4f}".format(res1 / res2))
    f.write("\n")

    f.write("||•||2")
    for i in range(len(n_vals)-1):
        u_ = u(a_0, b_0, n_vals[i])
        v_ = v(a_0, b_0, n_vals[i])

        res1 = np.linalg.norm(u_-v_, ord=2)

        u_ = u(a_0, b_0, n_vals[i+1])
        v_ = v(a_0, b_0, n_vals[i+1])

        res2 = np.linalg.norm(u_-v_, ord=2)

        f.write(", {:.4f}".format(res1 / res2))
    f.write("\n")

    f.close()

with open("u.csv", "w") as f:

    n = 2000
    x = grid(a_0, b_0, n)
    u_ = u(a_0, b_0, n)

    for i in range(n):
        f.write("{}, {} \n".format(x[i], u_[i] + psi_b))
    
    f.close()
