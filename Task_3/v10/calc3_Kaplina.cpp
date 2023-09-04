#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define lambda_a 1
#define mu_a 0
#define psi_a 0

#define lambda_b 0
#define mu_b 1
#define psi_b 0.5 * sin(0.5)

#define a_0 0
#define b_0 0.5
const int n_values[] = {25, 50, 100, 200, 500, 1000, 2000, 4000, 8000};

// u''(x) + p(x)u'+q(x)u=g(x)

double p(double x)
{
    return 2. * x;
}

double q(double x)
{
    return -sin(x);
}

double g(double x)
{
    return 2. * (x * x + 1) * cos(M_PI * x);
}

double *grid(double a, double b, int n)
{ //создание сетки
    double *grid = new double[n + 1];
    grid[0] = a;
    double h = (b - a) / n;
    for (int i = 1; i < n; i++)
        grid[i] = grid[i - 1] + h;
    grid[n] = b;
    return grid;
}

double *solver(double *a, double *b, double *c, double *d, int n)
{ //метод прогонки
    double *alpha = new double[n];
    double *beta = new double[n];
    double *y = new double[n];

    // прямой ход
    y[0] = b[0];
    alpha[0] = -c[0] / y[0];
    beta[0] = d[0] / y[0];
    for (int i = 1; i < n - 1; i++)
    {
        y[i] = b[i] + a[i] * alpha[i - 1];
        alpha[i] = -c[i] / y[i];
        beta[i] = (d[i] - a[i] * beta[i - 1]) / y[i];
    }
    y[n - 1] = b[n - 1] + a[n - 1] * alpha[n - 2];
    beta[n - 1] = (d[n - 1] - a[n - 1] * beta[n - 2]) / y[n - 1];

    // обратный ход
    double *x = new double[n];
    x[n - 1] = beta[n - 1];

    for (int i = n - 2; i >= 0; i--)
        x[i] = alpha[i] * x[i + 1] + beta[i];

    return x;
}

double norma_l1(double *v, int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
        norm += fabs(v[i]);
    return norm;
}

double norma_l2(double *v, int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
        norm += v[i] * v[i];
    return norm;
}

double norm_uniform(double *v, int n)
{ //равномерная норма
    double norm = fabs(v[0]);
    for (int i = 1; i < n; i++)
        norm = fmax(norm, fabs(v[i]));
    return norm;
}

double *decriment(double *a, double *b, int n)
{
    double *t = new double[n];
    for (int i = 0; i < n; i++)
        t[i] = a[i] - b[i];

    return t;
}

double *u(double a_, double b_, int n)
{ // метод несимметричной производной

    double *A = new double[n + 1];
    double *B = new double[n + 1];
    double *C = new double[n + 1];
    double *D = new double[n + 1];

    double h = (b_ - a_) / n;
    double *x = grid(a_, b_, n);

    for (int i = 1; i < n; i++)
        D[i] = g(x[i]);

    B[0] = -2 * lambda_a + 2 * mu_a * h - 2 * p(x[1]) * h * lambda_a + mu_a * p(x[1]) * h * h;
    C[0] = 2 * lambda_a + 2 * p(x[1]) * h * lambda_a + q(x[1]) * h * h * lambda_a;
    D[0] = 2 * psi_a * h + psi_a * p(x[1]) * h * h + g(x[1]) * h * h * lambda_a;

    for (int i = 1; i < n; i++)
    {
        A[i] = -(1 / (h * h) - p(x[i]) / (2 * h));
        B[i] = 2 / (h * h) - q(x[i]);
        C[i] = -(1 / (h * h) + p(x[i]) / (2 * h));
    }

    B[n] = -2 * lambda_b - 2 * h * mu_b - p(x[n - 1]) * h - h * h * p(x[n - 1]) * mu_b;
    A[n] = 2 * lambda_b - 2 * p(x[n - 1]) * h + lambda_b * h * h * q(x[n - 1]);
    D[n] = -2 * psi_b * h - p(x[n - 1]) * h * h * psi_b + lambda_b * h * h * g(x[n - 1]);

    return solver(A, B, C, D, n);
}

double *v(double a_, double b_, int n)
{ // метод фиктивного узла

    double *A = new double[n + 1];
    double *B = new double[n + 1];
    double *C = new double[n + 1];
    double *D = new double[n + 1];

    double h = (b_ - a_) / n;
    double *x = grid(a_, b_, n);

    for (int i = 1; i < n; i++)
        D[i] = g(x[i]);

    B[0] = (2 / (h * h) - q(x[0])) * lambda_a / 2 - (1 / h - p(x[0]) / 2);
    C[0] = -lambda_a / (h * h);
    D[0] = -(psi_a * (1 / h - p(x[0]) / 2) + g(x[0]) * lambda_a / 2);

    for (int i = 1; i < n; i++)
    {
        A[i] = -(1 / (h * h) - p(x[i]) / (2 * h));
        B[i] = 2 / (h * h) - q(x[i]);
        C[i] = -(1 / (h * h) + p(x[i]) / (2 * h));
    }

    A[n] = lambda_b / (h * h);
    B[n] = -(lambda_b * (1 / (h * h) - q(x[n]) / 2) + mu_b * (1 / h + p(x[n]) / (2 * h)));
    D[n] = -((1 / h + p(x[n]) / (2 * h)) * psi_b - g(x[n]) * lambda_b / 2);

    return solver(A, B, C, D, n);
}

int main()
{
    ofstream fout; // Открываем поток
    fout.open("calc3_Kaplina.csv");

    int num_vals = sizeof(n_values) / sizeof(*n_values);
    fout << "n1/n2";

    for (int i = 0; i < num_vals - 1; i++)
        fout << ", " << n_values[i] << "/" << n_values[i + 1];

    fout << endl;
    fout << "||•||";
    for (int i = 0; i < num_vals - 1; i++)
    {
        double *u_ = u(a_0, b_0, n_values[i]);
        double *v_ = v(a_0, b_0, n_values[i]);
        double calc1 = norm_uniform(decriment(u_, v_, n_values[i]), n_values[i]);

        free(u_);
        free(v_);

        u_ = u(a_0, b_0, n_values[i + 1]);
        v_ = v(a_0, b_0, n_values[i + 1]);
        double calc2 = norm_uniform(decriment(u_, v_, n_values[i + 1]), n_values[i + 1]);

        fout << ", " << calc1 / calc2;
    }

    fout << endl;

    fout << "||•||1";
    for (int i = 0; i < num_vals - 1; i++)
    {
        double *u_ = u(a_0, b_0, n_values[i]);
        double *v_ = v(a_0, b_0, n_values[i]);
        double calc1 = norma_l1(decriment(u_, v_, n_values[i]), n_values[i]);

        free(u_);
        free(v_);

        u_ = u(a_0, b_0, n_values[i + 1]);
        v_ = v(a_0, b_0, n_values[i + 1]);
        double calc2 = norma_l1(decriment(u_, v_, n_values[i + 1]), n_values[i + 1]);

        fout << ", " << calc1 / calc2;
    }

    fout << endl;

    fout << "||•||2";

    for (int i = 0; i < num_vals - 1; i++)
    {
        double *u_ = u(a_0, b_0, n_values[i]);
        double *v_ = v(a_0, b_0, n_values[i]);
        double calc1 = norma_l2(decriment(u_, v_, n_values[i]), n_values[i]);

        free(u_);
        free(v_);

        u_ = u(a_0, b_0, n_values[i + 1]);
        v_ = v(a_0, b_0, n_values[i + 1]);
        double calc2 = norma_l2(decriment(u_, v_, n_values[i + 1]), n_values[i + 1]);

        fout << ", " << calc1 / calc2;
    }

    fout << endl;
    fout.close();

    fout.open("calc3_Kaplina_u_vs_x.csv");

    int num_points = 2000;
    double *x = grid(a_0, b_0, num_points);
    double *u_ = u(a_0, b_0, num_points);

    for (int i = 0; i < num_points; i++)
        fout << x[i] << ", " << u_[i] + psi_b << endl;

    fout.close();
    return 0;
}
