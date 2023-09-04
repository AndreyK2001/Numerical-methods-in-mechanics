#include <math.h>  // Библиотека станлартных математических функций
#include <fstream> // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define a_0 0.2
#define b_0 1.

const int N[] = {25, 50, 100, 200, 500, 1000, 2000, 4000, 8000};

#define mu_a 1
#define lambda_a 1 // u(0.2)+u'(0.2)=1
#define psi_a 1

#define mu_b 0
#define lambda_b 1 // u(1)=0.5
#define psi_b 0.5

double *progonka(double *a, double *b, double *c, double *d, int n)
{
    double *Alpha = new double[n]; // коэффициенты
    double *Beta = new double[n];
    double *u = new double[n];

    u[0] = b[0];
    Alpha[0] = -c[0] / u[0];
    Beta[0] = d[0] / u[0];
    for (int i = 1; i < n - 1; i++)
    {
        u[i] = b[i] + a[i] * Alpha[i - 1];
        Alpha[i] = -c[i] / u[i];
        Beta[i] = (d[i] - a[i] * Beta[i - 1]) / u[i];
    }
    u[n - 1] = b[n - 1] + a[n - 1] * Alpha[n - 2];
    Beta[n - 1] = (d[n - 1] - a[n - 1] * Beta[n - 2]) / u[n - 1];

    double *x = new double[n]; // решения
    x[n - 1] = Beta[n - 1];
    for (int i = n - 2; i >= 0; i--)
        x[i] = Alpha[i] * x[i + 1] + Beta[i];

    free(Alpha);
    free(Beta);
    free(u);
    return x;
}

// u''(x) + p(x)u'+q(x)u = g(x)

double p(double x) { return -sqrt(2 / x); }

double q(double x) { return 4 / (x * x + 2); }

double g(double x) { return 2 / (x * x * x - 7); }

double *grid(double a, double b, int n)
{
    double *grid = new double[n + 1];
    grid[0] = a;
    double h = (b - a) / n;
    for (int i = 1; i < n; i++)
    {
        grid[i] = grid[i - 1] + h;
    }

    grid[n] = b;
    return grid;
}

double *u(double a, double b, int n)
{ // метод несимметричной производной

    double *A = new double[n + 1];
    double *B = new double[n + 1];
    double *C = new double[n + 1];
    double *D = new double[n + 1];

    double h = (b - a) / n;
    double *x = grid(a, b, n);

    for (int i = 1; i < n; i++)
        D[i] = g(x[i]);

    B[0] = -2 * mu_a + 2 * lambda_a * h - 2 * p(x[1]) * h * mu_a + lambda_a * p(x[1]) * h * h;
    C[0] = 2 * mu_a + 2 * p(x[1]) * h * mu_a + q(x[1]) * h * h * mu_a;
    D[0] = 2 * psi_a * h + psi_a * p(x[1]) * h * h + g(x[1]) * h * h * mu_a;

    for (int i = 1; i < n; i++)
    {
        A[i] = -(1 / (h * h) - p(x[i]) / (2 * h));
        B[i] = 2 / (h * h) - q(x[i]);
        C[i] = -(1 / (h * h) + p(x[i]) / (2 * h));
    }
    A[n] = 2 * mu_b - 2 * p(x[n - 1]) * h + mu_b * h * h * q(x[n - 1]);
    B[n] = -2 * mu_b - 2 * h * lambda_b - p(x[n - 1]) * h - h * h * p(x[n - 1]) * lambda_b;
    D[n] = -2 * psi_b * h - p(x[n - 1]) * h * h * psi_b + mu_b * h * h * g(x[n - 1]);

    return progonka(A, B, C, D, n);
}

double *v(double a, double b, int n)
{ // метод фиктивного узла

    double *A = new double[n + 1];
    double *B = new double[n + 1];
    double *C = new double[n + 1];
    double *D = new double[n + 1];

    double h = (b - a) / n;
    double *x = grid(a, b, n);

    for (int i = 1; i < n; i++)
        D[i] = g(x[i]);

    B[0] = (2 / (h * h) - q(x[0])) * mu_a / 2 - (1 / h - p(x[0]) / 2);
    C[0] = -mu_a / (h * h);
    D[0] = -(psi_a * (1 / h - p(x[0]) / 2) + g(x[0]) * mu_a / 2);

    for (int i = 1; i < n; i++)
    {
        A[i] = -(1 / (h * h) - p(x[i]) / (2 * h));
        B[i] = 2 / (h * h) - q(x[i]);
        C[i] = -(1 / (h * h) + p(x[i]) / (2 * h));
    }

    A[n] = mu_b / (h * h);
    B[n] = -(mu_b * (1 / (h * h) - q(x[n]) / 2) + lambda_b * (1 / h + p(x[n]) / (2 * h)));
    D[n] = -((1 / h + p(x[n]) / (2 * h)) * psi_b - g(x[n]) * mu_b / 2);

    return progonka(A, B, C, D, n);
}

double *decriment(double *a, double *b, int n)
{
    double *t = new double[n];
    for (int i = 0; i < n; i++)
        t[i] = a[i] - b[i];

    return t;
}

double l_1(double *v, int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
        norm += fabs(v[i]);
    return norm;
}

double l_2(double *v, int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
        norm += v[i] * v[i];
    return sqrt(norm);
}

double l_inf(double *v, int n)
{ // равномерная норма
    double norm = fabs(v[0]);
    for (int i = 1; i < n; i++)
        norm = fmax(norm, fabs(v[i]));
    return norm;
}

int main()
{
    FILE *pFile;
    pFile = fopen("eq_diff.csv", "w");

    int size_ = sizeof(N) / sizeof(*N);

    fprintf(pFile, "%s", "n1/n2");

    for (int i = 0; i < size_ - 1; i++)
    {
        fprintf(pFile, ",%d / %d", N[i], N[i + 1]);
    }

    fprintf(pFile, "\n %s", "||•||");

    for (int i = 0; i < size_ - 1; i++)
    {
        double *u_ = u(a_0, b_0, N[i]);
        double *v_ = v(a_0, b_0, N[i]);
        double calc1 = l_inf(decriment(u_, v_, N[i]), N[i]);
        free(u_);
        free(v_);
        u_ = u(a_0, b_0, N[i + 1]);
        v_ = v(a_0, b_0, N[i + 1]);
        double calc2 = l_inf(decriment(u_, v_, N[i + 1]), N[i + 1]);

        fprintf(pFile, ", %g", calc1 / calc2);
    }

    fprintf(pFile, "\n %s", "||•||1");

    for (int i = 0; i < size_ - 1; i++)
    {
        double *u_ = u(a_0, b_0, N[i]);
        double *v_ = v(a_0, b_0, N[i]);
        double calc1 = l_1(decriment(u_, v_, N[i]), N[i]);
        free(u_);
        free(v_);
        u_ = u(a_0, b_0, N[i + 1]);
        v_ = v(a_0, b_0, N[i + 1]);
        double calc2 = l_1(decriment(u_, v_, N[i + 1]), N[i + 1]);

        fprintf(pFile, ", %g", calc1 / calc2);
    }

    fprintf(pFile, "\n %s", "||•||2");

    for (int i = 0; i < size_ - 1; i++)
    {
        double *u_ = u(a_0, b_0, N[i]);
        double *v_ = v(a_0, b_0, N[i]);
        double calc1 = l_2(decriment(u_, v_, N[i]), N[i]);
        free(u_);
        free(v_);
        u_ = u(a_0, b_0, N[i + 1]);
        v_ = v(a_0, b_0, N[i + 1]);
        double calc2 = l_2(decriment(u_, v_, N[i + 1]), N[i + 1]);

        fprintf(pFile, ", %g ", calc1 / calc2);
    }

    fclose(pFile);

    FILE *csv;
    csv = fopen("u(x).csv", "w");

    int num_points = 2000;
    double *x = grid(a_0, b_0, num_points);
    double *u_ = u(a_0, b_0, num_points);

    for (int i = 0; i < num_points; i++)
    {
        fprintf(csv, "%g, %g \n", x[i], u_[i] / 1.622 + psi_b);
    }

    fclose(csv);
    return 0;
}
