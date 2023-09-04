#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define alpha_a 0             // коэффициент при производной в левом граничном условии
#define beta_a 1              // коэффициент при функции в левом граничном условии
#define L_a -1 / (2 * log(2)) // значение левого граничного условия

#define alpha_b 1 // коэффициент при производной в правом граничном условии
#define beta_b 0  // коэффициент при функции в правом граничном условии
#define L_b 0     // значение правого граничного условия

#define left 0.5 // левый конец интервала
#define right 1  // правый конец интервала
const int ns[] = {25, 50, 100, 200, 500, 1000, 2000, 4000, 8000};

long double b(long double x)
{ // коэффициент перед первой произвордной
    return 1.;
}

long double c(long double x)
{ // коэффициент перед функцией
    return -1 / x;
}

long double d(long double x)
{ // то, что после =
    return (x + 1) / x;
}

long double *create_grid(double a, double b, int n)
{
    long double *grid = new long double[n + 1];
    grid[0] = a;
    long double h = (b - a) / n;
    for (int i = 1; i < n; i++)
    {
        grid[i] = grid[i - 1] + h;
    }
    grid[n] = b;
    return grid;
}

long double norm_1(long double *v, int size)
{
    int n = size;
    long double ans = 0;
    for (int i = 0; i < n; i++)
    {
        ans += fabs(v[i]);
    }
    return ans;
}

long double norm_2(long double *v, int size)
{
    int n = size;
    long double ans = 0;
    for (int i = 0; i < n; i++)
    {
        ans += v[i] * v[i];
    }
    return ans;
}

long double norm_inf(long double *v, int size)
{
    int n = size;
    long double ans = fabs(v[0]);
    for (int i = 1; i < n; i++)
    {
        ans = fmax(ans, fabs(v[i]));
    }
    return ans;
}

long double *progonka(long double *A, long double *B, long double *C, long double *D, int size)
{
    int n = size;

    long double *y = new long double[n];
    long double *alpha = new long double[n];
    long double *beta = new long double[n];
    long double *x = new long double[n];

    // прямой ход
    y[0] = B[0];
    alpha[0] = -C[0] / y[0];
    beta[0] = D[0] / y[0];

    for (int i = 1; i < n - 1; i++)
    {
        y[i] = B[i] + A[i] * alpha[i - 1];
        alpha[i] = -C[i] / y[i];
        beta[i] = (D[i] - A[i] * beta[i - 1]) / y[i];
    }

    y[n - 1] = B[n - 1] + A[n - 1] * alpha[n - 2];
    beta[n - 1] = (D[n - 1] - A[n - 1] * beta[n - 2]) / y[n - 1];

    // обратный ход

    x[n - 1] = beta[n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    return x;
}

long double *u_1st_method(double low, double high, int n)
{ // метод несимметричной производной
    long double *A = new long double[n + 1];
    long double *B = new long double[n + 1];
    long double *C = new long double[n + 1];
    long double *D = new long double[n + 1];

    long double h = (high - low) / n;
    long double *x = create_grid(low, high, n);

    for (int i = 1; i < n; i++)
    {
        D[i] = d(x[i]);
    }

    B[0] = -2 * alpha_a + 2 * beta_a * h - 2 * b(x[1]) * h * alpha_a + beta_a * b(x[1]) * h * h;
    C[0] = 2 * alpha_a + 2 * b(x[1]) * h * alpha_a + c(x[1]) * h * h * alpha_a;
    D[0] = 2 * L_a * h + L_a * b(x[1]) * h * h + d(x[1]) * h * h * alpha_a;

    for (int i = 1; i < n; i++)
    {
        A[i] = -(1 / (h * h) - b(x[i]) / (2 * h));
        B[i] = 2 / (h * h) - c(x[i]);
        C[i] = -(1 / (h * h) + b(x[i]) / (2 * h));
    }

    B[n] = -2 * alpha_b - 2 * h * beta_b - b(x[n - 1]) * h - h * h * b(x[n - 1]) * beta_b;
    A[n] = 2 * alpha_b - 2 * b(x[n - 1]) * h + alpha_b * h * h * c(x[n - 1]);
    D[n] = -2 * L_b * h - b(x[n - 1]) * h * h * L_b + alpha_b * h * h * d(x[n - 1]);

    /*
    for (int i = 0; i < n + 1; i++)
    {
        cout << A[i] << " " << B[i] << " " << C[i] << " " << D[i] << endl;
    }
    */

    long double *u = progonka(A, B, C, D, n);

    return u;
}

long double *v_2nd_method(double low, double high, int n)
{ // метод фиктивного узла
    long double *A = new long double[n + 1];
    long double *B = new long double[n + 1];
    long double *C = new long double[n + 1];
    long double *D = new long double[n + 1];

    long double h = (high - low) / n;
    long double *x = create_grid(low, high, n);

    for (int i = 1; i < n; i++)
    {
        D[i] = d(x[i]);
    }

    B[0] = (2 / (h * h) - c(x[0])) * alpha_a / 2 - (1 / h - b(x[0]) / 2);
    C[0] = -alpha_a / (h * h);
    D[0] = -(L_a * (1 / h - b(x[0]) / 2) + d(x[0]) * alpha_a / 2);

    for (int i = 1; i < n; i++)
    {
        A[i] = -(1 / (h * h) - b(x[i]) / (2 * h));
        B[i] = 2 / (h * h) - c(x[i]);
        C[i] = -(1 / (h * h) + b(x[i]) / (2 * h));
    }

    A[n] = alpha_b / (h * h);
    B[n] = -(alpha_b * (1 / (h * h) - c(x[n]) / 2) + beta_b * (1 / h + b(x[n]) / (2 * h)));
    D[n] = -((1 / h + b(x[n]) / (2 * h)) * L_b - d(x[n]) * alpha_b / 2);

    /*
    for (int i = 0; i < n + 1; i++)
    {
        cout << A[i] << " " << B[i] << " " << C[i] << " " << D[i] << endl;
    }
    */

    long double *v = progonka(A, B, C, D, n);

    return v;
}

long double *subs(long double *a, long double *b, int n)
{
    long double *ans = new long double[n];
    for (int i = 0; i < n; i++)
    {
        ans[i] = a[i] - b[i];
    }
    return ans;
}

int main()
{
    FILE *pFile;
    pFile = fopen("otchet.tex", "w");
    fprintf(pFile, "%s \n", "\\documentclass[a4paper]{article}");
    fprintf(pFile, "%s \n", "\\usepackage[12pt]{extsizes}");
    fprintf(pFile, "%s \n", "\\usepackage[utf8]{inputenc}");
    fprintf(pFile, "%s \n", "\\usepackage[russian]{babel}");
    fprintf(pFile, "%s \n", "\\usepackage{setspace,amsmath}");
    fprintf(pFile, "%s \n", "\\usepackage[left=20mm, top=20mm, right=20mm, bottom=20mm, nohead, footskip=10mm, landscape]{geometry}");
    fprintf(pFile, "%s \n", "\\usepackage{indentfirst}");
    fprintf(pFile, "%s \n", "\\usepackage{multirow}");
    fprintf(pFile, "%s \n", "\\usepackage{fourier}");
    fprintf(pFile, "%s \n", "\\usepackage{array}");
    fprintf(pFile, "%s \n \n", "\\usepackage{makecell}");

    fprintf(pFile, "%s \n", "\\newcolumntype{P}[1]{>{\\centering\\arraybackslash}p{#1}} \n");
    fprintf(pFile, "%s \n", "\\newcolumntype{M}[1]{>{\\centering\\arraybackslash}m{#1}} \n");

    fprintf(pFile, "%s \n", "\\begin{document}");
    fprintf(pFile, "%s \n", "\\renewcommand{\\arraystretch}{2}");
    fprintf(pFile, "%s \n", "\\begin{center}");
    fprintf(pFile, "%s \n", "{\\Large Решение краевой задачи ОДУ сеточным методом с Использовнием метода прогонки.}");
    fprintf(pFile, "%s \n", "\\bigskip");
    fprintf(pFile, "%s \n", "\\end{center}");
    fprintf(pFile, "%s \n \n", "Выполнил: Капелюшников Андрей Сергеевич");
    fprintf(pFile, "%s \n", "\\begin{equation*}");
    fprintf(pFile, "%s \n", "\\begin{cases}");
    fprintf(pFile, "%s \n", "u\'\' (x) + u\' (x) - \\frac{1}{x}u=\\frac{1+x}{x}");
    fprintf(pFile, "%s \n", "\\\\");
    fprintf(pFile, "%s \n", "0,5 < x <1");
    fprintf(pFile, "%s \n", "\\\\");
    fprintf(pFile, "%s \n", "u(0,5)=\\frac{-1}{2 \\ln2}");
    fprintf(pFile, "%s \n", "\\\\");
    fprintf(pFile, "%s \n", "u\'(1)=0");
    fprintf(pFile, "%s \n", "\\end{cases}");
    fprintf(pFile, "%s \n", "\\end{equation*}");

    fprintf(pFile, "%s \n", "\\begin{table}[h!]");
    fprintf(pFile, "%s \n", "\\begin{tabular}{ccccccccc}");
    fprintf(pFile, "%s \n", "\\hline");

    int n_size = sizeof(ns) / sizeof(*ns);

    fprintf(pFile, "$n_1/n_2$ ");
    for (int i = 0; i < n_size - 1; i++)
    {
        fprintf(pFile, " & %d / %d ", ns[i], ns[i + 1]);
    }
    fprintf(pFile, "\\\\ \n");
    fprintf(pFile, "%s \n", "\\hline");

    fprintf(pFile, " $||\\cdot||_{\\infty}$ ");

    for (int i = 0; i < n_size - 1; i++)
    {
        long double *u1 = u_1st_method(left, right, ns[i]);
        long double *v1 = v_2nd_method(left, right, ns[i]);

        long double *u2 = u_1st_method(left, right, ns[i + 1]);
        long double *v2 = v_2nd_method(left, right, ns[i + 1]);

        long double ans = norm_inf(subs(u1, v2, ns[i]), ns[i]) / norm_inf(subs(u2, v2, ns[i + 1]), ns[i + 1]);

        fprintf(pFile, " & %g", ans);
    }
    fprintf(pFile, "\\\\ \n");
    fprintf(pFile, "%s \n", "\\hline");

    fprintf(pFile, " $||\\cdot||_{2}$ ");
    for (int i = 0; i < n_size - 1; i++)
    {
        long double *u1 = u_1st_method(left, right, ns[i]);
        long double *v1 = v_2nd_method(left, right, ns[i]);

        long double *u2 = u_1st_method(left, right, ns[i + 1]);
        long double *v2 = v_2nd_method(left, right, ns[i + 1]);

        long double ans = norm_2(subs(u1, v1, ns[i]), ns[i]) / norm_2(subs(u2, v2, ns[i + 1]), ns[i + 1]);

        fprintf(pFile, " & %g", ans);
    }

    fprintf(pFile, "\\\\ \n");
    fprintf(pFile, "%s \n", "\\hline");

    fprintf(pFile, " $||\\cdot||_{1}$ ");
    for (int i = 0; i < n_size - 1; i++)
    {
        long double *u1 = u_1st_method(left, right, ns[i]);
        long double *v1 = v_2nd_method(left, right, ns[i]);

        long double *u2 = u_1st_method(left, right, ns[i + 1]);
        long double *v2 = v_2nd_method(left, right, ns[i + 1]);

        long double ans = norm_1(subs(u1, v1, ns[i]), ns[i]) / norm_1(subs(u2, v2, ns[i + 1]), ns[i + 1]);
        fprintf(pFile, " & %g", ans);
    }
    fprintf(pFile, "\\\\ \n");
    fprintf(pFile, "%s \n", "\\hline");

    fprintf(pFile, "%s \n", "\\end{tabular}");
    fprintf(pFile, "%s \n", "\\end{table}");

    /*
    fprintf(pFile, "%s \n", );

    fprintf(pFile, "%s \n", );
    fprintf(pFile, "%s \n", );
    fprintf(pFile, "%s \n", );
    fprintf(pFile, "%s \n", );
    fprintf(pFile, "%s \n", );
    fprintf(pFile, "%s \n", );
    */

    /*
    // Проверка работостпсобности проконки
    double A[3] = {0, 5, 1};
    double B[3] = {2, 4, -3};
    double C[3] = {-1, 2, 0};
    double D[3] = {3, 6, 2};

    double *x = progonka(A, B, C, D, 3);

    cout << x[0] << " " << x[1] << " " << x[2];

    // верный ответ 1.48837 -0.0232558 -0.674419
    // спасибо https://pro-prof.com/forums/topic/sweep-method-for-solving-systems-of-linear-algebraic-equations?ysclid=lcp09zvrpb351052544
    */
    fprintf(pFile, "%s \n", "\\end{document}");
    fclose(pFile);

    long double *u = u_1st_method(left, right, ns[6]);
    long double *x = create_grid(left, right, ns[6]);

    FILE *csvFile;
    csvFile = fopen("n2000.csv", "w");

    for (int i = 0; i < ns[6]; i++)
    {
        fprintf(csvFile, "%e, %e \n", (double)x[i], (double)u[i]);
        // cout << x[i] << ", " << u[i] << endl;
    }

    fclose(csvFile);
    return 0;
}