#include <math.h>  // Библиотека станлартных математических функций
#include <fstream> // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001 // Требуемая точность вычисления
#define infty 1 / (eps * eps)

double si(double x)
{
    return sin(log(sqrt(1 + pow(1 - exp(-fabs(x)), 2))));
}

double co(double x) { return cos(log(sqrt(1 + pow(1 - exp(-fabs(x)), 2)))); }

double f(double x, double a) { return pow(si(x), 2) + a * log(x); }

double sign(double x)
{
    // знак х
    if (x == 0)
        return 0;
    else if (x > 0)
        return 1;
    else
        return -1;
}

double delta(double x)
{ // дельта-функция Дирака
    if (x == 0)
        return infty;
    else
        return 0;
}

double dfdx(double x, double a)
{ // первая производная
    return a / x + 2 * (1 - exp(-fabs(x))) * exp(-fabs(x)) * si(x) * co(x) * sign(x) / (1 + pow(1 - exp(-fabs(x)), 2));
}

double d2fdx2(double x, double a)
{ // вторая производная
    return -a / (x * x) -
           2 * pow(1 - exp(-fabs(x)), 2) * exp(-2 * fabs(x)) * pow(si(x), 2) / pow(1 + pow(1 - exp(-fabs(x)), 2), 2) -
           4 * pow(1 - exp(-fabs(x)), 2) * exp(-2 * fabs(x)) * si(x) * co(x) / pow(1 + pow(1 - exp(-fabs(x)), 2), 2) +
           2 * pow(1 - exp(-fabs(x)), 2) * exp(-2 * fabs(x)) * pow(co(x), 2) / pow(1 + pow(1 - exp(-fabs(x)), 2), 2) +
           4 * (1 - exp(-fabs(x))) * exp(fabs(x)) * si(x) * co(x) * delta(x) / (1 + pow(1 - exp(-fabs(x)), 2)) -
           2 * (1 - exp(-fabs(x))) * exp(fabs(x)) * si(x) * co(x) / (1 + pow(1 - exp(-fabs(x)), 2)) +
           2 * exp(-2 * fabs(x)) * si(x) * co(x) / (1 + pow(1 - exp(-fabs(x)), 2));
}

void Newtone(FILE *pFile, double (*f)(double, double), double x0, double alpha, double a, double b) // Метод Ньютона (метод касательных)
{
    int i = 0; // Число витков цикла n, счетчик i
    double x_2 = x0, x_1 = x0;

    double m1 = min(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha)));
    double M1 = max(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha)));

    double m = (eps + m1) / 2;
    double M = M1 + eps;

    while ((fabs(x_2 - x_1) > eps * m / (M - m)) || i == 0)
    {
        x_1 = x_2;
        i++;
        x_2 = x_1 - f(x_1, alpha) / dfdx(x_1, alpha);
    }

    fprintf(pFile, "%g, %g, %g, %d, %g, %g, - \n", x_2, f(x_2, alpha), x0, i, M, m);
}

void Half_div(FILE *pFile, double (*f)(double, double), double x0, double alpha, double a, double b)
{                                         // Метод деления отрезка пополам (метод половинного деления)
    int n, i;                             // Число витков цикла n, счетчик i
    double x = (a + b) / 2, l = a, r = b; // Вычисляем середину начального отрезка

    n = int(log((r - l) / eps) / log(2)) + 1;

    for (i = 1; i <= n && fabs(f(x, alpha)) > eps; i++)
    {
        if (f(x, alpha) * f(l, alpha) < 0)
            r = x;
        else
            l = x;

        x = (l + r) / 2; // Вычисляем середину нового отрезка
    }

    fprintf(pFile, "%g, %g, %g, %d, -, -, - \n", x, f(x, alpha), x0, i);
}

double phi(double (*f)(double, double), double x, double alpha, double a, double b)
{
    double M = max(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha))) + 3 * eps;
    return x - sign(dfdx(x, alpha)) * f(x, alpha) / M;
}

double phi_dfdx(double (*f)(double, double), double x, double alpha, double a, double b)
{
    double M = max(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha))) + 3 * eps;
    return 1 - delta(dfdx(x, alpha)) * d2fdx2(x, alpha) * f(x, alpha) / M - sign(dfdx(x, alpha)) * dfdx(x, alpha) / M;
}

void Simple_iteration(FILE *pFile, double (*f)(double, double), double x0, double alpha, double a, double b) // Метод простой итерации
{
    double x_2 = x0, x_1 = x0;
    int i = 0; // Число витков цикла n, счетчик i

    double M1 = max(fabs(phi_dfdx(f, a, alpha, a, b)), fabs(phi_dfdx(f, b, alpha, a, b))) + eps;
    double q = (1 - eps * 2 + M1) / 2;

    while (fabs(x_2 - x_1) > eps * (1 - q) / q / 10 || i == 0)
    {
        x_1 = x_2;
        i++;
        x_2 = phi(f, x_1, alpha, a, b);
    }

    fprintf(pFile, "%g, %g, %g, %d, -, -, %g \n", x_2, f(x_2, alpha), x0, i, q);
}

void Chord(FILE *pFile, double (*f)(double, double), double x0, double alpha, double a, double b)
{ // Метод хорд
    double x_1 = x0, x_2 = x0;

    int i = 0; // счётчик итераций алгоритма

    double m1 = min(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha)));
    double M1 = max(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha)));

    double m = (eps + m1) / 2;
    double M = M1 + eps;

    if (f(a, alpha) * d2fdx2(x0, alpha) > 0)
    {
        while ((fabs(x_2 - x_1) > eps * m / (M - m)) || i == 0)
        {
            x_1 = x_2;
            i++;
            x_2 = x_1 - f(x_1, alpha) / (f(x_1, alpha) - f(a, alpha)) * (x_1 - a);
        }
    }
    else
    {

        while ((fabs(x_2 - x_1) > eps * m / (M - m)) || i == 0)
        {
            x_1 = x_2;
            i++;
            x_2 = x_1 - f(x_1, alpha) / (f(b, alpha) - f(x_1, alpha)) * (b - x_1);
        }
    }

    fprintf(pFile, "%g, %g, %g, %d, %g, %g, - \n", x_2, f(x_2, alpha), x0, i, M, m);
}

void Etkin(FILE *pFile, double (*f)(double, double), double x0, double alpha, double a, double b) // Метод Ньютона (метод касательных)
{
    int i = 0; // Счетчик i
    double x = x0, x_2 = x0, x_1 = x0, x_12 = x0;

    double M1 = max(fabs(phi_dfdx(f, a, alpha, a, b)), fabs(phi_dfdx(f, b, alpha, a, b))) + 2 * eps;
    double q = (1 - eps * 2 + M1) / 2;

    x = phi(f, x0, alpha, a, b);
    x_2 = phi(f, x, alpha, a, b);

    while (fabs(x_2 - x_1) > eps * (1 - q) / q / 10 || i == 0)
    {
        x_1 = x;
        x = x_2;
        i++;
        x_12 = phi(f, x, alpha, a, b);
        x_2 = (x_1 * x_12 - x * x) / (x_1 - 2 * x + x_12);
    }

    fprintf(pFile, "%g, %g, %g, %d, -, -, %g \n", x_2, f(x_2, alpha), x0, i, q);
}

int main()
{
    FILE *pFile;
    pFile = fopen("report.csv", "w");

    for (int alpha = 1; alpha <= 10; alpha++)
    {
        double a = 0.95, b = 1, x0 = 1 - eps; // eps добавлен так как нельзя начинать метод хорд от какой-либо грнаницы
        fprintf(pFile, "α = %d, отрезок локализации: = [%g %g] \n", alpha, a, b);
        fprintf(pFile, "%s \n", ", x, f(x), x0, N, M, m, q");

        fprintf(pFile, "%s, ", "Половинное деление");
        Half_div(pFile, f, x0, alpha, a, b);

        fprintf(pFile, "%s, ", "Метод хорд");
        Chord(pFile, f, x0, alpha, a, b);

        fprintf(pFile, "%s, ", "Простая итерация");
        Simple_iteration(pFile, f, x0, alpha, a, b);

        fprintf(pFile, "%s, ", "Метод Эткена");
        Etkin(pFile, f, x0, alpha, a, b);

        fprintf(pFile, "%s, ", "Метод Ньютона");
        Newtone(pFile, f, x0, alpha, a, b);

        fprintf(pFile, "%s", "\n");
    }

    fclose(pFile);
    return 0;
}