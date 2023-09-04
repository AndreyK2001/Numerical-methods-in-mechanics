#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001   // Требуемая точность вычисления
#define inf 10000000 // Бесконечность (необходима для дельта-функции)

double sign(double x)
{
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
        return inf;
    else
        return 0;
}

double f(double x, double a)
{
    return fabs(x * x - a) - exp(a * fabs(x));
}

double diff(double x, double a)
{
    return 2 * x * sign(x * x - a) - a * exp(a * fabs(x)) * sign(x);
}

double diff2(double x, double a)
{
    return 2 * sign(x * x - a) + 2 * x * delta(x * x - a) - a * a * exp(a * fabs(x)) - a * exp(a * fabs(x)) * delta(x);
}

void Metod_chord(FILE *pFile, double (*f)(double, double), double x_0, double alpha, double a, double b)
{ // Метод хорд
    double x_1 = x_0, x_2 = x_0;

    int i = 0; // счётчик итераций алгоритма

    double m1 = min(fabs(diff(a, alpha)), fabs(diff(b, alpha)));
    double M1 = max(fabs(diff(a, alpha)), fabs(diff(b, alpha)));

    double m = (eps + m1) / 2;
    double M = M1 + eps;

    if (f(a, alpha) * diff2(x_0, alpha) > 0)
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

    fprintf(pFile, "%g, %g, %g, %d, %g, %g, - \n", x_2, f(x_2, alpha), x_0, i, M, m);
}

void Half_devision(FILE *pFile, double (*f)(double, double), double x_0, double alpha, double a, double b)
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

    fprintf(pFile, "%g, %g, %g, %d, -, -, - \n", x, f(x, alpha), x_0, i);
}

void Method_Newtona(FILE *pFile, double (*f)(double, double), double x_0, double alpha, double a, double b) // Метод Ньютона (метод касательных)
{
    int i = 0; // Число витков цикла n, счетчик i
    double x_2 = x_0, x_1 = x_0;

    double m1 = min(fabs(diff(a, alpha)), fabs(diff(b, alpha)));
    double M1 = max(fabs(diff(a, alpha)), fabs(diff(b, alpha)));

    double m = (eps + m1) / 2;
    double M = M1 + eps;

    while ((fabs(x_2 - x_1) > eps * m / (M - m)) || i == 0)
    {
        x_1 = x_2;
        i++;
        x_2 = x_1 - f(x_1, alpha) / diff(x_1, alpha);
    }

    fprintf(pFile, "%g, %g, %g, %d, %g, %g, - \n", x_2, f(x_2, alpha), x_0, i, M, m);
}

double phi(double (*f)(double, double), double x, double alpha, double a, double b)
{
    double M = max(fabs(diff(a, alpha)), fabs(diff(b, alpha))) + 3 * eps;
    return x - sign(diff(x, alpha)) * f(x, alpha) / M;
}

double phi_diff(double (*f)(double, double), double x, double alpha, double a, double b)
{
    double M = max(fabs(diff(a, alpha)), fabs(diff(b, alpha))) + 3 * eps;
    return 1 - delta(diff(x, alpha)) * diff2(x, alpha) * f(x, alpha) / M - sign(diff(x, alpha)) * diff(x, alpha) / M;
}

void Simple_iteration(FILE *pFile, double (*f)(double, double), double x_0, double alpha, double a, double b) // Метод простой итерации
{
    double x_2 = x_0, x_1 = x_0;
    int i = 0; // Число витков цикла n, счетчик i

    double M1 = max(fabs(phi_diff(f, a, alpha, a, b)), fabs(phi_diff(f, b, alpha, a, b))) + eps;
    double q = (1 - eps * 2 + M1) / 2;

    while (fabs(x_2 - x_1) > eps * (1 - q) / q / 10 || i == 0)
    {
        x_1 = x_2;
        i++;
        x_2 = phi(f, x_1, alpha, a, b);
    }

    fprintf(pFile, "%g, %g, %g, %d, -, -, %g \n", x_2, f(x_2, alpha), x_0, i, q);
}

void Method_Etkina(FILE *pFile, double (*f)(double, double), double x_0, double alpha, double a, double b) // Метод Ньютона (метод касательных)
{
    int i = 0; // Счетчик i
    double x = x_0, x_2 = x_0, x_1 = x_0, x_12 = x_0;

    double M1 = max(fabs(phi_diff(f, a, alpha, a, b)), fabs(phi_diff(f, b, alpha, a, b))) + 2 * eps;
    double q = (1 - eps * 2 + M1) / 2;

    x = phi(f, x_0, alpha, a, b);
    x_2 = phi(f, x, alpha, a, b);

    while (fabs(x_2 - x_1) > eps * (1 - q) / q / 10 || i == 0)
    {
        x_1 = x;
        x = x_2;
        i++;
        x_12 = phi(f, x, alpha, a, b);
        x_2 = (x_1 * x_12 - x * x) / (x_1 - 2 * x + x_12);
    }

    fprintf(pFile, "%g, %g, %g, %d, -, -, %g \n", x_2, f(x_2, alpha), x_0, i, q);
}

int main()
{
    FILE *pFile;
    pFile = fopen("answer.csv", "w");

    for (int alpha = 1; alpha <= 10; alpha++)
    {
        double a = -0.1, b = 0.5, x_0 = 0.5 - eps; // eps добавлен так как нельзя начинать метод хорд от какой-либо грнаницы
        fprintf(pFile, "alpha = %d, [a_%d b_%d] = [%g %g] \n", alpha, alpha, alpha, a, b);
        fprintf(pFile, "%s \n", ", x, f(x), x_0, N, M, m, q");

        fprintf(pFile, "%s, ", "Половинное деление");
        Half_devision(pFile, f, x_0, alpha, a, b);

        fprintf(pFile, "%s, ", "Метод хорд");
        Metod_chord(pFile, f, x_0, alpha, a, b);

        fprintf(pFile, "%s, ", "Простая итерация");
        Simple_iteration(pFile, f, x_0, alpha, a, b);

        fprintf(pFile, "%s, ", "Метод Эткена");
        Method_Etkina(pFile, f, x_0, alpha, a, b);

        fprintf(pFile, "%s, ", "Метод Ньютона");
        Method_Newtona(pFile, f, x_0, alpha, a, b);

        fprintf(pFile, "%s\n", "\n");
    }

    fclose(pFile);
    return 0;
}