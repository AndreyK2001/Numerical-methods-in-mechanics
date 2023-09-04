#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001 // Требуемая точность вычисления
#define inf 10000000

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
{
    if (x == 0)
        return inf;
    else
        return 0;
}

double f(double x, double a)
{
    return 10 * x - 0.5 * pow(x, 3) + pow(x, 5) - pow(a, 3) / (25 * exp(pow(x, 2)));
}

double dfdx(double x, double a)
{
    return 10 + 5 * pow(x, 4) - 3 * pow(x, 2) / 2 + 2 * x * pow(a, 3) * exp(-pow(x, 2)) / 25;
}

double d2fdx2(double x, double a)
{
    return -3 * x + 20 * pow(x, 3) + 2 * pow(a, 3) * exp(-pow(x, 2)) / 25 - 4 * pow(a, 3) * pow(x, 2) * exp(-pow(x, 2)) / 25;
}

int HalfMethod(double (*f)(double, double), double x0, double alpha, double a, double b) // Метод деления отрезка пополам (метод половинного деления)
{
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
    cout << x << ", " << f(x, alpha) << ", " << x0 << ", " << i << ", "
         << "-"
         << ", "
         << "-"
         << ", "
         << "-" << endl;
    return i;
}

int ChordMethod(double (*f)(double, double), double x0, double alpha, double a, double b) // Метод хорд
{
    double x_priv, x_next = x0;

    int i = 0; // счётчик итераций алгоритма

    double m1 = min(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha)));
    double M1 = max(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha)));

    double m = (eps + m1) / 2;
    double M = M1 + eps;

    if (f(a, alpha) * d2fdx2(x0, alpha) > 0)
    {
        while ((fabs(x_next - x_priv) > eps * m / (M - m)) || i == 0)
        {
            x_priv = x_next;
            i++;
            x_next = x_priv - f(x_priv, alpha) / (f(x_priv, alpha) - f(a, alpha)) * (x_priv - a);
        }
    }
    else
    {

        while ((fabs(x_next - x_priv) > eps * m / (M - m)) || i == 0)
        {
            x_priv = x_next;
            i++;
            x_next = x_priv - f(x_priv, alpha) / (f(b, alpha) - f(x_priv, alpha)) * (b - x_priv);
        }
    }

    cout << x_next << ", " << f(x_next, alpha) << ", " << x0 << ", " << i << ", " << M << ", " << m << ", "
         << "-" << endl;
    return i;
}

double phi(double (*f)(double, double), double x, double alpha, double a, double b)
{
    double M = max(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha))) + 4 * eps;
    return x - sign(dfdx(x, alpha)) * f(x, alpha) / M;
}

double dphidx(double (*f)(double, double), double x, double alpha, double a, double b)
{
    double M = max(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha))) + 4 * eps;
    return 1 - delta(dfdx(x, alpha)) * d2fdx2(x, alpha) * f(x, alpha) / M - sign(dfdx(x, alpha)) * dfdx(x, alpha) / M;
}

int SimpleIteration(double (*f)(double, double), double x0, double alpha, double a, double b) // Метод простой итерации
{
    double x_next = x0, x_priv = x0;
    int i = 0; // Число витков цикла n, счетчик i

    double M1 = max(fabs(dphidx(f, a, alpha, a, b)), fabs(dphidx(f, b, alpha, a, b))) + eps;
    double q = (1 - eps * 2 + M1) / 2;

    while ((fabs(x_next - x_priv) > eps * (1 - q) / q) || i == 0)
    {
        x_priv = x_next;
        i++;
        x_next = phi(f, x_priv, alpha, a, b);
    }

    cout << x_next << ", " << f(x_next, alpha) << ", " << x0 << ", " << i << ", "
         << "-"
         << ", "
         << "-"
         << ", "
         << q << endl;
    return i;
}

int Etkin(double (*f)(double, double), double x0, double alpha, double a, double b) // Метод Ньютона (метод касательных)
{
    int i = 0; // Счетчик i
    double x = x0, x_next = x0, x_priv = x0, x_mid = x0;

    double M1 = max(fabs(dphidx(f, a, alpha, a, b)), fabs(dphidx(f, b, alpha, a, b))) + eps;
    double q = (1 - eps * 2 + M1) / 2;

    x = phi(f, x0, alpha, a, b);
    x_next = phi(f, x, alpha, a, b);

    while ((fabs(x_next - x_priv) > eps * (1 - q) / q) || i == 0)
    {
        x_priv = x;
        x = x_next;
        i++;
        x_mid = phi(f, x, alpha, a, b);
        x_next = (x_priv * x_mid - pow(x, 2)) / (x_priv - 2 * x + x_mid);
    }

    cout << x_next << ", " << f(x_next, alpha) << ", " << x0 << ", " << i << ", "
         << "-"
         << ", "
         << "-"
         << ", "
         << q << endl;
    return i;
}

int Newtone(double (*f)(double, double), double x0, double alpha, double a, double b) // Метод Ньютона (метод касательных)
{
    int i = 0; // Число витков цикла n, счетчик i
    double x_next = x0, x_priv = x0;

    double m1 = min(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha)));
    double M1 = max(fabs(dfdx(a, alpha)), fabs(dfdx(b, alpha)));

    double m = (eps + m1) / 2;
    double M = M1 + eps;

    while ((fabs(x_next - x_priv) > eps * m / (M - m)) || i == 0)
    {
        x_priv = x_next;
        i++;
        x_next = x_priv - f(x_priv, alpha) / dfdx(x_priv, alpha);
    }

    cout  << x_next << ", " << f(x_next, alpha) << ", " << x0 << ", " << i << ", " << M << ", " << m << ", "
         << "-" << endl;

    return i;
}

int main()
{
    for (int alpha = 1; alpha <= 10; alpha++)
    {

        double x0 = 0.1356 * alpha;
        double a = 0;
        double b = x0 - 2 * eps;

        cout << "alpha = " << alpha << ", [a; b] = [" << a << "; " << b << "]" << endl;
        cout << ", x, f(x), x0, N, M, m, q" << endl;

        cout << "Половинное деление, ";
        HalfMethod(f, x0, alpha, a, b);

        cout << "Метод хорд, ";
        ChordMethod(f, x0, alpha, a, b);

        cout << "Проятая итерация, ";
        SimpleIteration(f, x0, alpha, a, b);

        cout << "Метод Эткена, ";
        Etkin(f, x0, alpha, a, b);

        cout << "Метод Ньютона, ";
        Newtone(f, x0, alpha, a, b);

        cout << "\n"
             << endl;
    }
    return 0;
}