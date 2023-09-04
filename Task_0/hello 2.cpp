#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001 // Требуемая точность вычисления
#define e 0.00005  // Минимальное значение невязки e
#define a 0.       // Левая и правая граница отрезка локализации корня
#define b 5.
#define x_0 2.    // Начально приближение корня
#define x_real 1. //Корень, рассчитанный аналитически

double f(double x) //функция (x-1)^5/3+x-1
{
    if (x > 1)
        return pow(x - 1, 5. / 3) + x - 1.;
    else
        return -pow(1 - x, 5. / 3) + x - 1.;
}

double dfdx(double (*f)(double), double x) //Производная
{
    return (f(x + eps / 1000) - f(x)) / (eps / 1000);
}

int Chord_method(double (*f)(double), double left, double right, double x_start) //Метод хорд (метод секущих)
{
    double x = x_start;

    ofstream fout; // Открываем поток
    fout.open("Chord.txt");

    fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом хорд на отрезке "
         << "[" << left << ", " << right << "], " << endl;
    fout << "Начальное приближение корня: " << x << ". Корень вычисляется с точностью: " << eps << endl;

    double d2fdx2 = (dfdx(f, x_start + eps / 1000) - dfdx(f, x_start)) / (eps / 1000); //вторая производная

    int i = 0; //счётчик итераций алгоритма

    if (f(left) * d2fdx2 > 0)
    {
        fout << "Левая граница неподвижна" << endl; /* Вывод результатов в файл*/

        while (fabs(x - x_real) > eps)
        {
            i++;
            fout << "i = " << i << " x = " << x << endl; /* Вывод результатов в файл*/
            x = x - f(x) / (f(x) - f(left)) * (x - left);
        }
    }
    else
    {

        fout << "Правая граница неподвижна" << endl; /* Вывод результатов в файл*/

        while (fabs(x - x_real) > eps)
        {
            i++;
            fout << "i = " << i << " x = " << x << endl; /* Вывод результатов в файл*/
            x = x - f(x) / (f(right) - f(x)) * (right - x);
        }
    }

    fout << "После " << i << " итераций получен корень x = " << x << endl;
    fout.close(); // Закрываем поток
    return i;
}

int Newtones_method(double (*f)(double), double left, double right, double x_start) //Метод Ньютона (метод касательных)
{
    int i = 0; // Число витков цикла n, счетчик i
    double x = x_start;
    ofstream fout; // Открываем поток

    fout.open("Newtones_method.txt"); // Связываем поток с файлом

    fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом Ньютона на отрезке "
         << "[" << left << ", " << right << "], " << endl;
    fout << "Начальное приближение x = " << x << ". Корень вычисляется с точностью: " << eps << endl;

    while (fabs(x - x_real) > eps)
    {
        i++;
        fout << "i = " << i << " x = " << x << endl; /* Вывод результатов в файл*/
        x = x - f(x) / dfdx(f, x);
    }

    fout << "После " << i << " итераций получен корень x = " << x << endl;
    fout.close(); // Закрываем поток
    return i;
}

double sign(double x)
{
    if (x == 0)
        return 0;
    else
        return x / fabs(x);
}

double phi(double (*f)(double), double x)
{
    double M = fmax(fabs(dfdx(f, a)), fabs(dfdx(f, b)));
    return x - sign(dfdx(f, x)) * f(x) / M;
}

int SimpleIteration_method(double (*f)(double), double left, double right, double x_start) //Метод простой итерации
{
    double x = x_start;
    int i = 0;     // Число витков цикла n, счетчик i
    ofstream fout; // Открываем поток

    fout.open("SimpleIteration.txt");

    fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом простой итерации на отрезке "
         << "[" << a << ", " << b << "], " << endl;
    fout << "Начальное приближение корня: " << x << ". Корень вычисляется с точностью: " << eps << endl;

    while (fabs(x - x_real) > eps)
    {
        i++;
        fout << "i = " << i << " x = " << x << endl; /* Вывод результатов в файл*/
        x = phi(f, x);
    }

    fout << "После " << i << " итераций получен корень x = " << x << endl;
    fout.close(); // Закрываем поток
    return i;
}

int Half_method(double (*f)(double), double left, double right, double x_start) //Метод деления отрезка пополам (метод половинного деления)
{
    int n, i;                                             // Число витков цикла n, счетчик i
    double x = (left + right) / 2, a_ = left, b_ = right; // Вычисляем середину начального отрезка
    ofstream fout;                                        // Открываем поток

    fout.open("Half_method.txt");

    fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом деления отрезка пополам "
         << "[" << a_ << ", " << b_ << "], " << endl;
    fout << "Начальное приближение корня: " << x << ". Корень вычисляется с точностью: " << eps << endl;

    n = int(log((b_ - a_) / eps) / log(2)) + 1;

    for (i = 1; i <= n && fabs(x - x_real) > eps; i++)
    {
        if (f(x) * f(a_) < 0)
            b_ = x;
        else
            a_ = x;

        fout << "i = " << i << " [a,b]=[" << a_ << "," << b_ << "]" << endl; /* Вывод результатов в файл*/
        x = (a_ + b_) / 2.;                                                  // Вычисляем середину нового отрезка
    }

    fout << "После " << i << " итераций получен корень x = " << x << endl;
    fout.close(); // Закрываем поток
    return i;
}

int main()
{
    cout << "Количество итераций метода хорд " << Chord_method(f, a, b, x_0) << endl;
    cout << "Количество итераций метода Ньютона " << Newtones_method(f, a, b, x_0) << endl;
    cout << "Количество итераций метода простой итерации " << SimpleIteration_method(f, a, b, x_0) << endl;
    cout << "Количество итераций метода деления отрезка пополам " << Half_method(f, a, b, x_0) << endl;

    return 0;
}