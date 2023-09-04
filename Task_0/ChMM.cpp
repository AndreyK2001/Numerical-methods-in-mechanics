#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001    // Требуемая точность вычисления
#define e 0.00005     // Минимальное значение невязки e
#define left_bound 0. // Левая и правая граница отрезка локализации корня
#define right_bound 5.
#define ksi0 2.   // Начально приближение корня
#define x_real 1. //Корень, рассчитанный аналитически

double f(double x) //функция (x-1)^5/3+x-1
{
    if (x > 1)
        return pow(x - 1, 5. / 3) + x - 1.;
    else
        return -pow(1 - x, 5. / 3) + x - 1.;
}

double ddx(double x, double (*f)(double)) //Производная
{
    return (f(x + eps / 1000) - f(x)) / (eps / 1000);
}

int horda(char verbose_mode = 0, char *file_name = "algorithms.txt", bool rewrite = false) //Метод хорд (метод секущих)
/*
0 - нет вывода ни в консоль ни в файл
1 - вывод только в файл
2 - вывод только в консоль
3 - вывод и в консоль и в файл
*/
{
    bool verbose_file, verbose_console;
    switch (verbose_mode)
    {
    case 1:
        verbose_file = true;
        verbose_console = false;
        break;
    case 2:
        verbose_file = false;
        verbose_console = true;
        break;
    case 3:
        verbose_file = true;
        verbose_console = true;
        break;

    default:
        verbose_file = false;
        verbose_console = false;
        break;
    }
    double a = left_bound, b = right_bound, ksi = ksi0; // Концы отрезка и искомый корень
    int n, i = 0;                                       // Число витков цикла n, счетчик i

    ofstream fout; // Открываем поток

    if (verbose_file)
    {
        if (rewrite)fout.open(file_name);
        else fout.open(file_name, ios_base::out | ios_base::app); // Связываем поток с файлом
        fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом хорд" << endl;
        fout << "Начальные приближения: "
             << "отрезок: [" << a << ", " << b << "], "
             << "Начальное приближение корня: " << ksi0 << ". Корень вычисляется с точностью: " << eps << endl;
    }

    double derivative2nd = (ddx(ksi + eps / 1000, f) - ddx(ksi, f)) / (eps / 1000); //вторая производная
    if (f(a) * derivative2nd > 0)
    {
        if (verbose_file)
            fout << "Левая граница неподвижна" << endl; /* Вывод результатов в файл*/
        if (verbose_console)
            cout << "Левая граница неподвижна" << endl;
        while (fabs(ksi - x_real) > eps)
        {
            i++;
            if (verbose_file)
                fout << "i=" << i << " ksi=" << ksi << endl; /* Вывод результатов в файл*/
            if (verbose_console)
                cout << "i=" << i << " ksi=" << ksi << endl;
            ksi = ksi - f(ksi) / (f(ksi) - f(a)) * (ksi - a);
        }
    }
    else
    {
        if (verbose_file)
            fout << "Правая граница неподвижна" << endl; /* Вывод результатов в файл*/
        if (verbose_console)
            cout << "Правая граница неподвижна" << endl;
        while (fabs(ksi - x_real) > eps)
        {
            i++;
            if (verbose_file)
                fout << "i=" << i << " ksi=" << ksi << endl; /* Вывод результатов в файл*/
            if (verbose_console)
                cout << "i=" << i << " ksi=" << ksi << endl;
            ksi = ksi - f(ksi) / (f(b) - f(ksi)) * (b - ksi);
        }
    }

    if (verbose_console)
        cout << "ksi = " << ksi << endl;
    if (verbose_file)
        fout << "ksi = " << ksi << endl;
    return i;
}

int newtone(char verbose_mode = 0, char *file_name = "algorithms.txt", bool rewrite = false) //Метод Ньютона (метод касательных)
/*
0 - нет вывода ни в консоль ни в файл
1 - вывод только в файл
2 - вывод только в консоль
3 - вывод и в консоль и в файл
*/
{
    bool verbose_file, verbose_console;
    switch (verbose_mode)
    {
    case 1:
        verbose_file = true;
        verbose_console = false;
        break;
    case 2:
        verbose_file = false;
        verbose_console = true;
        break;
    case 3:
        verbose_file = true;
        verbose_console = true;
        break;

    default:
        verbose_file = false;
        verbose_console = false;
        break;
    }
    double a = left_bound, b = right_bound, ksi = ksi0; // Концы отрезка и искомый корень
    int n, i = 0;                                       // Число витков цикла n, счетчик i

    ofstream fout; // Открываем поток
    if (verbose_file)
    {
        if (rewrite)fout.open(file_name);
        else fout.open(file_name, ios_base::out | ios_base::app); // Связываем поток с файлом
        fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом Ньютона" << endl;
        fout << "Начальные приближения: "
             << "отрезок: [" << a << ", " << b << "], "
             << "Начальное приближение корня: " << ksi0 << ". Корень вычисляется с точностью: " << eps << endl;
    }
    while (fabs(ksi - x_real) > eps)
    {
        i++;
        if (verbose_file)
            fout << "i=" << i << " ksi=" << ksi << endl; /* Вывод результатов в файл*/
        if (verbose_console)
            cout << "i=" << i << " ksi=" << ksi << endl;
        ksi = ksi - f(ksi) / ddx(ksi, f);
    }

    if (verbose_console)
        cout << "ksi = " << ksi << endl;
    if (verbose_file)
        fout << "ksi = " << ksi << endl;
    return i;
}

double sign(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

double phi(double x, double (*f)(double))
{
    double M = fmax(fabs(ddx(left_bound, f)), fabs(ddx(right_bound, f)));
    return x - sign(ddx(x, f)) * f(x) / M;
}

int simple_iteration(char verbose_mode = 0, char *file_name = "algorithms.txt", bool rewrite = false) //Метод простой итерации
/*
0 - нет вывода ни в консоль ни в файл
1 - вывод только в файл
2 - вывод только в консоль
3 - вывод и в консоль и в файл
*/
{
    bool verbose_file, verbose_console;
    switch (verbose_mode)
    {
    case 1:
        verbose_file = true;
        verbose_console = false;
        break;
    case 2:
        verbose_file = false;
        verbose_console = true;
        break;
    case 3:
        verbose_file = true;
        verbose_console = true;
        break;

    default:
        verbose_file = false;
        verbose_console = false;
        break;
    }

    double a = left_bound, b = right_bound, ksi = ksi0; // Концы отрезка и искомый корень
    int n, i = 0;                                       // Число витков цикла n, счетчик i
    ofstream fout;                                      // Открываем поток
    if (verbose_file)
    {
        if (rewrite)fout.open(file_name);
        else fout.open(file_name, ios_base::out | ios_base::app); // Связываем поток с файлом
        fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом простой итерации" << endl;
        fout << "Начальные приближения: "
             << "отрезок: [" << a << ", " << b << "], "
             << "Начальное приближение корня: " << ksi0 << ". Корень вычисляется с точностью: " << eps << endl;
    }

    while (fabs(ksi - x_real) > eps)
    {
        i++;
        if (verbose_file)
            fout << "i=" << i << " ksi=" << ksi << endl; /* Вывод результатов в файл*/
        if (verbose_console)
            cout << "i=" << i << " ksi=" << ksi << endl;
        ksi = phi(ksi, f);
    }

    if (verbose_console)
        cout << "ksi = " << ksi << endl;
    if (verbose_file)
        fout << "ksi = " << ksi << endl;
    fout.close(); // Закрываем поток
    return i;
}

int half(char verbose_mode = 0, char *file_name = "algorithms.txt", bool rewrite = false) //Метод деления отрезка пополам (метод половинного деления)
/*
0 - нет вывода ни в консоль ни в файл
1 - вывод только в файл
2 - вывод только в консоль
3 - вывод и в консоль и в файл
*/
{
    bool verbose_file, verbose_console;
    switch (verbose_mode)
    {
    case 1:
        verbose_file = true;
        verbose_console = false;
        break;
    case 2:
        verbose_file = false;
        verbose_console = true;
        break;
    case 3:
        verbose_file = true;
        verbose_console = true;
        break;

    default:
        verbose_file = false;
        verbose_console = false;
        break;
    }

    double a = left_bound, b = right_bound, ksi; // Концы отрезка и искомый корень
    int n, i;                                    // Число витков цикла n, счетчик i
    ofstream fout;                               // Открываем поток
    if (verbose_file)
    {
        if (rewrite)fout.open(file_name);
        else fout.open(file_name, ios_base::out | ios_base::app); // Связываем поток с файлом
        fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом деления отрезка пополам" << endl;
        fout << "Начальные приближения: "
             << "отрезок: [" << a << ", " << b << "], "
             << "Начальное приближение корня: " << ksi0 << ". Корень вычисляется с точностью: " << eps << endl;
    }
    n = int(log((b - a) / eps) / log(2)) + 1; // n D log2 " C 1
    ksi = (a + b) / 2.;                       // Вычисляем середину начального отрезка
    for (i = 1; i <= n && fabs(f(ksi)) >= e; i++)
    {
        if (f(ksi) * f(a) < 0)
            b = ksi;
        else
            a = ksi;

        if (verbose_console)
            cout << "i=" << i << " [a,b]=[" << a << "," << b << "]" << endl; /* Вывод результатов на консоль*/
        if (verbose_file)
            fout << "i=" << i << " [a,b]=[" << a << "," << b << "]" << endl; /* Вывод результатов в файл*/
        ksi = (a + b) / 2.;                                                  // Вычисляем середину нового отрезка
    }

    if (verbose_console)
        cout << "ksi = " << ksi << endl;
    if (verbose_file)
        fout << "ksi = " << ksi << endl;
    fout.close(); // Закрываем поток
    return i;
}

int main()
{
    horda(1);
    newtone(1);
    simple_iteration(1);
    half(1);
    return 0;
}