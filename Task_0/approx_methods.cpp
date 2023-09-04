#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001 // Требуемая точность вычисления
#define e 0.00005  // Минимальное значение невязки e
#define left 0.    // Левый и правый конец отрезка, на котором ведётся поиск
#define right 5.
#define ksi 2.    // Начально приближение корня
#define x_real 1. //Корень, рассчитанный аналитически

double f(double x) //функция (x-1)^5/3+x-1
{
    if (x > 1)
        return pow(x - 1, 5. / 3) + x - 1.;
    else
        return -pow(1 - x, 5. / 3) + x - 1.;
}

double dfdx(double x) // 1я производная
{
    if (x > 1)
        return 5. / 3 * pow(x - 1, 2. / 3) + 1.;
    else
        return -5. / 3 * pow(1 - x, 2. / 3) + 1.;
}

double d2fdx2(double x) // 2я производная
{
    if (x > 1)
        return 5. / 3 * 2. / 3 * pow(x - 1, -1. / 3);
    else
        return -5. / 3 * 2. / 3 * pow(1 - x, -1. / 3);
}

int ChordMethod(char *file_name = "Chord_output.txt") //Метод хорд (метод секущих)
{
    double a = left, b = right, x = ksi; // Концы отрезка и искомый корень
    int n, i = 0;                        // Число витков цикла n, счетчик i

    ofstream fout; // Открываем поток

    fout.open(file_name);

    fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом хорд на отрезке "
         << "[" << a << ", " << b << "], " << endl;
    fout << "Начальное приближение корня: " << x << ". Корень вычисляется с точностью: " << eps << endl;

    if (f(a) * d2fdx2(x) > 0)
    {
        fout << "Левая граница неподвижна" << endl; /* Вывод результатов в файл*/

        while (fabs(x - x_real) > eps)
        {
            i++;
            fout << "i = " << i << " x = " << x << endl; /* Вывод результатов в файл*/
            x = x - f(x) / (f(x) - f(a)) * (x - a);
        }
    }
    else
    {

        fout << "Правая граница неподвижна" << endl; /* Вывод результатов в файл*/

        while (fabs(x - x_real) > eps)
        {
            i++;
            fout << "i = " << i << " x = " << x << endl; /* Вывод результатов в файл*/
            x = x - f(x) / (f(b) - f(x)) * (b - x);
        }
    }

    fout << "После " << i << " итераций получен корень x = " << x << endl;
    fout.close(); // Закрываем поток
    return i;
}

int NewtoneMethod(char *file_name = "Newtone_output.txt") //Метод Ньютона (метод касательных)
{
    double a = left, b = right, x = ksi; // Концы отрезка и искомый корень
    int n, i = 0;                        // Число витков цикла n, счетчик i

    ofstream fout; // Открываем поток

    fout.open(file_name); // Связываем поток с файлом

    fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом Ньютона на отрезке "
         << "[" << a << ", " << b << "], " << endl;
    fout << "Начальное приближение корня: " << x << ". Корень вычисляется с точностью: " << eps << endl;

    while (fabs(x - x_real) > eps)
    {
        i++;
        fout << "i = " << i << " x = " << x << endl; /* Вывод результатов в файл*/
        x = x - f(x) / dfdx(x);
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
    {
        if (x > 0)
            return 1;
        else
            return -1;
    }
}

double phi(double x, double (*f)(double))
{
    double M = fmax(fabs(dfdx(left)), fabs(dfdx(right)));
    return x - sign(dfdx(x)) * f(x) / M;
}

int SimpleIteration(char *file_name = "SimpleIteration_output.txt") //Метод простой итерации
{
    double a = left, b = right, x = ksi; // Концы отрезка и искомый корень
    int n, i = 0;                        // Число витков цикла n, счетчик i
    ofstream fout;                       // Открываем поток

    fout.open(file_name);

    fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом простой итерации на отрезке "
         << "[" << a << ", " << b << "], " << endl;
    fout << "Начальное приближение корня: " << x << ". Корень вычисляется с точностью: " << eps << endl;

    while (fabs(x - x_real) > eps)
    {
        i++;
        fout << "i = " << i << " x = " << x << endl; /* Вывод результатов в файл*/
        x = phi(x, f);
    }

    fout << "После " << i << " итераций получен корень x = " << x << endl;
    fout.close(); // Закрываем поток
    return i;
}

int HalfMethod(char *file_name = "Half_output.txt") //Метод деления отрезка пополам (метод половинного деления)
{
    double a = left, b = right; // Концы отрезка и искомый корень
    int n, i;               // Число витков цикла n, счетчик i
    double x = (a + b) / 2;     // Вычисляем середину начального отрезка
    ofstream fout;              // Открываем поток

    fout.open(file_name);

    fout << "Поиск корня уравнения (x-1)^5/3+x-1=0 методом деления отрезка пополам "
         << "[" << a << ", " << b << "], " << endl;
    fout << "Начальное приближение корня: " << x << ". Корень вычисляется с точностью: " << eps << endl;

    n = int(log((b - a) / eps) / log(2)) + 1;

    for (i = 1; i <= n && fabs(x - x_real) > eps; i++)
    {
        if (f(x) * f(a) < 0)
            b = x;
        else
            a = x;

        fout << "i = " << i << " [a,b]=[" << a << "," << b << "]" << endl; /* Вывод результатов в файл*/
        x = (a + b) / 2.;                                                  // Вычисляем середину нового отрезка
    }

    fout << "После " << i << " итераций получен корень x = " << x << endl;
    fout.close(); // Закрываем поток
    return i;
}

int main()
{
    cout << "Количество итераций метода хорд " << ChordMethod() << endl;
    cout << "Количество итераций метода Ньютона " << NewtoneMethod() << endl;
    cout << "Количество итераций метода простой итерации " << SimpleIteration() << endl;
    cout << "Количество итераций метода деления отрезка пополам " << HalfMethod() << endl;
    
    return 0;
}