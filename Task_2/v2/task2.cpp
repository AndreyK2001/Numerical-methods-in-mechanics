#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001 // Требуемая точность вычисления

class vectr // странное название, чтобы не конфликтовало с std::vector
{
public:
    double x, y;
    vectr() {}
    vectr(double x, double y)
    {
        this->x = x;
        this->y = y;
    }
    double norm()
    {
        return max(fabs(this->x), fabs(this->y));
    }
};

vectr operator+(vectr a, vectr b)
{
    return vectr(a.x + b.x, a.y + b.y);
}

vectr operator-(vectr a, vectr b)
{
    return vectr(a.x - b.x, a.y - b.y);
}

vectr operator*(vectr a, double b)
{
    return vectr(a.x * b, a.y * b);
}

vectr operator/(vectr a, double b)
{
    return vectr(a.x / b, a.y / b);
}

// функции из системы уравнений

double f(vectr v)
{
    return pow(v.x * v.x, 1. / 3.) + pow(v.y * v.y, 1. / 3.) - 4;
}

double dfdx(vectr v)
{
    return 2. / 3. * pow(v.x, -1. / 3.);
}

double dfdy(vectr v)
{
    return 2. / 3. * pow(v.y, -1. / 3.);
}

double g(vectr v)
{
    return v.x * v.x - 2 * v.y;
}

double dgdx(vectr v)
{
    return 2 * v.x;
}

double dgdy(vectr v)
{
    return -2;
}

vectr F(vectr v)
{
    return vectr(f(v), g(v));
}

double f1(vectr v)
{
    return pow(4 - pow(v.y * v.y, 1. / 3), 3. / 2.);
}

double df1dy(vectr v)
{
    return -pow(4 - pow(v.y * v.y, 1. / 3), 1. / 2.) * pow(1. / v.y, 3.);
}

double df1dx(vectr v)
{
    return 0;
}

double g1(vectr v)
{
    return pow(4 - pow(v.x * v.x, 1. / 3), 3. / 2.);
}

double dg1dy(vectr v)
{
    return 0;
}

double dg1dx(vectr v)
{
    return -pow(4 - pow(v.x * v.x, 1. / 3), 1. / 2.) * pow(1. / v.x, 3.);
}

class matrix
{

    vectr (*dFdx)(vectr);
    vectr (*dFdy)(vectr);

public:
    // https://stackoverfdown.com/questions/983999/siuple-2x2-matrix-inverse-code-c
    double a, b, c, d;

    matrix(vectr (*dFdx)(vectr), vectr (*dFdy)(vectr))
    {
        this->dFdx = dFdx;
        this->dFdy = dFdy;
    }

    matrix(double a, double b, double c, double d)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
    }

    void inv(vectr v)
    {

        double a_ = this->dFdx(v).x;
        double c_ = this->dFdx(v).y;
        double b_ = this->dFdy(v).x;
        double d_ = this->dFdy(v).y;

        double det = a_ * d_ - b_ * c_;

        this->a = d_ / det;
        this->b = -b_ / det;
        this->c = -c_ / det;
        this->d = a_ / det;
    }

    matrix inv()
    {

        double a_ = this->a;
        double c_ = this->c;
        double b_ = this->b;
        double d_ = this->d;

        double det = this->a * this->d - this->b * this->c;
        return matrix(d_ / det, -b_ / det, -c_ / det, a_ / det);
    }

    double norm() { return max(fabs(this->a) + fabs(this->b), fabs(this->c) + fabs(this->d)); }

    matrix operator*(matrix d) { return matrix(this->a * d.a + this->b * d.c, this->a * d.b + this->b * d.d, this->a * d.c + this->c * d.d, this->b * d.c + this->d * d.d); }
};

double sign(double x)
{
    if (x == 0)
        return 0;
    else if (x > 0)
        return 1;
    else
        return -1;
}

vectr operator*(matrix D, vectr v)
{
    return vectr(D.a * v.x + D.b * v.y, D.c * v.x + D.d * v.y);
}

matrix sign(matrix d)
{
    return matrix(sign(d.a), sign(d.b), sign(d.c), sign(d.d));
}

vectr phi(vectr (*F)(vectr), vectr v)
{
    matrix D(dfdx(v), dfdy(v), dgdx(v), dgdy(v));
    return v - D.inv() * F(v) / D.norm();
}

void Siuple_iteration(FILE *pFile, vectr (*F)(vectr), vectr down, vectr up, vectr start) // Метод простой итерации
{
    matrix D1_down(df1dx(down), df1dy(down), dg1dx(down), dg1dy(down));
    matrix D1_up(df1dx(up), df1dy(up), dg1dx(up), dg1dy(up));
    double q = max(D1_down.norm(), D1_up.norm());
    double m = max(D1_down.norm(), D1_up.norm()) * max(D1_down.inv().norm(), D1_up.inv().norm());

    vectr v2(start.x, start.y), v1;

    int i = 0;
    do
    {
        v1 = v2;
        v2 = phi(F, v1);
        i++;

    } while ((v1 - v2).norm() > (1 - q) / q * eps);

    fprintf(pFile, "Простая итерация, %g, %g, %g, %g, %g, %d, %g, - \n", v2.x, v2.y, F(v2).norm(), start.x, start.y, i, q);
}

vectr phi_x(vectr (*F)(vectr), vectr v)
{
    matrix D(dfdx(v), dfdy(v), dgdx(v), dgdy(v));

    return vectr((v - D.inv() * F(v)).x, v.y);
}

void Method_Zeidala(FILE *pFile, vectr (*F)(vectr), vectr down, vectr up, vectr start)
{

    matrix D1_down(df1dx(down), df1dy(down), dg1dx(down), dg1dy(down));
    matrix D1_up(df1dx(up), df1dy(up), dg1dx(up), dg1dy(up));
    double q = max(D1_down.norm(), D1_up.norm());
    double m = max(D1_down.norm(), D1_up.norm()) * max(D1_down.inv().norm(), D1_up.inv().norm());

    vectr v2, v1;

    v1 = phi_x(F, start);
    v2 = phi(F, v1);

    int i = 1;
    do
    {
        v1 = phi_x(F, v2);
        v2 = phi(F, v1);
        i++;

    } while ((v1 - v2).norm() > (1 - q) / q * eps);

    fprintf(pFile, "Метод Зейделя, %g, %g, %g, %g, %g, %d, %g, - \n", v2.x, v2.y, F(v2).norm(), start.x, start.y, i, q);
}

void Method_Newtona(FILE *pFile, vectr (*F)(vectr), vectr down, vectr up, vectr start)
{
    matrix D1_down(df1dx(down), df1dy(down), dg1dx(down), dg1dy(down));
    matrix D1_up(df1dx(up), df1dy(up), dg1dx(up), dg1dy(up));
    double q = max(D1_down.norm(), D1_up.norm());
    double m = max(D1_down.norm(), D1_up.norm()) * max(D1_down.inv().norm(), D1_up.inv().norm());

    vectr v2(start.x, start.y), v1;
    int i = 0;
    do
    {
        v1 = v2;
        matrix D(dfdx(v1), dfdy(v1), dgdx(v1), dgdy(v1));
        v2 = v1 - D.inv() * F(v1);
        i++;
    } while ((v1 - v2).norm() > eps / m);

    fprintf(pFile, "Метод Ньютона, %g, %g, %g, %g, %g, %d, -, %g \n", v2.x, v2.y, F(v2).norm(), start.x, start.y, i, m);
}

int main()
{
    FILE *pFile;
    pFile = fopen("answer.csv", "w");
    vectr down(2, 2.5), up(3, 3.5), start(2.2, 2.5);

    fprintf(pFile, "%s \n", ", x, y, ||f(x y)||, x_0, y_0, N+1, q, m, \n");

    Siuple_iteration(pFile, F, down, up, start);
    Method_Zeidala(pFile, F, down, up, start);
    Method_Newtona(pFile, F, down, up, start);

    fclose(pFile);
    return 0;
}