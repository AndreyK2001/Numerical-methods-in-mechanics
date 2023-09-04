#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001 // Требуемая точность вычисления

class vector2
{
public:
    double x, y;

    vector2(double x, double y)
    {
        this->x = x;
        this->y = y;
    }

    vector2() {}

    double norm() { return max(fabs(this->x), fabs(this->y)); }
};

vector2 operator+(vector2 a, vector2 b) { return vector2(a.x + b.x, a.y + b.y); }

vector2 operator-(vector2 a, vector2 b) { return vector2(a.x - b.x, a.y - b.y); }

vector2 operator*(vector2 a, double b) { return vector2(a.x * b, a.y * b); }

vector2 operator/(vector2 a, double b) { return vector2(a.x / b, a.y / b); }

class Dmatrix
{

    vector2 (*dFdx)(vector2);
    vector2 (*dFdy)(vector2);

public:
    double a, b, c, d;

    Dmatrix(vector2 (*dFdx)(vector2), vector2 (*dFdy)(vector2))
    {
        this->dFdx = dFdx;
        this->dFdy = dFdy;
    }

    Dmatrix(double a, double b, double c, double d)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
    }

    void point(vector2 v)
    {
        this->a = this->dFdx(v).x;
        this->c = this->dFdx(v).y;
        this->b = this->dFdy(v).x;
        this->d = this->dFdy(v).y;
    }

    void inverse(vector2 v)
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

    Dmatrix inverse()
    {

        double a_ = this->a;
        double c_ = this->c;
        double b_ = this->b;
        double d_ = this->d;

        double det = this->a * this->d - this->b * this->c;
        return Dmatrix(d_ / det, -b_ / det, -c_ / det, a_ / det);
    }

    double norm() { return max(fabs(this->a) + fabs(this->b), fabs(this->c) + fabs(this->d)); }

    Dmatrix operator*(Dmatrix d) { return Dmatrix(this->a * d.a + this->b * d.c, this->a * d.b + this->b * d.d, this->a * d.c + this->c * d.d, this->b * d.c + this->d * d.d); }
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

vector2 operator*(Dmatrix D, vector2 v) { return vector2(D.a * v.x + D.b * v.y, D.c * v.x + D.d * v.y); }

double f(vector2 v) { return tan(v.x) - cos(1.5 * v.y); }

double dfdx(vector2 v) { return pow(cos(v.x), -2); }

double dfdy(vector2 v) { return 1.5 * sin(1.5 * v.y); }

double g(vector2 v) { return 2 * v.y * v.y - v.x * v.x + 4 * v.x - 3; }

double dgdx(vector2 v) { return -2 * v.x + 4; }

double dgdy(vector2 v) { return 4 * v.y; }

vector2 F(vector2 v) { return vector2(f(v), g(v)); }

double f1(vector2 v) { return atan(cos(1.5 * v.y)); }

double df1dy(vector2 v) { return -1.5 * sin(1.5 * v.y) / (1 + cos(1.5 * v.y) * cos(1.5 * v.y)); }

double df1dx(vector2 v) { return 0; }

double g1(vector2 v) { return sqrt((pow(v.x - 2, 2) - 1) / 2); }

double dg1dy(vector2 v) { return 0; }

double dg1dx(vector2 v) { return (v.x - 2) / (sqrt(2) * sqrt((v.x - 2) * (v.x - 2) - 1)); }

Dmatrix sign(Dmatrix d) { return Dmatrix(sign(d.a), sign(d.b), sign(d.c), sign(d.d)); }

vector2 phi(vector2 (*F)(vector2), vector2 v)
{
    Dmatrix D(dfdx(v), dfdy(v), dgdx(v), dgdy(v));
    return v - F(v) / D.norm();
}

void SimpleIteration(vector2 (*F)(vector2), vector2 low, vector2 up, vector2 start) // Метод простой итерации
{
    Dmatrix D1_low(df1dx(low), df1dy(low), dg1dx(low), dg1dy(low));
    Dmatrix D1_up(df1dx(up), df1dy(up), dg1dx(up), dg1dy(up));
    double q = max(D1_low.norm(), D1_up.norm());
    double u = max(D1_low.norm(), D1_up.norm()) * max(D1_low.inverse().norm(), D1_up.inverse().norm());

    vector2 v2(start.x, start.y), v1;

    int i = 0;
    do
    {
        v1 = v2;
        v2 = phi(F, v1);
        i++;

    } while ((v1 - v2).norm() > (1 - q) / q * eps);

    cout << "Простая итерация, " << v2.x << ", " << v2.y << ", " << F(v2).norm() << ", " << start.x << ", " << start.y << ", "
         << i << ", " << q << ", " << u << endl;
}

vector2 phi_(vector2 (*F)(vector2), vector2 v)
{
    Dmatrix D(dfdx(v), dfdy(v), dgdx(v), dgdy(v));

    return vector2((v - F(v) / D.norm()).x, v.y);
}

void Zeidels_method(vector2 (*F)(vector2), vector2 low, vector2 up, vector2 start)
{

    Dmatrix D1_low(df1dx(low), df1dy(low), dg1dx(low), dg1dy(low));
    Dmatrix D1_up(df1dx(up), df1dy(up), dg1dx(up), dg1dy(up));
    double q = max(D1_low.norm(), D1_up.norm());
    double u = max(D1_low.norm(), D1_up.norm()) * max(D1_low.inverse().norm(), D1_up.inverse().norm());

    vector2 v2, v1;

    v1 = phi_(F, start);
    v2 = phi(F, v1);

    int i = 1;
    do
    {
        v1 = phi_(F, v2);
        v2 = phi(F, v1);
        i++;
        // cout << i << " " << v2.x << " " << v2.y << " " << F(v2).norm() << " " << q << endl;
    } while ((v1 - v2).norm() > (1 - q) / q * eps);
    cout << "Метод Зейделя, " << v2.x << ", " << v2.y << ", " << F(v2).norm() << ", " << start.x << ", " << start.y << ", "
         << i << ", " << q << ", " << u << endl;
}

void Newtones_method(vector2 (*F)(vector2), vector2 low, vector2 up, vector2 start)
{
    Dmatrix D1_low(df1dx(low), df1dy(low), dg1dx(low), dg1dy(low));
    Dmatrix D1_up(df1dx(up), df1dy(up), dg1dx(up), dg1dy(up));
    double q = max(D1_low.norm(), D1_up.norm());
    double u = max(D1_low.norm(), D1_up.norm()) * max(D1_low.inverse().norm(), D1_up.inverse().norm());

    vector2 v2(start.x, start.y), v1;
    int i = 0;
    do
    {
        v1 = v2;
        Dmatrix D(dfdx(v1), dfdy(v1), dgdx(v1), dgdy(v1));
        v2 = v1 - D.inverse() * F(v1);
        i++;
        // cout << i << " " << v2.x << " " << v2.y << " " << F(v2).norm() << " " << q << endl;
    } while ((v1 - v2).norm() > eps / u);
    cout << "Метод Ньютона, " << v2.x << ", " << v2.y << ", " << F(v2).norm() << ", " << start.x << ", " << start.y << ", "
        << i << ", " << q << ", " << u << endl;
}

int main()
{
    vector2 low(-0.5, 1.5), up(-1, 2), start(-1, 1);

    cout << ", x, y, норма невязки ||f(x; y)||, x0, y0, количество итераций  N+1, q, m" << endl;

    SimpleIteration(F, low, up, start);
    Zeidels_method(F, low, up, start);
    Newtones_method(F, low, up, start);
    cout << endl;

    low = vector2(-0.5, -2);
    up = vector2(-1, -1.5);
    start = vector2(-1, -1.5);

    SimpleIteration(F, low, up, start);
    Zeidels_method(F, low, up, start);
    Newtones_method(F, low, up, start);
    cout << endl;

    return 0;
}