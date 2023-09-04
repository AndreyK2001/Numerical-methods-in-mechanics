#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001 // Требуемая точность вычисления

double sign(double x)
{
    if (x == 0)
        return 0;
    else if (x > 0)
        return 1;
    else
        return -1;
}

typedef struct _vector
{ // вектор
    double x;
    double y;
} _vector;

double f(_vector v) { return pow(v.x + 1, 2) + pow(v.y + 1, 2) - 1; }

double dfdx(_vector v) { return 2 * (v.x + 1); }

double dfdy(_vector v) { return 2 * (v.y + 1); }

double g(_vector v) { return sqrt(v.x + 1) - v.y - 1; }

double dgdx(_vector v) { return 1 / (2 * sqrt(v.x + 1)); }

double dgdy(_vector v) { return -1; }

double f1(_vector v) { return sqrt(1 - pow(v.y + 1, 2)) - 1; }

double df1dx(_vector v) { return 0; }

double df1dy(_vector v) { return -(v.y + 1) / sqrt(1 - pow(v.y + 1, 2)); }

double g1(_vector v) { return sqrt(v.x + 1) - 1; }

double dg1dx(_vector v) { return 1 / (2 * sqrt(v.x + 1)); }

double dg1dy(_vector v) { return 0; }

typedef struct Matrix
{
    double D[2][2];
} Matrix;

Matrix D(_vector v)
{
    Matrix matrix = {.D = {{dfdx(v), dfdy(v)}, {dgdx(v), dgdy(v)}}};
    return matrix;
};

Matrix invD(_vector v)
{
    double a = dfdx(v);
    double b = dfdy(v);
    double c = dgdx(v);
    double d = dgdy(v);
    double det = (a * d - b * c);
    a /= det;
    b /= det;
    c /= det;
    d /= det;

    Matrix matrix = {.D = {{d, -b}, {-c, a}}};
    return matrix;
};

double norm(Matrix matrix)
{
    return max(fabs(matrix.D[0][0]) + fabs(matrix.D[0][1]), fabs(matrix.D[1][0]) + fabs(matrix.D[1][1]));
}

double norm(_vector v)
{
    return max(fabs(v.x), fabs(v.y));
}

_vector F1(_vector v)
{
    _vector res = {.x = f1(v), .y = g1(v)};
    return res;
}

_vector F(_vector v)
{
    _vector res = {.x = f(v), .y = g(v)};
    return res;
}

_vector operator+(_vector a, _vector b)
{
    return _vector{.x = a.x + b.x, .y = a.y + b.y};
}

_vector operator-(_vector a, _vector b)
{
    return _vector{.x = a.x - b.x, .y = a.y - b.y};
}

_vector operator*(Matrix b, _vector a)
{
    double x = b.D[0][0] * a.x + b.D[0][1] * a.y;
    double y = b.D[1][0] * a.x + b.D[1][1] * a.y;
    return _vector{.x = x, .y = y};
}

_vector operator/(_vector v, double b)
{
    return _vector{.x = v.x / b, .y = v.y / b};
}

_vector operator*(_vector v, double b)
{
    return _vector{.x = v.x * b, .y = v.y * b};
}

typedef struct result
{
    _vector root;
    double normF;
    _vector v0;
    int n_iterations;
    double q;
    double mu;
} result;

Matrix sign(Matrix d)
{
    Matrix matrix = {.D = {{sign(d.D[0][0]), sign(d.D[0][1])}, {sign(d.D[1][0]), sign(d.D[1][1])}}};
    return matrix;
};

_vector phi(_vector (*F)(_vector), _vector v, double q)
{
    double sx = sign((invD(v) * F(v)).x);
    double sy = sign((invD(v) * (v)).y);
    _vector delta{.x = F(v).x * sx, .y = F(v).y * sy};

    return v - invD(v) * F(v) / norm(D(v));
}

_vector phi(_vector (*F)(_vector), _vector v)
{
    double sx = sign((invD(v) * F(v)).x);
    double sy = sign((invD(v) * (v)).y);
    _vector delta{.x = F(v).x * sx, .y = F(v).y * sy};

    return v - invD(v) * F(v);
}

Matrix D_(_vector v)
{
    Matrix matrix = {.D = {{df1dx(v), df1dy(v)}, {dg1dx(v), dg1dy(v)}}};
    return matrix;
};

result Simple_iteration(_vector (*F)(_vector), _vector left, _vector right, _vector v_0) // Метод простой итерации
{
    double q = max(norm(D_(left)), norm(D_(right))) / 10;
    _vector v_new, v;
    v_new = v_0;
    int i = 0;
    do
    {
        v = v_new;
        v_new = phi(F, v, q);
        i++;

    } while (norm(v - v_new) > eps);

    result params = {.root = v_new, .normF = norm(F(v_new)), .v0 = v_0, .n_iterations = i, .q = q, .mu = norm(D(v_new)) * norm(invD(v_new))};
    return params;
}

_vector phi1(_vector (*F)(_vector), _vector v, double q)
{
    double sx = sign((invD(v) * F(v)).x);
    double sy = sign((invD(v) * (v)).y);
    _vector delta{.x = F(v).x * sx, .y = F(v).y * sy};
    _vector v_res = v - sign(D(v)) * F(v);

    return _vector{.x = v_res.x, .y = v.y};
}

result Zeydal(_vector (*F)(_vector), _vector left, _vector right, _vector v_0) // Метод Зейделя
{
    double q = max(norm(D_(left)), norm(D_(right))) / 10;
    _vector v1, v2;
    v1 = phi1(F, v_0, q);
    v2 = phi(F, v1);
    int i = 1;
    do
    {
        v1 = phi1(F, v2, q);
        v2 = phi(F, v1);
        i++;
        // cout << i << " " << v2.x << " " << v2.y << endl;
    } while (norm(v1 - v2) > eps);

    result params = {.root = v2, .normF = norm(F(v2)), .v0 = v_0, .n_iterations = i, .q = q, .mu = norm(D(v2)) * norm(invD(v2))};
    return params;
}

result Newtones(_vector (*F)(_vector), _vector left, _vector right, _vector v_0) // Метод Ньютона (метод касательных)
{
    _vector v1, v2 = v_0;
    double q = max(norm(D_(left)), norm(D_(right)));
    double mu = norm(D(v2)) * norm(invD(v2));
    int i = 0;
    do
    {
        v1 = v2;
        i++;
        v2 = v2 - invD(v2) * F(v2);

    } while (norm(v1 - v2) > eps / mu);

    result params = {.root = v2, .normF = norm(F(v2)), .v0 = v_0, .n_iterations = i, .q = q, .mu = mu};
    return params;
}

int main()
{
    FILE *pFile;
    pFile = fopen("result.csv", "w");

    fprintf(pFile, "%s \n", ", x, y, ||f(x y)||, x_0, y_0, N+1, q, m, \n");

    _vector left_bottom_corner{.x = -0.5, .y = -0.5}, right_top_corner{.x = -0.1, .y = -0.1}, v_0{.x = -0.5, .y = -0.5};

    result simple = Simple_iteration(F, left_bottom_corner, right_top_corner, v_0);
    fprintf(pFile, "Простая итерация, %g, %g, %g, %g, %g, %d, %g,  \n", simple.root.x, simple.root.y, simple.normF, simple.v0.x, simple.v0.y, simple.n_iterations, simple.q, simple.mu);

    result zeyd = Zeydal(F, left_bottom_corner, right_top_corner, v_0);
    fprintf(pFile, "Метод Зейделя, %g, %g, %g, %g, %g, %d, %g,  \n", zeyd.root.x, zeyd.root.y, zeyd.normF, zeyd.v0.x, zeyd.v0.y, zeyd.n_iterations, zeyd.q, zeyd.mu);

    result newt = Newtones(F, left_bottom_corner, right_top_corner, v_0);
    fprintf(pFile, "Метод Ньютона, %g, %g, %g, %g, %g, %d, , %g \n", newt.root.x, newt.root.y, newt.normF, newt.v0.x, newt.v0.y, newt.n_iterations, newt.q, newt.mu);

    fclose(pFile);
    return 0;
}