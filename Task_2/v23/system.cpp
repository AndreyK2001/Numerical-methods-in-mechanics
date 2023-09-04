#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 1e-4 // Требуемая точность вычисления

typedef struct xy
{ //для удобства работы с векторами
    double x;
    double y;
} xy;

xy operator+(xy a, xy b)
{
    return xy{.x = a.x + b.x, .y = a.y + b.y};
}

double f(xy v)
{ //f(x, y)=0
    return v.x / (1 + v.x * v.x) - v.y;
}

double g(xy v)
{ //g(x, y)=0
    return v.x * v.x + v.y * v.y - 1;
}

double dfdx(xy v)
{
    return (1 - v.x * v.x) / pow(1 + v.x * v.x, 2);
}

double dgdx(xy v)
{
    return 2 * v.x;
}

double dfdy(xy v)
{
    return -1;
}

double dgdy(xy v)
{
    return 2 * v.y;
}

double ix(xy v)
{ //x=phi(y)
    return sqrt(1 - v.y * v.y);
}

double dixdx(xy v)
{
    return 0;
}

double dixdy(xy v)
{
    return -v.y / sqrt(1 - v.y * v.y);
}

double iy(xy v)
{//y=ksi(x)
    return v.x / (1 + v.x * v.x);
}

double diydx(xy v)
{
    return (1 - v.x * v.x) / pow(1 + v.x * v.x, 2);
}

double diydy(xy v)
{
    return 0;
}

typedef struct Matrix
{
    double D[2][2];
} Matrix;

Matrix D(xy v)
{
    Matrix matrix = {.D = {{dfdx(v), dfdy(v)}, {dgdx(v), dgdy(v)}}};
    return matrix;
};
xy operator*(Matrix b, xy a)
{
    double x = b.D[0][0] * a.x + b.D[0][1] * a.y;
    double y = b.D[1][0] * a.x + b.D[1][1] * a.y;
    return xy{.x = x, .y = y};
}

double sign(double x)
{
    if (x == 0)
        return 0;
    else if (x > 0)
        return 1;
    else
        return -1;
}

Matrix invD(xy v)
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

double norm(xy v)
{
    return max(fabs(v.x), fabs(v.y));
}

xy F1(xy v)
{
    xy res = {.x = ix(v), .y = iy(v)};
    return res;
}

xy F(xy v)
{
    xy res = {.x = f(v), .y = g(v)};
    return res;
}

xy operator-(xy a, xy b)
{
    return xy{.x = a.x - b.x, .y = a.y - b.y};
}

xy operator/(xy v, double b)
{
    return xy{.x = v.x / b, .y = v.y / b};
}

xy operator*(xy v, double b)
{
    return xy{.x = v.x * b, .y = v.y * b};
}

typedef struct result
{
    xy root;
    double normF;
    xy v0;
    int n_iterations;
    double q;
    double mu;
} result;

Matrix sign(Matrix d)
{
    Matrix matrix = {.D = {{sign(d.D[0][0]), sign(d.D[0][1])}, {sign(d.D[1][0]), sign(d.D[1][1])}}};
    return matrix;
};

xy phi(xy (*F)(xy), xy v, double q)
{
    double sx = sign((invD(v) * F(v)).x);
    double sy = sign((invD(v) * (v)).y);
    xy delta{.x = F(v).x * sx, .y = F(v).y * sy};

    return v - invD(v) * F(v) / norm(D(v));
}

xy phi(xy (*F)(xy), xy v)
{
    double sx = sign((invD(v) * F(v)).x);
    double sy = sign((invD(v) * (v)).y);
    xy delta{.x = F(v).x * sx, .y = F(v).y * sy};

    return v - invD(v) * F(v);
}

Matrix D_(xy v)
{
    Matrix matrix = {.D = {{dixdx(v), dixdy(v)}, {diydx(v), diydy(v)}}};
    return matrix;
};

result Iteration(xy (*F)(xy), xy left, xy right, xy X0) // Метод простой итерации
{
    double q = max(norm(D_(left)), norm(D_(right))) / 10;
    xy v_new, v;
    v_new = X0;
    int i = 0;
    do
    {
        v = v_new;
        v_new = phi(F, v, q);
        i++;

    } while (norm(v - v_new) > eps);

    result params = {.root = v_new, .normF = norm(F(v_new)), .v0 = X0, .n_iterations = i, .q = q, .mu = norm(D(v_new)) * norm(invD(v_new))};
    return params;
}

xy phi1(xy (*F)(xy), xy v, double q)
{
    double sx = sign((invD(v) * F(v)).x);
    double sy = sign((invD(v) * (v)).y);
    xy delta{.x = F(v).x * sx, .y = F(v).y * sy};
    xy v_res = v - sign(D(v)) * F(v);

    return xy{.x = v_res.x, .y = v.y};
}

result Zeydel(xy (*F)(xy), xy left, xy right, xy X0) // Метод Зейделя
{
    double q = max(norm(D_(left)), norm(D_(right))) / 10;
    xy v1, v2;
    v1 = phi1(F, X0, q);
    v2 = phi(F, v1);
    int i = 1;
    do
    {
        v1 = phi1(F, v2, q);
        v2 = phi(F, v1);
        i++;
        // cout << i << " " << v2.x << " " << v2.y << endl;
    } while (norm(v1 - v2) > eps);

    result params = {.root = v2, .normF = norm(F(v2)), .v0 = X0, .n_iterations = i, .q = q, .mu = norm(D(v2)) * norm(invD(v2))};
    return params;
}

result Newton(xy (*F)(xy), xy left, xy right, xy X0) // Метод Ньютона (метод касательных)
{
    xy v1, v2 = X0;
    double q = max(norm(D_(left)), norm(D_(right)));
    double mu = norm(D(v2)) * norm(invD(v2));
    int i = 0;
    do
    {
        v1 = v2;
        i++;
        v2 = v2 - invD(v2) * F(v2);

    } while (norm(v1 - v2) > eps / mu);

    result params = {.root = v2, .normF = norm(F(v2)), .v0 = X0, .n_iterations = i, .q = q, .mu = mu};
    return params;
}

int main()
{
    FILE *pFile;
    pFile = fopen("result.csv", "w");

    fprintf(pFile, "%s \n", ", x, y, ||f(x y)||, x0, y0, количество итераций, q, m, \n");

    xy A{.x = 0.8, .y = 0.4}, B{.x = 0.9, .y = 0.6}, X0{.x = 0.85, .y = 0.45};
    result S = Iteration(F, A, B, X0);
    fprintf(pFile, "Простая итерация, %g, %g, %g, %g, %g, %d, %g,  \n", S.root.x, S.root.y, S.normF, S.v0.x, S.v0.y, S.n_iterations, S.q, S.mu);
    result Z = Zeydel(F, A, B, X0);
    fprintf(pFile, "Метод Зейделя, %g, %g, %g, %g, %g, %d, %g,  \n", Z.root.x, Z.root.y, Z.normF, Z.v0.x, Z.v0.y, Z.n_iterations, Z.q, Z.mu);
    result N = Newton(F, A, B, X0);
    fprintf(pFile, "Метод Ньютона, %g, %g, %g, %g, %g, %d, , %g \n", N.root.x, N.root.y, N.normF, N.v0.x, N.v0.y, N.n_iterations, N.q, N.mu);

    fclose(pFile);
    return 0;
}