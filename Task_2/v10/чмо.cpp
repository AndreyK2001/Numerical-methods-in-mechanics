#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

/**
 * f1(x, y) = x + y^(2/3)
 * f2(x, y) = x^2 + 2y^2 - 1
 * phi1(x,y) = -y^(2/3) = x
 * phi2(x,y) = sqrt((1 - x^2)/2) = y
 * J =
 *  1 2/3y^(1/3)
 *  2x 4y
 *  J^-1 =
 *  1/(x-3y^(4/3)) *
 *  -3y^(4/3)  1/2
 *  3xy^(1/3)/2  -3y^(1/3) / 4
 */

using namespace std;

double norm(double x, double y)
{
    return sqrt(x * x + y * y);
}

double *simpleIteration(double x0, double y0, double eps)
{
    double *ans = new double[6];
    double x = x0, y = y0;
    double xk = -pow(fabs(y), 2. / 3);
    double yk = -sqrt((1 - x * x) / 2);
    int num_iter = 1;
    while (norm(xk - x, yk - y) > eps)
    {
        double tempx = xk, tempy = yk;
        xk = -pow(fabs(tempy), 2. / 3);
        yk = -sqrt((1 - tempx * tempx) / 2);
        x = tempx;
        y = tempy;
        num_iter += 1;
    }

    ans[0] = xk;
    ans[1] = yk;
    ans[2] = norm(xk + pow(yk * yk, 1. / 3), xk * xk + 2 * yk * yk - 1);
    ans[3] = x0;
    ans[4] = y0;
    ans[5] = num_iter;

    // printf("simple iteration num of iterations: %d\n", num_iter);
    return ans;
}

double *zeidelIreration(double x0, double y0, double eps)
{
    double *ans = new double[6];
    double x = x0, y = y0;
    double xk = -pow(fabs(y), 2. / 3);
    double yk = -sqrt((1 - xk * xk) / 2);
    int num_iter = 1;
    while (norm(xk - x, yk - y) > eps)
    {
        double tempx = xk, tempy = yk;
        xk = -pow(fabs(tempy), 2. / 3);
        yk = -sqrt((1 - xk * xk) / 2);
        x = tempx;
        y = tempy;
        num_iter += 1;
    }
    ans[0] = xk;
    ans[1] = yk;
    ans[2] = norm(xk + pow(yk * yk, 1. / 3), xk * xk + 2 * yk * yk - 1);
    ans[3] = x0;
    ans[4] = y0;
    ans[5] = num_iter;

    // printf("Zeidel's method num of iteration: %d\n", num_iter);
    return ans;
}

double JinverseX(double x0, double y0, double xk, double yk)
{
    double y43 = pow(y0 * y0, 2. / 3);
    return (xk * (-3 * y43) + 0.5 * yk) / (x0 - 3 * y43);
}

double JinverseY(double x0, double y0, double xk, double yk)
{
    double y43 = pow(y0 * y0, 2. / 3);
    double y13 = -pow(fabs(y0), 1. / 3);
    return (3 * x0 * y13 / 2 * xk - 3 * y13 / 4 * yk) / (x0 - 3 * y43);
}

double fx(double x, double y)
{
    return x + pow(y * y, 1. / 3);
}

double fy(double x, double y)
{
    return x * x + 2 * y * y - 1;
}

double *newtonIteration(double x0, double y0, double eps)
{
    double *ans = new double[2];
    double x = x0, y = y0;
    double xk, yk;
    double fkx = fx(x0, y0), fky = fy(x0, y0);
    xk = x0 - JinverseX(x0, y0, fkx, fky), yk = fky - JinverseY(x0, y0, fkx, fky);
    int num_iter = 1;
    while (norm(xk - x, yk - y) > eps)
    {
        double tempx = xk, tempy = yk;
        // printf("%lf %lf\n", xk, yk);
        fkx = fx(xk, yk), fky = fy(xk, yk);
        xk -= JinverseX(tempx, tempy, fkx, fky);
        yk -= JinverseY(tempx, tempy, fkx, fky);
        x = tempx;
        y = tempy;
        num_iter += 1;
    }

    ans[0] = xk;
    ans[1] = yk;
    ans[2] = norm(xk + pow(yk * yk, 1. / 3), xk * xk + 2 * yk * yk - 1);
    ans[3] = x0;
    ans[4] = y0;
    ans[5] = num_iter;

    // printf("Newton's method num of iterations: %d\n", num_iter);
    return ans;
}

int main(void)
{
    ofstream fout; // Открываем поток
    fout.open("otchet.csv");

    fout << "x, y, Норма невязки, x0, y0, Простая итерация, Метод Зейделя, Метод Ньютона" << endl;
    double *ans = simpleIteration(0.0, -0.2, 1e-4);
    fout << ans[0] << ", " << ans[1] << ", " << ans[2] << ", " << ans[3] << ", " << ans[4] << ", " << ans[5] << ", , " << endl;
    free(ans);
    ans = zeidelIreration(0.0, -0.2, 1e-4);
    fout << ans[0] << ", " << ans[1] << ", " << ans[2] << ", " << ans[3] << ", " << ans[4] << ", , " << ans[5] << ", " << endl;
    free(ans);
    ans = newtonIteration(0.0, -0.2, 1e-4);
    fout << ans[0] << ", " << ans[1] << ", " << ans[2] << ", " << ans[3] << ", " << ans[4] << ", , , " << ans[5] << endl;
    free(ans);

    fout.close();
}
