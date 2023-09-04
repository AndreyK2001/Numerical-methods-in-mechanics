#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001  // Требуемая точность вычисления
#define inf 100000 // Требуемая точность вычисления

double sign(double x)
{
    if (x == 0)
        return 0;
    else
        return x / fabs(x);
}

typedef struct vec
{
    double x;
    double y;
} vec;

double f(vec v)
{
    return 0.6 * v.x + 7.5 * v.y + v.x * v.x * v.y;
}

double dfdx(vec v)
{
    return 0.6 + 2 * v.x * v.y;
}

double dfdy(vec v)
{
    return 7.5 + v.x * v.x;
}

double g(vec v)
{
    return 6 * v.x + cos(v.y);
}

double dgdx(vec v)
{
    return 6;
}

double dgdy(vec v)
{
    return -sin(v.y);
}

double f1(vec v)
{
    return -cos(v.y) / 6;
}

double df1dx(vec v)
{
    return 0;
}

double df1dy(vec v)
{
    return sin(v.y) / 6;
}

double g1(vec v)
{
    return -0.6 * v.x / (7.5 + v.x * v.x);
}

double dg1dx(vec v)
{
    return -(0.6 * (7.5 + v.x * v.x) - 1.2 * v.x * v.x) / ((7.5 + v.x * v.x) * (7.5 + v.x * v.x));
}

double dg1dy(vec v)
{
    return 0;
}

typedef struct Matrix4
{
    double D[2][2];
} Matrix4;

Matrix4 D(vec v)
{
    Matrix4 matrix = {.D = {{dfdx(v), dfdy(v)}, {dgdx(v), dgdy(v)}}};
    return matrix;
};

Matrix4 invD(vec v)
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

    Matrix4 matrix = {.D = {{d, -b}, {-c, a}}};
    return matrix;
};

double norm(Matrix4 matrix)
{
    return max(fabs(matrix.D[0][0]) + fabs(matrix.D[0][1]), fabs(matrix.D[1][0]) + fabs(matrix.D[1][1]));
}

double norm(vec v)
{
    return max(fabs(v.x), fabs(v.y));
}

vec F1(vec v)
{
    vec res = {.x = f1(v), .y = g1(v)};
    return res;
}

vec F(vec v)
{
    vec res = {.x = f(v), .y = g(v)};
    return res;
}

vec operator+(vec a, vec b)
{
    return vec{.x = a.x + b.x, .y = a.y + b.y};
}

vec operator-(vec a, vec b)
{
    return vec{.x = a.x - b.x, .y = a.y - b.y};
}

vec operator*(Matrix4 b, vec a)
{
    double x = b.D[0][0] * a.x + b.D[0][1] * a.y;
    double y = b.D[1][0] * a.x + b.D[1][1] * a.y;
    return vec{.x = x, .y = y};
}

vec operator/(vec v, double b)
{
    return vec{.x = v.x / b, .y = v.y / b};
}

vec operator*(vec v, double b)
{
    return vec{.x = v.x * b, .y = v.y * b};
}

typedef struct ans
{
    vec root;
    double normF;
    vec v0;
    int n_iterations;
    double q;
    double mu;
} ans;

Matrix4 sign(Matrix4 d)
{
    Matrix4 matrix = {.D = {{sign(d.D[0][0]), sign(d.D[0][1])}, {sign(d.D[1][0]), sign(d.D[1][1])}}};
    return matrix;
};

vec phi(vec (*F)(vec), vec v, double q)
{
    double sx = sign((invD(v) * F(v)).x);
    double sy = sign((invD(v) * (v)).y);
    vec delta{.x = F(v).x * sx, .y = F(v).y * sy};

    return v - sign(D(v)) * F(v) * q;
}
Matrix4 D_(vec v)
{
    Matrix4 matrix = {.D = {{df1dx(v), df1dy(v)}, {dg1dx(v), dg1dy(v)}}};
    return matrix;
};

ans SimpleIteration_method(vec (*F)(vec), vec left, vec right, vec v_0) // Метод простой итерации
{
    double q = max(norm(D_(left)), norm(D_(right)));
    vec v_new, v;
    v_new = v_0;
    int i = 0;
    do
    {
        v = v_new;
        v_new = phi(F, v, q);
        i++;
        // cout << i << " " << v_new.x << " " << v_new.y << endl;
    } while (norm(v - v_new) > eps);

    ans params = {.root = v_new, .normF = norm(F(v_new)), .v0 = v_0, .n_iterations = i, .q = q, .mu = norm(D(v_new)) * norm(invD(v_new))};
    return params;
}

vec phi1(vec (*F)(vec), vec v, double q)
{
    double sx = sign((invD(v) * F(v)).x);
    double sy = sign((invD(v) * (v)).y);
    vec delta{.x = F(v).x * sx, .y = F(v).y * sy};
    vec v_res = v - sign(D(v)) * F(v) * q;

    return vec{.x = v_res.x, .y = v.y};
}

ans Zeydels_method(vec (*F)(vec), vec left, vec right, vec v_0) // Метод Зейделя
{
    double q = max(norm(D_(left)), norm(D_(right)));
    vec v1, v2;
    v1 = phi1(F, v_0, q);
    v2 = phi(F, v1, q);
    int i = 1;
    do
    {
        v1 = phi1(F, v2, q);
        v2 = phi(F, v1, q);
        i++;
        // cout << i << " " << v2.x << " " << v2.y << endl;
    } while (norm(v1 - v2) > eps);

    ans params = {.root = v2, .normF = norm(F(v2)), .v0 = v_0, .n_iterations = i, .q = q, .mu = norm(D(v2)) * norm(invD(v2))};
    return params;
}

ans Newtones_method(vec (*F)(vec), vec left, vec right, vec v_0) // Метод Ньютона (метод касательных)
{
    vec v1, v2 = v_0;
    double q = max(norm(D_(left)), norm(D_(right)));
    double mu = norm(D(v2)) * norm(invD(v2));
    int i = 0;
    do
    {
        v1 = v2;
        i++;
        v2 = v2 - invD(v2) * F(v2);
        // cout << i << " " << v2.x << " " << v2.y << endl;
    } while (norm(v1 - v2) > eps / mu);

    ans params = {.root = v2, .normF = norm(F(v2)), .v0 = v_0, .n_iterations = i, .q = q, .mu = mu};
    return params;
}

int main()
{
    ofstream fout; // Открываем поток

    fout.open("otchet.tex");

    // Преамбула

    fout << "\\documentclass[a4paper]{article}" << endl;
    fout << "\\usepackage[12pt]{extsizes}" << endl;
    fout << "\\usepackage[utf8]{inputenc}" << endl;
    fout << "\\usepackage[russian]{babel}" << endl;
    fout << "\\usepackage{setspace,amsmath}" << endl;
    fout << "\\usepackage[left=20mm, top=20mm, right=20mm, bottom=20mm, nohead, footskip=10mm, landscape]{geometry}" << endl;
    fout << "\\usepackage{indentfirst}" << endl;
    fout << "\\usepackage{multirow}" << endl;
    fout << "\\usepackage{fourier} " << endl;
    fout << "\\usepackage{array}" << endl;
    fout << "\\usepackage{makecell}" << endl;
    fout << "\n"
         << endl;

    fout << "\\newcolumntype{P}[1]{>{\\centering\\arraybackslash}p{#1}} \n"
         << endl;
    fout << "\\newcolumntype{M}[1]{>{\\centering\\arraybackslash}m{#1}} \n"
         << endl;

    // fout << "\\renewcommand\\theadalign{bc}" << endl;
    // fout << "\\renewcommand\\theadfont{\\bfseries}" << endl;
    // fout << "\\renewcommand\\theadgape{\\Gape[4pt]}" << endl;
    // fout << "\\renewcommand\\cellgape{\\Gape[4pt]}" << endl;

    fout << "\\begin{document}" << endl;

    fout << "\\renewcommand{\\arraystretch}{2}" << endl;
    fout << "\\begin{center}" << endl;
    fout << "{\\Large Решение нелинейных алгебраических систем}" << endl;
    fout << "\\bigskip" << endl;
    fout << "\\end{center}" << endl;
    fout << "\n"
         << endl;

    fout << "Выполнил: Капелюшников Андрей Сергеевич"
         << "\n"
         << endl;

    fout << "\\begin{equation*}" << endl;
    fout << "\\begin{cases}" << endl;
    fout << "0.6x+7.5y+x^2y=0" << endl;
    fout << "\\\\" << endl;
    fout << "6x+\\cos(y)=0" << endl;
    fout << "\\end{cases}" << endl;
    fout << "\\end{equation*}" << endl;
    fout << "\n"
         << endl;

    fout << "$\\varepsilon = $" << eps << endl;
    fout << "\n"
         << endl;

    //   fout << "\\pagebreak" << endl;

    fout << "\\begin{table}[h!]" << endl;
    fout << "\\centering" << endl;
    fout << "\\begin{tabular}{P{2cm} ccccccccc}" << endl;
    fout << "\\hline" << endl;
    fout << "\\multicolumn{2}{c}{Пара корней $\\vec{\\xi} ^* _r$} & \\multirow{2}{*}{\\thead{Норма невязки \\\\ $||\\vec{f}(\\xi ^* _r)||$}} & \\multicolumn{2}{c}{Начальный вектор} & \\multicolumn{3}{c}{Число итераций $N+1$} & \\multirow{2}{*}{$q_r$} & \\multirow{2}{*}{$\\mu _r$} \\\\[2em] \\cline{1-2} \\cline{4-8}" << endl;
    fout << " $\\xi ^* _r$ & $\\eta ^* _r$ &                   &            $x_0 ^{(r)}$     &    $y_0 ^{(r)}$   &  \\thead{Простая \\\\ итерация}  & \\thead{Метод \\\\ Зейделя}  &   \\thead{Метод \\\\ Ньютона}   &                   &                   \\\\ \\hline" << endl;

    vec left{.x = -0.2, .y = 0.01}, right{.x = -0.15, .y = 0.02}, v_0{.x = -0.155, .y = 0.01};
    ans simple = SimpleIteration_method(F, left, right, v_0);
    fout << fixed << setprecision(6) << simple.root.x << " & " << simple.root.y << " & " << simple.normF << " & " << simple.v0.x << "&" << simple.v0.y << "&" << simple.n_iterations << "&       &      &" << simple.q << "&" << simple.mu << "\\\\ \\hline" << endl;

    ans zeyd = Zeydels_method(F, left, right, v_0);
    fout << fixed << setprecision(6) << zeyd.root.x << " & " << zeyd.root.y << " & " << zeyd.normF << " & " << zeyd.v0.x << "&" << zeyd.v0.y << "& & " << zeyd.n_iterations << "& &" << zeyd.q << "&" << zeyd.mu << "\\\\ \\hline" << endl;

    ans newt = Newtones_method(F, left, right, v_0);
    fout << fixed << setprecision(6) << newt.root.x << " & " << newt.root.y << " & " << newt.normF << " & " << newt.v0.x << "&" << newt.v0.y << "& & & " << newt.n_iterations << " & " << newt.q << "&" << newt.mu << "\\\\ \\hline" << endl;

    fout << "\\end{tabular}" << endl;
    fout << "\\end{table}" << endl;
    fout << "\\bigskip" << endl;

    fout << "\\end{document}" << endl;

    // для формирования отчёта необходимо запустить в компиляторе LaTex файл otchet.tex
    return 0;
}