#include <iostream> // Стандартная библиотека ввода-вывода
#include <math.h>   // Библиотека станлартных математических функций
#include <fstream>  // Библиотека функций для работы с потоками (файлами)

using namespace std;

#define eps 0.0001  // Требуемая точность вычисления
#define inf 1000000 // Требуемая точность вычисления

double sign(double x)
{
    if (x == 0)
        return 0;
    else
        return x / fabs(x);
}

double DiracDelta(double x)
{
    if (x == 0)
        return inf;
    else
        return 0;
}

double f(double x, double a)
{
    return cos(exp(fabs(sin(x)))) - a * x;
}

double dfdx(double x, double a)
{
    return -a - exp(fabs(sin(x))) * sin(exp(fabs(sin(x)))) * cos(x) * sign(sin(x));
}

double d2fdx2(double x, double a)
{
    return (-(-exp(fabs(sin(x))) * sin(x) + pow(cos(x), 2) * exp(fabs(sin(x))) * sign(sin(x))) * sign(sin(x)) -
            2 * pow(cos(x), 2) * DiracDelta(sin(x)) * exp(fabs(sin(x)))) *
               sin(exp(fabs(sin(x)))) -
           pow(cos(x), 2) * pow(sign(sin(x)), 2) * cos(exp(fabs(sin(x)))) * exp(2 * fabs(sin(x)));
}

typedef struct ans
{
    double root;
    double f;
    double start;
    int n_iterations;
    double M;
    double m;
    double q;
} ans;

ans Chord_method(double (*f)(double, double), double left, double right, double x_0, double a) // Метод хорд (метод секущих)
{
    double x, x_new = x_0;
    double d2 = d2fdx2(x_0, a); // вторая производная

    int i = 0; // счётчик итераций алгоритма

    double m1 = min(fabs(dfdx(left, a)), fabs(dfdx(right, a)));
    double M1 = max(fabs(dfdx(left, a)), fabs(dfdx(right, a)));

    double m = (0 + m1) / 2;
    double M = M1 + 2 * eps;

    // cout << eps * m / (M - m);

    if (f(left, a) * d2fdx2(x_0, a) > 0)
    {
        do
        {
            x = x_new;
            i++;
            x_new = x - f(x, a) / (f(x, a) - f(left, a)) * (x - left);
        } while (fabs(x_new - x) > eps * m / (M - m));
    }
    else
    {
        do
        {
            x = x_new;
            i++;
            x_new = x - f(x, a) / (f(right, a) - f(x, a)) * (right - x);
            // cout << i << " " << x_new << endl;
        } while (fabs(x_new - x) > eps * m / (M - m));
    }

    ans params = {.root = x_new, .f = f(x_new, a), .start = x_0, .n_iterations = i, .M = M, .m = m};
    return params;
}

ans Newtones_method(double (*f)(double, double), double left, double right, double x_0, double a) // Метод Ньютона (метод касательных)
{
    int i = 0; // Число витков цикла n, счетчик i
    double x = x_0, x_new = x_0;

    double m1 = min(fabs(dfdx(left, a)), fabs(dfdx(right, a)));
    double M1 = max(fabs(dfdx(left, a)), fabs(dfdx(right, a)));

    double m = (0 + m1) / 2;
    double M = M1 + 2 * eps;

    do
    {
        x = x_new;
        i++;
        x_new = x - f(x, a) / dfdx(x, a);
        // cout << i << " " << x_new << endl;
    } while (fabs(x_new - x) > fabs(eps * m / (M - m)));

    ans params = {.root = x_new, .f = f(x_new, a), .start = x_0, .n_iterations = i, .M = M, .m = m};
    return params;
}

double phi(double (*f)(double, double), double left, double right, double x, double a)
{
    double M = max(fabs(dfdx(left, a)), fabs(dfdx(right, a))) + 2 * eps;
    return x - sign(dfdx(x, a)) * f(x, a) / M;
}

double dphidx(double (*f)(double, double), double left, double right, double x, double a)
{
    double M = max(fabs(dfdx(left, a)), fabs(dfdx(right, a))) + 2 * eps;
    return 1 - DiracDelta(dfdx(x, a)) * d2fdx2(x, a) * f(x, a) / M - sign(dfdx(x, a)) * dfdx(x, a) / M;
}

ans SimpleIteration_method(double (*f)(double, double), double left, double right, double x_0, double a) // Метод простой итерации
{
    double x = x_0, x_new = x_0;
    int i = 0; // Число витков цикла n, счетчик i

    double M1 = max(fabs(dphidx(f, left, right, left, a)), fabs(dphidx(f, left, right, right, a))) + eps;
    double q = (1 + M1) / 2;

    do
    {
        x = x_new;
        i++;
        x_new = phi(f, left, right, x, a);
        // cout << i << " " << x_new << endl;
    } while (fabs(x_new - x) > fabs(eps * (1 - q) / q));

    ans params = {.root = x_new, .f = f(x_new, a), .start = x_0, .n_iterations = i, .q = q};
    return params;
}

ans half_method(double (*f)(double, double), double left, double right, double x_0, double a) // Метод деления отрезка пополам (метод половинного деления)
{
    int n, i;                                             // Число витков цикла n, счетчик i
    double x = (left + right) / 2, a_ = left, b_ = right; // Вычисляем середину начального отрезка

    n = int(log((b_ - a_) / eps) / log(2)) + 1;

    for (i = 1; i <= n && fabs(f(x, a)) > eps; i++)
    {
        if (f(x, a) * f(a_, a) < 0)
            b_ = x;
        else
            a_ = x;

        x = (a_ + b_) / 2.; // Вычисляем середину нового отрезка
    }

    ans params = {.root = x, .f = f(x, a), .start = x_0, .n_iterations = i};
    return params;
}

ans Etkins_method(double (*f)(double, double), double left, double right, double x_0, double a) // Метод Ньютона (метод касательных)
{
    int i = 1; // Счетчик i
    double x_n = x_0, x_n_next = x_0, x_n_priv = x_0, x_n_mid = x_0;

    double M1 = max(fabs(dphidx(f, left, right, left, a)), fabs(dphidx(f, left, right, right, a))) + eps;
    double q = (1 + M1) / 2;

    x_n = phi(f, left, right, x_0, a);
    x_n_next = phi(f, left, right, x_n, a);
    do
    {
        x_n_priv = x_n;
        x_n = x_n_next;
        i++;
        x_n_mid = phi(f, left, right, x_n, a);
        x_n_next = (x_n_priv * x_n_mid - x_n * x_n) / (x_n_priv - 2 * x_n + x_n_mid);

        // cout << i << " " << x_n_priv << " " << x_n_mid << " " << x_n << " " << x_n_next << endl;

    } while (fabs(x_n_next - x_n) > eps);

    ans params = {.root = x_n_next, .f = f(x_n_next, a), .start = x_0, .n_iterations = i, .q = q};
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
    fout << "\\usepackage[left=20mm, top=20mm, right=20mm, bottom=20mm, nohead, footskip=10mm]{geometry}" << endl;
    fout << "\\usepackage{indentfirst}" << endl;
    fout << "\\usepackage[margin=10pt,font=small,labelfont=bf, labelsep=endash]{caption}" << endl;
    fout << "\\usepackage{multirow}" << endl;
    fout << "\\usepackage{fourier} " << endl;
    fout << "\\usepackage{array}" << endl;
    fout << "\\usepackage{makecell}" << endl;
    fout << "\\usepackage{lscape}" << endl;
    fout << "\n"
         << endl;
    // fout << "\\renewcommand\\theadalign{bc}" << endl;
    // fout << "\\renewcommand\\theadfont{\\bfseries}" << endl;
    // fout << "\\renewcommand\\theadgape{\\Gape[4pt]}" << endl;
    // fout << "\\renewcommand\\cellgape{\\Gape[4pt]}" << endl;

    fout << "\\begin{document}" << endl;
    fout << "\\begin{landscape}" << endl;
    fout << "\\begin{center}" << endl;
    fout << "{\\Large Решение нелинейных алгебраических уравнений}" << endl;
    fout << "\\bigskip" << endl;
    fout << "\\end{center}" << endl;
    fout << "\n"
         << endl;

    fout << "Выполнил: Капелюшников Андрей Сергеевич"
         << "\n"
         << endl;
    fout << "Функция: $\\cos(e^{|\\sin(x)|}) - \\alpha x$" << endl;
    fout << "\n"
         << endl;

    fout << "Первая производная: $-a - \\cos(x)e^{|\\sin(x)|}sign(\\sin(x))\\sin(e^{|\\sin(x)|})$" << endl;
    fout << "\n"
         << endl;

    fout << "$\\varepsilon = $" << eps << endl;
    fout << "\n"
         << endl;

 //   fout << "\\pagebreak" << endl;


    int a = 1;

    for (int a = 1; a <= 10; a++)
    {

        fout << "\\begin{center}" << endl;
        fout << "{\\huge $\\alpha = " << a << "$, $[a^{" << a << "} _r, b^{" << a << "} _r] = [0, 0.5]$}" << endl;

        fout << "\n"
             << endl;

        fout << "\\begin{table}[h!]" << endl;
        fout << "\\centering" << endl;
        fout << "\\begin{tabular}{ccccccccccc}" << endl;
        fout << "\\hline" << endl;

        fout << "\\multirow{2}{*}{\\thead{Корень \\\\ $\\xi_r ^*$}} & \\multirow{2}{*}{\\thead{Невязка \\\\ $f_{" << a << "} (\\xi_r ^*)$}} & \\multirow{2}{*}{$x_0 ^{(r)}$} & \\multicolumn{5}{c}{Число итераций $N+1$} & \\multirow{2}{*}{$M ^* _{(r)}$} & \\multirow{2}{*}{$m ^* _{(r)}$} & \\multirow{2}{*}{$q ^* _{(r)}$} \\\\ \\cline{4-8}" << endl;
        fout << "                  &                   &                   &  \\thead{Половинное \\\\ деление}  &  \\thead{Метод \\\\ хорд}  &  \\thead{Простая \\\\ итерация}  &  \\thead{Метод \\\\ Эткена} &  \\thead{Метод \\\\ Ньютона} &                   &                   &                   \\\\ \\hline" << endl;
        ans chord = Chord_method(f, 0, 0.5, 0.5 - eps, a);
        fout << fixed << setprecision(6) << chord.root << " & " << chord.f << " & " << chord.start << " & & " << chord.n_iterations << " & & & & " << chord.M << " & " << chord.m << " \\\\ \\hline" << endl;

        ans newton = Newtones_method(f, 0, 0.5, 0.5 - eps, a);
        fout << fixed << setprecision(6) << newton.root << " & " << newton.f << " & " << newton.start << " & & & & & " << newton.n_iterations << " & " << newton.M << " & " << newton.m << " \\\\ \\hline" << endl;

        ans simple = SimpleIteration_method(f, 0, 0.5, 0.5 - eps, a);
        fout << fixed << setprecision(6) << simple.root << " & " << simple.f << " & " << simple.start << " & & & " << simple.n_iterations << " & & & & & " << simple.q << " \\\\ \\hline" << endl;

        ans half = half_method(f, 0, 0.5, 0.5 - eps, a);
        fout << fixed << setprecision(6) << half.root << " & " << half.f << " & " << half.start << " & " << half.n_iterations << " & & & & & & & "
             << " \\\\ \\hline" << endl;

        ans etkins = Etkins_method(f, 0, 0.5, 0.5 - eps, a);
        fout << fixed << setprecision(6) << etkins.root << " & " << etkins.f << " & " << etkins.start << " & & & & " << etkins.n_iterations << " & & & & " << etkins.q << " \\\\ \\hline" << endl;

        fout << "\\end{tabular}" << endl;
        fout << "\\end{table}" << endl;
        fout << "\\end{center}" << endl;
        fout << "\\bigskip" << endl;

        if (a % 2 == 0)
            fout << "\\pagebreak" << endl;
    }

    fout << "\\end{landscape}" << endl;
    fout << "\\end{document}" << endl;

    // для формирования отчёта необходимо запустить в компиляторе LaTex файл otchet.tex
    return 0;
}