#include <stdio.h>
#include <math.h>
#include <string.h>
#include "OtherFunction.c"
//Метод Рунге-Кутты 5-го порядка для двумерной вектор функции f(t, x, y) с динамическим изменением шага

double a[8][8] = {{0, 0, 0, 0, 0, 0, 0, 0},
                  {0, 0, 0, 0, 0, 0, 0, 0},
                  {0, 1.0/5, 0, 0, 0, 0, 0, 0},
                  {0, 3.0/40, 9.0/40, 0, 0, 0, 0, 0},
                  {0, 44.0/45, -56.0/15, 32.0/9, 0, 0, 0, 0},
                  {0, 19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729, 0, 0, 0},
                  {0, 9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656, 0, 0},
                  {0, 35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0}};
double b[8] = {0, 35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0};
double b_k[8] = {0, 5179.0/57600, 0, 7571.0/16695, 393.0/640, -92097.0/339200, 187.0/2100, 1.0/40};
double c[8] = {0, 0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1.0, 1.0};

int P = 8; //Порядок метода
int N = 2;

double Norma(double x[], double y[]){
    int i;
    double max = fabs(x[0] - y[0]);
    for(i = 1; i < N; ++i){
        if(max < fabs(x[i] - y[i])) max = fabs(x[i] - y[i]);
    }
    return max/(pow(2, P)-1);
}

int f(double t, double x[], double res[]){
    res[0] = x[1];
    res[1] = -x[0];
    return 1;
}

int Move(double from[], double to[]){
    int i;
    for(i = 0; i < N; ++i) to[i] = from[i];
    return 1;
}

void RungeKutta(double h, double x_0[], double t_0, double T, char str[], double tol, double y[]){
    FILE *fout, *fout_h;
    double x[N], dx[N], t = t_0, t_p;
    int i, j, l;
    double tmp_x[N], res[N];
    double kx[8][N];
    double x_k[N], x_p[N]; //x_k, y_k - значения с крышечкой, x_p, y_p - переменные для суммирования
    int flag = 1;
    double fac = 0.9, facmax = 1.5, facmin = 0.7;
    double h_0 = h;
    double delta = 0;
    int steps = 0;

    char strtol[20], strx0[20], strx1[20], strdx0[20], strdx1[20], strdelta[20];

    Move(x_0, x);
    for(i = 0; i < N; ++i) dx[i] = 0;
    fout = fopen(str, "a");
    while(t < T){
        if (t+h > T) h = T - t;
        t_p = t;

        for(i = 1; i < 8; ++i){
            Move(x, tmp_x);
            for(j = 1; j < i; ++j){
                for(l = 0; l < N; ++l) tmp_x[l] += h*a[i][j]*kx[j][l];
            }
            f(t + c[i]*h, tmp_x, res);
            for(l = 0; l < N; ++l) kx[i][l] = res[l];
        }

        Move(x, x_k);
        Move(x, x_p);
        for(i = 1; i < 8; ++i){
            for(l = 0; l < N; ++l){
                x_p[l] += h*b[i]*kx[i][l];
                x_k[l] += h*b_k[i]*kx[i][l];
            }
        }
        if (Norma(x_p, x_k) < tol){
            Move(x_p, x);
            t += h;

            steps += 1;

            if(dx[0] < fabs(x[0] - sin(t))) dx[0] = fabs(x[0] - sin(t));
            if(dx[1] < fabs(x[1] - cos(t))) dx[1] = fabs(x[1] - cos(t));

            delta += Norma(x_p, x_k);
        }
        h = h * Min(facmax, Max(facmin, fac * pow(tol/Norma(x_p, x_k), 1.0/(P + 1))));
    }
    //fprintf(fout_h, "__________________________\n");


    for(i = 0; i < N; ++i) y[i] = x[i];

    ConvertDoubleToLatex(tol, 1, strtol);
    ConvertDoubleToLatex(fabs(x[0]), 1, strx0);
    ConvertDoubleToLatex(fabs(x[1]-cos(T)), 1, strx1);
    ConvertDoubleToLatex(dx[0], 1, strdx0);
    ConvertDoubleToLatex(dx[1], 1, strdx1);
    ConvertDoubleToLatex(delta, 1, strdelta);

    fprintf(fout, "%s & %d & %s & %s & %s & %s & %s \\\\\n", strtol, steps, strx0, strx1, strdx0, strdx1, strdelta);
    printf("%s & %d & %s & %s & %s & %s & %s \\\\\n", strtol, steps, strx0, strx1, strdx0, strdx1, strdelta);


    fclose(fout);
}

double NumberRunge(double x1, double x2, double x3){
    return fabs((x1-x2)/(x2-x3));
}

int main(){
    double x[N];
    double r1[N], r2[N], r3[N];
    FILE *fout;
    x[0] = 0;
    x[1] = 1;

    remove("ResultOscillator.txt");
    fout = fopen("ResultOscillator.txt", "a");

    RungeKutta(0.1, x, 0, M_PI, "ResultOscillator.txt", pow(10, -8), r1);
    RungeKutta(0.1, x, 0, M_PI, "ResultOscillator.txt", pow(10, -10), r2);
    RungeKutta(0.1, x, 0, M_PI, "ResultOscillator.txt", pow(10, -12), r3);
    printf("Pi: %f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    fprintf(fout, "%f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    RungeKutta(0.1, x, 0, 10*M_PI, "ResultOscillator.txt", pow(10, -8), r1);
    RungeKutta(0.1, x, 0, 10*M_PI, "ResultOscillator.txt", pow(10, -10), r2);
    RungeKutta(0.1, x, 0, 10*M_PI, "ResultOscillator.txt", pow(10, -12), r3);
    printf("10Pi: %f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    fprintf(fout, "%f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    RungeKutta(0.1, x, 0, 100*M_PI, "ResultOscillator.txt", pow(10, -8), r1);
    RungeKutta(0.1, x, 0, 100*M_PI, "ResultOscillator.txt", pow(10, -10), r2);
    RungeKutta(0.1, x, 0, 100*M_PI, "ResultOscillator.txt", pow(10, -12), r3);
    printf("100Pi: %f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    fprintf(fout, "%f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    RungeKutta(0.1, x, 0, 1000*M_PI, "ResultOscillator.txt", pow(10, -8), r1);
    RungeKutta(0.1, x, 0, 1000*M_PI, "ResultOscillator.txt", pow(10, -10), r2);
    RungeKutta(0.1, x, 0, 1000*M_PI, "ResultOscillator.txt", pow(10, -12), r3);
    printf("1000Pi: %f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    fprintf(fout, "%f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    RungeKutta(0.1, x, 0, 10000*M_PI, "ResultOscillator.txt", pow(10, -8), r1);
    RungeKutta(0.1, x, 0, 10000*M_PI, "ResultOscillator.txt", pow(10, -10), r2);
    RungeKutta(0.1, x, 0, 10000*M_PI, "ResultOscillator.txt", pow(10, -12), r3);
    printf("10000Pi: %f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    fprintf(fout, "%f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    RungeKutta(0.1, x, 0, 100000*M_PI, "ResultOscillator.txt", pow(10, -8), r1);
    RungeKutta(0.1, x, 0, 100000*M_PI, "ResultOscillator.txt", pow(10, -10), r2);
    RungeKutta(0.1, x, 0, 100000*M_PI, "ResultOscillator.txt", pow(10, -12), r3);
    printf("100000Pi: %f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    fprintf(fout, "%f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    RungeKutta(0.1, x, 0, 1000000*M_PI, "ResultOscillator.txt", pow(10, -8), r1);
    RungeKutta(0.1, x, 0, 1000000*M_PI, "ResultOscillator.txt", pow(10, -10), r2);
    RungeKutta(0.1, x, 0, 1000000*M_PI, "ResultOscillator.txt", pow(10, -12), r3);
    printf("1000000Pi: %f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    fprintf(fout, "%f %f\n", NumberRunge(r1[0], r2[0], r3[0]), NumberRunge(r1[1], r2[1], r3[1]));
    return 0;
}
