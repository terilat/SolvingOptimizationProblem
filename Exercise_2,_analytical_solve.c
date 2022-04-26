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

int P = 5; //Порядок метода
int N = 4; //Размерность системы уравнений
double eps = pow(10, -8);  //Погрешность вычислений
//номера индексов для параметров пристрелки
int alpha_1 = 2;
int alpha_2 = 3;
//номера индексов для вектора невязок
int beta_1 = 0;
int beta_2 = 3;
//Норма для метода Рунге-Кутта
double Norma(double x[], double y[]){
    int i;
    double max = fabs(x[0] - y[0]);
    for(i = 1; i < N; ++i){
        if(max < fabs(x[i] - y[i])) max = fabs(x[i] - y[i]);
    }
    return max/(pow(2, P)-1);
}
//Обычная евклидова норма
double UsualNorma(double x[]){
    int i;
    double sum = 0;
    for(i = 0; i < N; ++i) sum += x[i]*x[i];
    return sqrt(sum);
}
//Векторное поле в точке x[]
int f(double t, double x[], double res[]){
    res[0] = x[1];
    if (x[3] > 1) res[1] = 1;
    else{
        if (x[3] < -1) res[1] = -1;
        else res[1] = x[3];
    }
    //res[1] = x[3];

    res[2] = -x[0];
    res[3] = -x[2]-x[1];
    return 1;
}
//Копирование элементов массива from в массив to
int Move(double from[], double to[]){
    int i;
    for(i = 0; i < N; ++i) to[i] = from[i];
    return 1;
}
//Вывод вектора x в консоль
void PrintVector(double x[]){
    int i;
    for(i = 0; i < N; ++i) printf("%e ", x[i]);
    printf("\n");
}
//Условие при пересечении +1
int Rule1(double x_l[], double x[], double t){
    return (x[3]-1)*(x_l[3]-1) < 0;
}
//Условие при пересечении -1
int Rule2(double x_l[], double x[], double t){
    return (x[3]+1)*(x_l[3]+1) < 0;
}
//Вычисление коэффициентов уравнения
//когда |p_y| <= 1
int EnumerateKoefAnalyticSolve1(double t, double x[], double k[]){
    double a = sqrt((sqrt(5)-1)/2), b = sqrt((sqrt(5)+1)/2);
    k[0] = (exp(-a*t)*(-x[2]+a*x[3]+a*b*b*x[0]-x[1]+b*b*x[1]))/(2*a*(a*a+b*b));
    k[1] = (exp(a*t)*(x[2]+a*x[3]+a*b*b*x[0]+x[1]-b*b*x[1]))/(2*a*(a*a+b*b));
    k[2] = -(b*(x[3]-a*a*x[0])*cos(b*t)+(x[2]+x[1]+a*a*x[1])*sin(b*t))/(b*(a*a+b*b));
    k[3] = ((b*(-x[3]+a*a*x[0])+(x[2]+x[1]+a*a*x[1])*cos(b*t)/sin(b*t))*sin(b*t))/(b*(a*a+b*b));
}
//когда p_y > 1
int EnumerateKoefAnalyticSolve2(double t, double x[], double k[]){
    k[0] = x[1] - t;
    k[1] = x[0] - t*t/2-k[0]*t;
    k[2] = -t*t*t/6-k[0]/2*t*t-k[1]*t-x[2];
    k[3] = x[3] - t*t*t*t/24-k[0]/6*t*t*t-(k[1]-1)*t*t/2-(k[2]-k[0])*t;
}
//когда p_y < -1
int EnumerateKoefAnalyticSolve3(double t, double x[], double k[]){
    k[0] = x[1] + t;
    k[1] = x[0] + t*t/2-k[0]*t;
    k[2] = t*t*t/6-k[0]/2*t*t-k[1]*t-x[2];
    k[3] = x[3] + t*t*t*t/24-k[0]/6*t*t*t-(k[1]+1)*t*t/2-(k[2]-k[0])*t;
}

//Метод Рунге-Кутта
void RungeKutta(double h, double x_0[], double t_0, double T, char str[], double tol, double y[]){
    FILE *fout, *fout_h;
    double x[N], dx[N], t = t_0; //dx[N] - массив, которых хранит набольшую разницу между численным решением в каком-то узле и аналитическим решением, t - переменная для времени
    double x_l[N]; //хранит значения x на предыдущей итерации
    int i, j, l; //итераторы для циклов
    double tmp_x[N], res[N]; //tmp_x - для вычисления промежуточных значений, res[N] - массив, для хранения компонент векторного поля в точке
    double kx[8][N]; //коэффициенты k при вычислении сумм в методе Рунге-Кутта
    double x_k[N], x_p[N]; //x_k[N] - значения с крышечкой, x_p[N] - переменные для суммирования
    int IsOpen = 0; //IsOpen - флаг открытия файла. если файл открыт, то запись в него идет промежуточных значений, если нет, то не идет)))
    double fac = 0.9, facmax = 1.5, facmin = 0.7; //параметры для автошага
    double delta = 0; //Переменная для вычисления погрешности
    int steps = 0; //число шагов
    int error = 0; //число для отслеживания зацикливания
    double t_n;
    int FlagPlus = 1, FlagMinus = 1; //флаг для метод хорд, чтобы не заходить повторно в цикл, когда вышли из него
    int count  = 0; //счетчик для метода хорд на случай закцикливания

    Move(x_0, x);
    Move(x_0, x_l);
    for(i = 0; i < N; ++i) dx[i] = 0;
    if (strcmp(str, "NAN")){
        fout = fopen(str, "a");
        IsOpen = 1;
    }

    while(t < T && steps < 10000){
        //printf("%e: ", t);
        //PrintVector(x);
        //printf("              ");
        //PrintVector(x_l);

        //чтобы последняя точка была вычислена ровно в последний момент времени
        if (t+h > T) h = T - t;

        //вычисление сумм аргументов
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
        //условие на то, допустима ли точность следующего шага
        if (Norma(x_p, x_k) < tol){
            Move(x, x_l);
            Move(x_p, x);

            /*
            if (Rule1(x_l, x, t) && FlagPlus){
                FlagPlus = 0;
                printf("+1: %e %e %e\n", x[3], x_l[3], t);
                Move(x_l, x);
                count = 0;
                while(fabs(x[3]-1) > eps && count < 1000){
                    //вычисление сумм аргументов
                    for(i = 1; i < 8; ++i){
                        Move(x, tmp_x);
                        for(j = 1; j < i; ++j){
                            for(l = 0; l < N; ++l) tmp_x[l] += h*a[i][j]*kx[j][l];
                        }
                        f(t + c[i]*h, tmp_x, res);
                        for(l = 0; l < N; ++l) kx[i][l] = res[l];
                    }
                    Move(x, x_p);
                    for(i = 1; i < 8; ++i){
                        for(l = 0; l < N; ++l) x_p[l] += h*b[i]*kx[i][l];
                    }
                    if(Rule1(x, x_p, t)){
                        //t_n = t;
                        //t = ((x_p[3]-1) * t - (x[3]-1)*(t+h))/(x_p[3] - x[3]);
                        //h = t - t_n;
                        h = h / 2;
                        //t = t_n;
                    }
                    else {
                        Move(x_p, x);
                        t += h;
                        delta += Norma(x_p, x_k);
                    }
                    count += 1;
                }
                if (count > 999) {
                    printf("U SUKA + | %e\n", t);
                    PrintVector(x);
                }
                count = 0;
                //Move(x, x_l);
                printf("%e %e %e\n", x[beta_1], x[beta_2], t);
            }
            else FlagPlus = 1;
            if (Rule2(x_l, x, t) && FlagMinus){
                FlagMinus = 0;
                printf("-1: %e %e %e\n", x[3], x_l[3], t);
                Move(x_l, x);
                count = 0;
                while(fabs(x[3]+1) > eps && count < 1000){
                    //вычисление сумм аргументов
                    for(i = 1; i < 8; ++i){
                        Move(x, tmp_x);
                        for(j = 1; j < i; ++j){
                            for(l = 0; l < N; ++l) tmp_x[l] += h*a[i][j]*kx[j][l];
                        }
                        f(t + c[i]*h, tmp_x, res);
                        for(l = 0; l < N; ++l) kx[i][l] = res[l];
                    }
                    Move(x, x_p);
                    for(i = 1; i < 8; ++i){
                        for(l = 0; l < N; ++l) x_p[l] += h*b[i]*kx[i][l];
                    }
                    if(Rule2(x, x_p, t)){
                        //t_n = t;
                        //t = ((x_p[3]+1) * t - (x[3]+1)*(t+h))/(x_p[3] - x[3]);
                        //h = t - t_n;
                        h = h / 2;
                        //t = t_n;
                    }
                    else {
                        Move(x_p, x);
                        t += h;
                        delta += Norma(x_p, x_k);
                    }
                    count += 1;
                    //printf("In -1: %e %e   ", t, h);
                    //PrintVector(x);
                }
                if (count > 999) {
                    printf("U SUKA - | %e\n", t);
                    PrintVector(x);
                }
                count = 0;
                //Move(x, x_l);
                printf("%e %e %e\n", x[beta_1], x[beta_2], t);

            }
            else FlagMinus = 1;

            if (FlagMinus && FlagMinus) */ t += h;
            steps += 1;

            if (IsOpen) fprintf(fout, "%e %e %e %e %e\n", t, x[0], x[1], x[2], x[3]);
            delta += Norma(x_p, x_k);
        }
        //автошаг
        h = h * Min(facmax, Max(facmin, fac * pow(tol/Norma(x, x_k), 1.0/(P + 1))));
        if (h > 0.005) h = 0.005;
        error += 1;
        //printf("%d | ", error);
        //PrintVector(x);
        //printf("%e %e %e %e\n_________________________\n", Norma(x_p, x_k), tol, h, t);
    }
    //fprintf(fout_h, "__________________________\n");

    printf("RK: ");
    PrintVector(x);
    printf("%e\n", t);

    for(i = 0; i < N; ++i) y[i] = x[i];

    //fprintf(fout, "%s & %d & %s & %s & %s & %s & %s \\\\\n", strtol, steps, strx0, strx1, strdx0, strdx1, strdelta);*/
    //printf("%s & %d & %s & %s & %s & %s & %s \\\\\n", strtol, steps, strx0, strx1, strdx0, strdx1, strdelta);
    //printf("RK: %d %f ||", steps, delta);
    //printf("RK: ");
    //PrintVector(x);
    if (IsOpen) fclose(fout);
}

int AnalyticalSolve1(double x[], double h, double *t, double T){
    FILE *fin;
    int i, j;
    double y[N], k[N];

    EnumerateKoefAnalyticSolve3(*t, x, k);
    PrintVector(k);
    printf("____________________________\n");
    printf("%e\n", *t);
    fin = fopen("AnalyticSolve.txt", "w");
    while(*t <= T && x[3] < -1){
        x[0] = -*t*(*t)/2+k[0]*(*t)+k[1];
        x[1] = -*t+k[0];
        x[2] = *t**t**t/6-k[0]/2**t**t-k[1]**t-k[2];
        x[3] = -*t**t**t**t/24+k[0]/6**t**t**t+(k[1]+1)**t**t/2+(k[2]-k[0])**t+k[3];
        fprintf(fin, "%e %e %e %e %e\n", *t, x[0], x[1], x[2], x[3]);
        printf("%e %e\n", *t, x[3]);
        *t += h;
        if (*t > T) *t = T;
        if (*t == T) *t += h;
    }
    printf("FINISH\n");
    fclose(fin);
}
int AnalyticalSolve2(double x[], double h, double *t, double T){
    FILE *fin;
    int i, j;
    double a = sqrt((sqrt(5)-1)/2), b = sqrt((sqrt(5)+1)/2);
    double y[N], k[N];

    PrintVector(x);
    EnumerateKoefAnalyticSolve1(*t, x, k);
    PrintVector(k);
    printf("____________________________\n");
    printf("%e\n", *t);
    fin = fopen("AnalyticSolve.txt", "a");
    while(*t <= T && x[3] < 1){
        x[0] = k[0]*exp(a**t) + k[1]*exp(-a**t)+k[2]*cos(b**t)+k[3]*sin(b**t);
        x[1] = k[0]*a*exp(a**t) - k[1]*a*exp(-a**t) - k[2]*b*sin(b**t) + k[3]*b*cos(b**t);
        x[2] = -k[0]*a*(a*a+1)*exp(a**t) + k[1]*a*(a*a+1)*exp(-a**t) + k[2]*b*(1-b*b)*sin(b**t) + k[3]*b*(b*b-1)*cos(b**t);
        x[3] = k[0]*a*a*exp(a**t) + k[1]*a*a*exp(-a**t) - k[2]*b*b*cos(b**t) - k[3]*b*b*sin(b**t);
        fprintf(fin, "%e %e %e %e %e\n", *t, x[0], x[1], x[2], x[3]);
        printf("%e %e\n", *t, x[3]);
        *t += h;
        if (*t > T) *t = T;
        if (*t == T) *t += h;
    }
    printf("FINISH\n");
    fclose(fin);
}
int AnalyticalSolve3(double x[], double h, double *t, double T){
    FILE *fin;
    int i, j;
    double y[N], k[N];
    printf("____________________________\n");
    PrintVector(x);
    EnumerateKoefAnalyticSolve2(*t, x, k);
    PrintVector(k);
    printf("____________________________\n");
    printf("%e\n", *t);
    fin = fopen("AnalyticSolve.txt", "a");
    while(*t <= T && x[3] > 1){
        x[0] = *t**t/2+k[0]**t+k[1];
        x[1] = *t+k[0];
        x[2] = -*t**t**t/6-k[0]/2**t**t-k[1]**t-k[2];
        x[3] = *t**t**t**t/24+k[0]/6**t**t**t+(k[1]-1)**t**t/2+(k[2]-k[0])**t+k[3];
        fprintf(fin, "%e %e %e %e %e\n", *t, x[0], x[1], x[2], x[3]);
        printf("%e %e\n", *t, x[3]);
        *t += h;
        if (*t > T) *t = T;
        if (*t == T) *t += h;
    }
    printf("FINISH\n");
    fclose(fin);
}
int AnalyticalSolve4(double x[], double h, double *t, double T){
    FILE *fin;
    int i, j;
    double a = sqrt((sqrt(5)-1)/2), b = sqrt((sqrt(5)+1)/2);
    double y[N], k[N];

    EnumerateKoefAnalyticSolve1(*t, x, k);
    PrintVector(k);
    printf("____________________________\n");
    printf("%e\n", *t);
    fin = fopen("AnalyticSolve.txt", "a");
    while(*t <= T && x[3] < 1){
        x[0] = k[0]*exp(a**t) + k[1]*exp(-a**t)+k[2]*cos(b**t)+k[3]*sin(b**t);
        x[1] = k[0]*a*exp(a**t) - k[1]*a*exp(-a**t) - k[2]*b*sin(b**t) + k[3]*b*cos(b**t);
        x[2] = -k[0]*a*(a*a+1)*exp(a**t) + k[1]*a*(a*a+1)*exp(-a**t) + k[2]*b*(1-b*b)*sin(b**t) + k[3]*b*(b*b-1)*cos(b**t);
        x[3] = k[0]*a*a*exp(a**t) + k[1]*a*a*exp(-a**t) - k[2]*b*b*cos(b**t) - k[3]*b*b*sin(b**t);
        fprintf(fin, "%e %e %e %e %e\n", *t, x[0], x[1], x[2], x[3]);
        printf("%e %e\n", *t, x[3]);
        *t += h;
        if (*t > T) *t = T;
        if (*t == T) *t += h;
    }
    printf("FINISH\n");
    fclose(fin);
}
int main(){
    double x[N];
    double y[N];
    double t = 0;
    remove("CheckParametrs.txt");

    x[0] = 0;
    x[1] = 0;
    x[2] = -233.4462;
    x[3] = -1335.591;
    RungeKutta(0.01, x, 0, 20, "CheckParametrs.txt", pow(10, -12), y);
    x[0] = 0;
    x[1] = 0;
    x[2] = -233.4462;
    x[3] = -1335.591;
    AnalyticalSolve1(x, 0.001, &t, 20);
    printf("%e\n", t);
    AnalyticalSolve2(x, 0.001, &t, 20);
    AnalyticalSolve3(x, 0.001, &t, 20);
    AnalyticalSolve3(x, 0.001, &t, 20);
    return 0;
}
