double fabs(double x){
    if (x > 0) return x;
    return -x;
}
double Max(double a, double b){
    if(a > b) return a;
    return b;
}
double Min(double a, double b){
    if(a > b) return b;
    return a;
}

/* reverse:  переворачиваем строку s на месте */
void reverse(char s[]){
    int i, j;
    char c;

    for(i = 0, j = strlen(s)-1; i<j; i++, j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}
void IntToStr(int n, char s[]){
     int i, sign;

     if ((sign = n) < 0)  /* записываем знак */
         n = -n;          /* делаем n положительным числом */
     i = 0;
     do {       /* генерируем цифры в обратном порядке */
         s[i++] = n % 10 + '0';   /* берем следующую цифру */
     } while ((n /= 10) > 0);     /* удаляем */
     if (sign < 0)
         s[i++] = '-';
     s[i] = '\0';
     reverse(s);
 }
int ConvertDoubleToLatex(double x, int accuracy, char str[]){
    double mant;
    int p, sign = 1;
    int q;
    int i, j;
    char s[20], ss[20], tmp;
    if (x < 0) str[0] = '-';
    if (x > 1){
        p = log10(x);
        q = x * pow(10, accuracy - p);
        IntToStr(q, s);
        i = 0;
        str[i++] = s[0];
        if(accuracy != 0) {
            str[i++] = '.';
            for(; i <= accuracy+1; ++i) str[i] = s[i-1];
        }
        str[i++] = '\\';
        str[i++] = 'c';
        str[i++] = 'd';
        str[i++] = 'o';
        str[i++] = 't';
        str[i++] = ' ';
        str[i++] = '1';
        str[i++] = '0';
        str[i++] = '^';
        str[i++] = '{';
        IntToStr(p, ss);
        j = i;
        for(; i < strlen(ss)+j; ++i) str[i] = ss[i-j];
        str[i++] = '}';
        str[i++] = '\0';
    }
    else{
        p = log10(1/x);
        q = x * pow(10, accuracy + 1 + p);
        IntToStr(q, s);
        i = 0;
        str[i++] = s[0];
        if(accuracy != 0) {
            str[i++] = '.';
            for(; i <= accuracy+1; ++i) str[i] = s[i-1];
        }
        str[i++] = '\\';
        str[i++] = 'c';
        str[i++] = 'd';
        str[i++] = 'o';
        str[i++] = 't';
        str[i++] = ' ';
        str[i++] = '1';
        str[i++] = '0';
        str[i++] = '^';
        str[i++] = '{';
        str[i++] = '-';
        IntToStr(p+1, ss);
        j = i;
        for(; i < strlen(ss)+j; ++i) str[i] = ss[i-j];
        str[i++] = '}';
        str[i++] = '\0';
    }
    return 1;
}
