#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>

double sample_interval(const double a, const double b)
{
    double x = ((double)rand()) / ((double)RAND_MAX);
    return (b - a) * x + a;
}

double func(double *x, int dim)
{
    double f = 1.0;
    for (int i = 0; i < dim; i++)
    {
        f += x[i];
    }
    return f;
}

int isPowerOfFour(int num)
{
    if (num == 1)
    {
        return 0;
    }
    if (num > 0 && (num & (num - 1)) == 0)
    {
        return (num & 0x55555555) != 0;
    }
    return 0;
}

int main(int argc, char **argv)
{
    const int dim = 10; // dimension of integration region R
    // a 2-by-2 10-dimensional box, each side is length 2
    const double V = 1024; // volume of R
    const double xL = -1.0;
    const double xR = 1.0;
    srand(time(NULL));
    long long int N = atoll(argv[1]);
    double x[dim];
    double integral = 0;
    double sum = 0;
    double previous_integral = 0;
    for (long long int i = 1; i <= N; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            x[j] = sample_interval(xL, xR);
        }
        double f_i = func(x, dim);
        sum += f_i;
        if (isPowerOfFour(i))
        {
            integral = (V / i) * sum;
            printf("N = %lld, integral = %1.5e error estimate=%1.5e\n", i, integral, fabs(previous_integral - integral));
            previous_integral = integral;
        }
    }

    // integral = (V / N) * integral;
    // printf("MC points = %lld, integral = %1.5e error estimate=%1.5e", N, integral, fabs(1024.0 - integral));
    return 0;
}