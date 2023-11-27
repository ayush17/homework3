#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>

void sample_rand_r(const double a, const double b, const int dim, unsigned int *seed, double *x)
{

    for (int i = 0; i < dim; ++i)
    {
        double tmp = ((double)rand_r(seed)) / ((double)RAND_MAX);
        tmp = (b - a) * tmp + a;
        x[i] = tmp;
    }
}

// for code testing (problems 1 and 2)
double func(double *x, int dim)
{
    double f = 1.0;

    for (int i = 0; i < dim; ++i)
    {
        f += x[i];
    }

    return f;
}

int main(int argc, char **argv)
{

    double startTime = omp_get_wtime();
    const int dim = 10;
    const double V = 1024;
    const double xL = -1.0;
    const double xR = 1.0;
    long long int N = atoll(argv[1]);
    const int NThreads = atoi(argv[2]);
    long long int start_loop = 1;
    long long int end_loop = 4;
    bool loopStoopperflag = true;
    double x[dim];
    double sum = 0;
    omp_set_num_threads(NThreads);
#pragma omp parallel firstprivate(x) shared(start_loop, end_loop, sum, loopStoopperflag)
    {
        double previous_integral = -100;
        unsigned int seed = time(NULL) * (int)(omp_get_thread_num() + 1);

        while (loopStoopperflag)
        {

#pragma omp for
            for (long long int i = start_loop; i <= end_loop; ++i)
            {
                // sample 10-dimensional box
                sample_rand_r(xL, xR, dim, &seed, x);
                double f_i = func(x, dim); // for code testing (problems 1 and 2)
#pragma omp atomic
                sum += f_i;
            }

#pragma omp barrier
#pragma omp master
            {
                double integral = V * sum / end_loop;
                double err_est = fabs(previous_integral - integral);
                printf("N= %lld Integral=%1.3e actual_estimate=%1.3e,estimated_error= %1.3e\n", end_loop, integral, fabs(1024.0 - integral), err_est);
                previous_integral = integral;
                start_loop = end_loop + 1;
                if (start_loop >= N)
                {
                    loopStoopperflag = false;
                }
                else
                {
                    if (end_loop > (int)(N / 4))
                    {
                        end_loop = N;
                    }
                    else
                    {
                        end_loop = 4 * end_loop;
                    }
                }
            }
#pragma omp barrier
        }
    }
    double endTime = omp_get_wtime();
    double walltime = endTime - startTime;
    printf("omp time is = %f\n", walltime);
    return 0;
}
