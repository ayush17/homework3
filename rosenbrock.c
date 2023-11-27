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

double rosenBrock(double *x, int dim)
{

    double f = 0.0;

    for (int i = 0; i < dim; ++i)
    {
        f += (1. - x[i]) * (1. - x[i]) + 100. * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]);
    }

    f = -f;
    return f;
}

int main(int argc, char **argv)
{

    double startTime = omp_get_wtime();
    const int dim = 10;
    const double xL = -1.0;
    double V = 1;
    const double xR = 1.0;
    long long int N = atoll(argv[1]);
    const int NThreads = atoi(argv[2]);
    long long int start_loop = 1;
    long long int end_loop = 4;
    bool stopper = true;
    double x[dim];
    double sum = 0;

    omp_set_num_threads(NThreads);

    unsigned int seed = time(NULL) * (int)(omp_get_thread_num() + 1);

//  parallel
#pragma omp parallel firstprivate(x) shared(start_loop, end_loop, stopper) reduction(+ : sum)
    {
        const int thread_rank = omp_get_thread_num();

        // printf("Thread %i \n", thread_rank);
        while (stopper)
        {
#pragma omp for
            for (long long int i = start_loop; i <= end_loop; ++i)
            {
                sample_rand_r(xL, xR, dim, &seed, x);
                double f_i = exp(rosenBrock(x, dim));
                // #pragma omp atomic
                sum += f_i;
            }

#pragma omp barrier
// thread 0  will execute this block
#pragma omp master
            {
                // This is just V/N *(rosenbrookfunction value)
                double integral = V * sum / end_loop;
                printf("%lld %1.3e\n", end_loop, integral);
                start_loop = end_loop + 1;
                if (start_loop >= N)
                {
                    stopper = false;
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
        }
    }

    double endTime = omp_get_wtime();
    double walltime = endTime - startTime;
    printf("%f", walltime);
    return 0;
}
