#!/bin/bash
#SBATCH -J problem3 # Job name
#SBATCH -o problem3_scaling_%j.stdout
#SBATCH -e problem3_scaling_%j.err
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH --mail-user=atiwari3@umassd.edu
#SBATCH --mail-type=all



# strong scaling #
run=1
threads="1 2 4 8 16 32 64 128"
while [[ run -le 3 ]]
do
	output_file="data/strong_scaling_result${run}.txt"
for thread in $threads
  do
   echo "$thread $(./rosenbrock 100000000 $thread)" >> $output_file
  done
((run = run + 1))
done


# weak scaling #
run=1
threads="1 2 4 8 16 32 64 128"
while [[ run -le 3 ]]
do
	output_file="data/weak_scaling_result${run}.txt"
for thread in $threads
    do
        ((N = 100000000 * $thread))
        echo "$thread $(./rosenbrock $N $thread)" >> $output_file
        ./rosenbrock $N $thread
    done
((run = run + 1))
done
wait 
