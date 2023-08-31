#!/bin/bash

# Script Description:
# This Bash script automates the parallel execution of a computational task by dividing it into smaller chunks (faster processing)
# and submitting each chunk as a separate job to the SLURM job scheduler. The script takes three arguments:
# the starting step (start), the ending step (end), and the number of chunks (num_chunks) to create.
# It then calculates the chunk size to distribute the workload evenly and processes each chunk independently.
# Example use for 6 jobs from step 100 to 200 and put output in folder named "test": sh create_all_gammas.sh /global/scratch/users/yvesrobert/HxF_dev/Cases/HTR10_restart_P1T09/wrk_Serpent/ test 100 200 6


# Input arguments: start, end, num_chunks
start=$3
end=$4
num_chunks=$5

# Calculate the chunk size based on the number of chunks
chunk_size=$(( (end - start + 1) / num_chunks ))

# Adjust the chunk size if the division isn't even
remainder=$(( (end - start + 1) % num_chunks ))
if [ $remainder -gt 0 ]; then
    chunk_size=$((chunk_size + 1))
fi

current=$start
mkdir -p $2
cp make_gammas.py $2
cp create_gammas.sub $2
cp sss_environment $2
cp -r detector_data $2
cd $2

while [ $current -le $end ]
do
    # Calculate the end value for this chunk
    current_end=$((current + chunk_size - 1))

    # Check if the current_end exceeds the end value
    if [ $current_end -gt $end ]; then
        current_end=$end
    fi

    # Create a new directory for the current chunk
    NAME="gammas_${current}_${current_end}"
    mkdir -p ${NAME}

    # Copy necessary files into the new directory
    cp make_gammas.py ${NAME}
    cp create_gammas.sub ${NAME}
    cp sss_environment ${NAME}
    cp -r detector_data ${NAME}

    # Change the working directory to the new directory
    cd ${NAME}

    # Submit the job using sbatch with the correct job name and pass the arguments to the main script
    mkdir -p "./Logs"
    sbatch -J ${NAME} -o "./Logs/${NAME}.out" -e "./Logs/${NAME}.err"  create_gammas.sub "$1" "$current" "$current_end"

    # Go back to the previous directory for the next iteration
    cd ..

    # Update the current value for the next iteration
    current=$((current_end + 1))
done