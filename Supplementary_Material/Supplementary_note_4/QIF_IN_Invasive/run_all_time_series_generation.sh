#!/bin/bash

Ns=(200 500 1000 10000 100000)
ix_start=3
ix_end=40

for N in "${Ns[@]}"; do
    for ((ix=ix_start; ix<=ix_end; ix++)); do
        echo "Running: ./qif_in_ext $N $ix"
        ./qif_in_ext "$N" "$ix"
    done
done
