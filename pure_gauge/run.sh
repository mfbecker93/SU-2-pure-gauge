#!/bin/bash

# Default values
beta=5.6
grid=8.8.8.8
trajectories=20
mdsteps=20
trajL=1.0
checkpoint_interval=5
threads=8


# Parse command line arguments
while [[ $# -gt 0 ]]
do
    case $1 in
        --beta)
            beta=$2
            shift 2
            ;;
        --grid)
            grid=$2
            shift 2
            ;;
        --Trajectories)
            trajectories=$2
            shift 2
            ;;
        --mdsteps)
            mdsteps=$2
            shift 2
            ;;
        --trajL)
            trajL=$2
            shift 2
            ;;
        --checkpoint_interval)
            checkpoint_interval=$2
            shift 2
            ;;
        --threads)
            threads=$2
            shift 2
            ;;
        *)
            echo "Unknown option $1"
            exit 1
            ;;
    esac
done


# Create readable name
grid_name=${grid//./x}

name="L${grid_name}_beta${beta}_md${mdsteps}_trajL${trajL}_traj${trajectories}_ckpt${checkpoint_interval}"



mkdir -p raw_data/logs


echo "================================="
echo "Running pure gauge simulation"
echo "================================="
echo "Grid:        ${grid}"
echo "Beta:        ${beta}"
echo "MD steps:    ${mdsteps}"
echo "trajL:       ${trajL}"
echo "Trajectories:${trajectories}"
echo "Checkpoint:  ${checkpoint_interval}"
echo "Threads:     ${threads}"
echo "Log:         raw_data/logs/${name}.log"  
echo "================================="


# Run simulation
./build/pure_gauge \
    --grid ${grid} \
    --beta ${beta} \
    --Trajectories ${trajectories} \
    --mdsteps ${mdsteps} \
    --trajL ${trajL} \
    --checkpoint_interval ${checkpoint_interval} \
    --threads ${threads} \
    2>&1 | tee raw_data/logs/${name}.log 
