#!/bin/sh
function confirm {
    #echo -n "$1 Confirm[y/n]?"
    #read response
    #response=`echo $response|tr '[:upper:]' '[:lower:]'`
    #if [ "x$response" == "xn" ]; then
    #    echo "Aborted!"
    #    exit 0
    #fi
    echo $1
}

param=`echo $1|tr '[:upper:]' '[:lower:]'`

case $param in 
    'testrecon'*)
    cmd="sbatch --array=0-3 --export=Step=4,Lc=10,Curv=800 slurm-recon"
    echo $cmd
    eval $cmd
    ;;
    'recon'*)
    Step=3
    # noiseLevel=1
    for Step in 2 3; do
        for noiseLevel in 1; do
            for Lc in 5 6 8 10 20 Inf; do
                cmd="sbatch --export=Step=$Step,Lc=$Lc,noiseLevel=$noiseLevel slurm-recon"
                echo $cmd
                eval $cmd
            done
        done
    done
esac

