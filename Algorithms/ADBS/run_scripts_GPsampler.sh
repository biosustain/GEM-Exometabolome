#!/bin/bash

counter=0
flag=true

while [ ${counter} -lt 8 ]; do
    if [ ${flag} = true ]; then
        if [ ${counter} -eq 0 ]; then # EColi raw
            sed -e "s/modelName/e_coli_core.mat/g" -e "s/looplessOption/0/g" -e "s/compactOption/0/g" < submit_GPsampler.sh | bsub
            flag=false
        elif [ ${counter} -eq 1 ]; then # EColi Loopless
            sed -e "s/modelName/e_coli_core.mat/g" -e "s/looplessOption/1/g" -e "s/compactOption/0/g" < submit_GPsampler.sh | bsub
            flag=false
        elif [ ${counter} -eq 2 ]; then # EColi Compact
            sed -e "s/modelName/e_coli_core.mat/g" -e "s/looplessOption/0/g" -e "s/compactOption/1/g" < submit_GPsampler.sh | bsub
            flag=false
        elif [ ${counter} -eq 3 ]; then # EColi Loopless + Compact
            sed -e "s/modelName/e_coli_core.mat/g" -e "s/looplessOption/1/g" -e "s/compactOption/1/g" < submit_GPsampler.sh | bsub
            flag=false
        elif [ ${counter} -eq 4 ]; then # ECells raw
            sed -e "s/modelName/iEC2997.mat/g" -e "s/looplessOption/0/g" -e "s/compactOption/0/g" < submit_GPsampler.sh | bsub
            flag=false
        elif [ ${counter} -eq 5 ]; then # ECells Loopless
            sed -e "s/modelName/iEC2997.mat/g" -e "s/looplessOption/1/g" -e "s/compactOption/0/g" < submit_GPsampler.sh | bsub
            flag=false
        elif [ ${counter} -eq 6 ]; then # ECells Compact
            sed -e "s/modelName/iEC2997.mat/g" -e "s/looplessOption/0/g" -e "s/compactOption/1/g" < submit_GPsampler.sh | bsub
            flag=false
        elif [ ${counter} -eq 7 ]; then # ECells Loopless + Compact
            sed -e "s/modelName/iEC2997.mat/g" -e "s/looplessOption/1/g" -e "s/compactOption/1/g" < submit_GPsampler.sh | bsub
            flag=false
        fi
    elif [ ${flag} = false ]; then
        sleep 10
	queueStat=$(bstat)
        if [ -z "${queueStat}" ]; then
            flag=true
            ((counter+=1))
        fi
    fi
done
