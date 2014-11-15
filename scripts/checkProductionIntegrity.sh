#!/bin/bash

dir=$1
#dir=/store/cmst3/group/hgcal/CMSSW/MinBias_v4_CMSSW_6_2_0_SLHC18 
#dir=/store/cmst3/group/hgcal/CMSSW/MinBias_CMSSW_6_2_0_SLHC18

a=(`cmsLs ${dir} | awk '{print $5}'`)
totEvents=0
for i in ${a[@]}; do
    iEvents=`edmFileUtil -P ${i} | grep -ir bytes | awk '{print $6}'`;
    if [ -z "$iEvents" ]; then
	echo "Removing corrupted file";
	cmsRm ${i};
	continue;
    fi
    echo "$i contains $iEvents events";
    totEvents=`expr $totEvents + $iEvents`
done 
echo "Total events available $totEvents"