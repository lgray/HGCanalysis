#!/bin/bash


step=$1

if [ -z ${step} ]; then
    echo ""
    echo -e "\e[7mrunCalibrationSequence.sh\e[27m wraps up the steps for a basic calibration of HGCal simulation"
    echo "The step number must be provided as argument to run"
    echo "   1  produce the relevant particle gun samples"
    echo "   2  convert the EDM files to simple trees for calibration"
    echo "   3  perform the calibration of the e.m. scale of the subdetectors"
    echo "   4  perform the calibration of the hadronic sub-detectors"
    echo ""
    exit -1
fi


#launch production
if [ "${step}" -eq "1" ]; then
    exit -1
    energies=(10 20 40 50 75 100 250)
    pids=(22 211)
    for pid in ${pids[@]}; do
	for en in ${energies[@]}; do
            python scripts/submitLocalHGCalProduction.py -q 1nd -n 50 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION} -p ${pid} -n 200 -e ${en}";
	    if [[ "${pid}" -eq "22" ]]; then
		python scripts/submitLocalHGCalProduction.py -q 2nd -n 50 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}_EE_AIR -p ${pid} -n 200 -e ${en} -x";
		python scripts/submitLocalHGCalProduction.py -q 8nh -n 50 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}_EE_HEF_AIR -p ${pid} -n 200 -e ${en} -x -z";
	    fi
	done
    done
fi

#create trees
if [ "${step}" -eq "2" ]; then
    pids=(22 211)
    for pid in ${pids[@]}; do
	cmsRun test/runHGCSimHitsAnalyzer_cfg.py Single${pid}_${CMSSW_VERSION};
	mv /tmp/`whoami`/Single${pid}_${CMSSW_VERSION}_SimHits_0.root ./
	if [[ "${pid}" -eq "22" ]]; then
	    cmsRun test/runHGCSimHitsAnalyzer_cfg.py Single${pid}_${CMSSW_VERSION}_EE_AIR;
	    mv /tmp/`whoami`/Single${pid}_${CMSSW_VERSION}_EE_AIR_SimHits_0.root ./
	    cmsRun test/runHGCSimHitsAnalyzer_cfg.py Single${pid}_${CMSSW_VERSION}_EE_HEF_AIR;
	    mv /tmp/`whoami`/Single${pid}_${CMSSW_VERSION}_EE_HEF_AIR_SimHits_0.root ./
	fi
    done
fi

#EM calibration
if [ "${step}" -eq "3" ]; then
    
    samples=("Single22_${CMSSW_VERSION}_SimHits_0" "Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0" "Single22_${CMSSW_VERSION}_EE_HEF_AIR_SimHits_0")       
    vars=("edep_sim" "edep_rec")
    extraOpts=("" "--vetoTrackInt")
    for sample in ${samples[@]}; do 
	for var in ${vars[@]}; do
	    for opts in "${extraOpts[@]}"; do
		python test/analysis/runEMCalibration.py ${opts} -i ${sample}.root -v ${var};
		outDir=${sample}/${var}${opts};
		mkdir -p ${outDir};
		mv ${sample}/*.* ${outDir};
		python test/analysis/runEMCalibration.py -w ${outDir}/workspace.root -c ${outDir}/calib_uncalib.root;
		rm core*;
	    done
	done
    done
fi

#Pion calibration
if [ "${step}" -eq "4" ]; then
    python test/analysis/runPionCalibration.py -i Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0.root \
        --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/edep_sim--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/edep_sim--vetoTrackInt/calib_uncalib.root,HEB:Single22_${CMSSW_VERSION}_EE_HEF_AIR_SimHits_0/edep_sim--vetoTrackInt/calib_uncalib.root
fi