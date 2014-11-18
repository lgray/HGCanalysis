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

    echo "********************************************"
    echo "launching production"
    echo "********************************************"    


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

    echo "********************************************"
    echo "creating analysis trees"
    echo "********************************************"    

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
    
    echo "********************************************"
    echo "e.m. calibration"
    echo "********************************************"

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

    echo "********************************************"
    echo "pion calibration"
    echo "********************************************"

    vars=("edep_sim" "edep_rec")
    for var in ${vars[@]}; do

        #HEF + HEB calibration
	python test/analysis/runPionCalibration.py --vetoTrackInt --vetoHEBLeaks -i Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0.root --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root,HEB:Single22_CMSSW_6_2_0_SLHC20_EE_HEF_AIR_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root --noEE -v ${var}

	python test/analysis/runPionCalibration.py --vetoTrackInt --vetoHEBLeaks -w Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/workspace_uncalib_pion.root --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root,HEB:Single22_CMSSW_6_2_0_SLHC20_EE_HEF_AIR_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root --hefhebComb Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/HEFHEB_comb.root --noEE --calib Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/calib_uncalib.root

	mkdir -p Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/${var}
	mv Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/*.* Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/${var}

        #EE + HE(F+B) calibration
	python test/analysis/runPionCalibration.py --vetoTrackInt --vetoHEBLeaks -i Single211_${CMSSW_VERSION}_SimHits_0.root --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root,HEB:Single22_CMSSW_6_2_0_SLHC20_EE_HEF_AIR_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root --hefhebComb Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/${var}/HEFHEB_comb.root -v ${var}

	python test/analysis/runPionCalibration.py --vetoTrackInt --vetoHEBLeaks -w Single211_${CMSSW_VERSION}_SimHits_0/workspace_uncalib_pion.root --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root,HEB:Single22_CMSSW_6_2_0_SLHC20_EE_HEF_AIR_SimHits_0/${var}--vetoTrackInt/calib_uncalib.root --hefhebComb Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/${var}/HEFHEB_comb.root --calib Single211_${CMSSW_VERSION}_SimHits_0/calib_uncalib.root

	mkdir -p Single211_${CMSSW_VERSION}_SimHits_0/${var}
	mv Single211_${CMSSW_VERSION}_SimHits_0/*.* Single211_${CMSSW_VERSION}_SimHits_0/${var}
    done

    rm core.*
fi

#pion calibration (for now feedback for Lindsey mostly)
if [ "${step}" -eq "5" ]; then

    echo "********************************************"
    echo "pion calibration (cluster level)"
    echo "********************************************"

    vars=("edep_clus")
    for var in ${vars[@]}; do

        #EE + HE(F+B) calibration
	python test/analysis/runPionCalibration.py --vetoTrackInt --vetoHEBLeaks -i Single211_${CMSSW_VERSION}_lgray_SimHits_0.root --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/edep_rec--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/edep_rec--vetoTrackInt/calib_uncalib.root,HEB:Single22_CMSSW_6_2_0_SLHC20_EE_HEF_AIR_SimHits_0/edep_rec--vetoTrackInt/calib_uncalib.root --hefhebComb Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/edep_rec/HEFHEB_comb.root --ehComb Single211_${CMSSW_VERSION}_SimHits_0/edep_rec/EEHEFHEB_comb.root  -v ${var}

	python test/analysis/runPionCalibration.py --vetoTrackInt --vetoHEBLeaks -w Single211_${CMSSW_VERSION}_lgray_SimHits_0/workspace_uncalib_pion.root --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/edep_rec--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/edep_rec--vetoTrackInt/calib_uncalib.root,HEB:Single22_CMSSW_6_2_0_SLHC20_EE_HEF_AIR_SimHits_0/edep_rec--vetoTrackInt/calib_uncalib.root --hefhebComb Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/edep_rec/HEFHEB_comb.root --ehComb Single211_${CMSSW_VERSION}_SimHits_0/edep_rec/EEHEFHEB_comb.root  -v ${var}  --calib Single211_${CMSSW_VERSION}_lgray_SimHits_0/calib_uncalib.root

    done
fi
