#!/bin/bash


step=$1

if [ -z ${step} ]; then
    echo ""
    echo -e "\e[7mrunCalibrationSequence.sh\e[27m wraps up the steps for a basic calibration of HGCal simulation"
    echo "The step number must be provided as argument to run"
    echo "   sim        produce the relevant particle gun samples SIM-DIGI-RECO"
    echo "   integ      check production integrity"
    echo "   ntuple     convert the EDM files to simple trees for calibration"
    echo "   emcalib    perform the calibration of the e.m. scale of the subdetectors"
    echo "   pioncalib  perform the calibration of the hadronic sub-detectors"
    echo ""
    exit -1
fi

energies=(10 20 40 50 75 100 125 175 250 400 500)
pids=(22 211 130)
prods=(RECO-PU0 RECO-PU0-EE_HEF_AIR RECO-PU0-EE_AIR)
WHOAMI=`whoami`

#launch production
if [ "$step" == "sim" ]; then

    echo "********************************************"
    echo "launching SIM production"
    echo "********************************************"    

    for pid in ${pids[@]}; do
	for en in ${energies[@]}; do

            python scripts/submitLocalHGCalProduction.py -q 1nw -n 100 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/RECO-PU0 -p ${pid} -n 250 -e ${en} -f";
	    
	    if [ "$pid" == "130" ]; then
		continue
	    fi

	    #the following simulations will remove EE and EE+HEF sequentially
	    python scripts/submitLocalHGCalProduction.py -q 2nd -n 100 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/RECO-PU0-EE_AIR -p ${pid} -n 250 -e ${en} -x -f";	    
	    python scripts/submitLocalHGCalProduction.py -q 8nh -n 100 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/RECO-PU0-EE_HEF_AIR -p ${pid} -n 250 -e ${en} -x -z -f";

	done
    done
fi

#integrity checks
if [ "$step" == "integ" ]; then

    echo "********************************************"
    echo "checking production integrity"
    echo "********************************************"    

    for pid in ${pids[@]}; do
	for prod in ${prods[@]}; do
	    python scripts/checkProductionIntegrity.py -d /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/${prod} &
	done
    done
fi

#create trees
if [ "${step}" == "ntuple" ]; then

    echo "********************************************"
    echo "ntuplizing for analysis"
    echo "********************************************"    

    for pid in ${pids[@]}; do
	for prod in ${prods[@]}; do
	    inputFiles=(`cmsLs /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/${prod} | awk '{print $5}'`);
	    nFiles=${#inputFiles[@]};
	    nJobs=$((nFiles/50));
	    echo "******* Starting with $Single${pid}_${CMSSW_VERSION}/${prod} ${nFiles} will be analyzed in ${nJobs} jobs"
	    for i in `seq 0 ${nJobs}`; do
		startFile=$((i*=50));
		cmsRun test/runHGCSimHitsAnalyzer_cfg.py Single${pid}_${CMSSW_VERSION}/${prod} ${startFile} 50 &
	    done
	    nrunning="`ps -u $WHOAMI | grep -ir cmsRun | wc -l`"
	    while [ $nrunning -gt 0 ]; do
		echo "(...waiting for $nrunning jobs to finish...)"
		sleep 30;
		nrunning="`ps -u $WHOAMI | grep -ir cmsRun | wc -l`"
	    done
	    echo "******* All jobs finished, calling hadd"
	    hadd Single${pid}_${CMSSW_VERSION}_${prod}_SimHits.root /tmp/$WHOAMI/Single${pid}_${CMSSW_VERSION}_${prod}_SimHits_*.root;
	    rm /tmp/$WHOAMI/Single${pid}_${CMSSW_VERSION}_SimHits_*.root;
	    echo "******* Single${pid}_${CMSSW_VERSION}_${prod}_SimHits.root is ready for analysis ******"
	done
    done
fi

#e.m. calibration
if [ "${step}" == "emcalib" ]; then
    
    echo "********************************************"
    echo "e.m. calibration"
    echo "********************************************"
    vars=("edep_sim" "edep_rec")
    baseOpts="--vetoTrackInt"
    extraOpts=("" "--lambdaWeighting")
    for prod in ${prods[@]}; do
	sample=Single22_${CMSSW_VERSION}_${prod}_SimHits
	if [ ! -f ${sample}.root ]; then
	    continue
	fi
	echo "Launching calibration for ${sample}"
	for opt in "${extraOpts[@]}"; do
	    for var in ${vars[@]}; do
		echo "[${var}${opt}]"
		python test/analysis/runEMCalibration.py ${opt} ${baseOpts} -i ${sample}.root -v ${var};
		outDir=${sample}/${var}${opt};
		mkdir -p ${outDir};
		mv ${sample}/*.* ${outDir};
		python test/analysis/runEMCalibration.py -w ${outDir}/workspace.root -c ${outDir}/calib_uncalib.root;
		rm core*;
	    done
	done
    done
fi

#Pion calibration
if [ "${step}" -eq "5" ]; then

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

#pion calibration
if [ "${step}" -eq "7" ]; then

    echo "********************************************"
    echo "pion calibration with software compensation"
    echo "********************************************"

    vars=("edep_rec" "edep_clus")
    baseDir="/store/cmst3/group/hgcal/CMSSW/ntuples/"
    #ntuples=("Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_RECO-PU0_SimHits" "Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_pandoraRECO_SimHits" "Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_mqRECO_SimHits")
    #ntuples=("Single211_CMSSW_6_2_0_SLHC23_patch1_RECO-PU0_SimHits" "Single211_CMSSW_6_2_0_SLHC23_patch1_pandoraRECO_SimHits" "Single211_CMSSW_6_2_0_SLHC23_patch1_mqRECO_SimHits")
    ntuples=("Single211_CMSSW_6_2_0_SLHC23_patch1_pandoraRECO_SimHits")
    noCompOption="--noComp"
    for var in ${vars[@]}; do
	for ntuple in ${ntuples[@]}; do
	  
	    outDir=${ntuple}
	    mkdir -p ${outDir}/${var};

	    calibVar=${var}
	    extraOpts=""
	    if [ "${var}" = "edep_clus" ]; then
		calibVar="edep_rec";
	    fi

	    #python test/analysis/runPionCalibration.py --vetoTrackInt --vetoHEBLeaks -i ${baseDir}/${ntuple}.root --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/${calibVar}--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/${calibVar}--vetoTrackInt/calib_uncalib.root,HEB:Single22_CMSSW_6_2_0_SLHC20_EE_HEF_AIR_SimHits_0/${calibVar}--vetoTrackInt/calib_uncalib.root --hefhebComb Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/${calibVar}/HEFHEB_comb.root ${extraOpts} -v ${var};
	    
	    #mv ${outDir}/*.* ${outDir}/${var}

	    #extraOpts="--ehComb ${ntuple}/edep_rec/EEHEFHEB_comb.root"
	    extraOpts="--ehComb Single211_CMSSW_6_2_0_SLHC23_patch1_RECO-PU0_SimHits/edep_rec/EEHEFHEB_comb.root ${noCompOption}"
	    
	    python test/analysis/runPionCalibration.py --vetoTrackInt --vetoHEBLeaks -w ${outDir}/${var}/workspace.root --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/${calibVar}--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/${calibVar}--vetoTrackInt/calib_uncalib.root,HEB:Single22_CMSSW_6_2_0_SLHC20_EE_HEF_AIR_SimHits_0/${calibVar}--vetoTrackInt/calib_uncalib.root --hefhebComb Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/${calibVar}/HEFHEB_comb.root ${extraOpts} -v ${var};

	    python test/analysis/runPionCalibration.py --vetoTrackInt --vetoHEBLeaks -w ${outDir}/${var}/workspace.root --emCalib EE:Single22_${CMSSW_VERSION}_SimHits_0/${calibVar}--vetoTrackInt/calib_uncalib.root,HEF:Single22_${CMSSW_VERSION}_EE_AIR_SimHits_0/${calibVar}--vetoTrackInt/calib_uncalib.root,HEB:Single22_CMSSW_6_2_0_SLHC20_EE_HEF_AIR_SimHits_0/${calibVar}--vetoTrackInt/calib_uncalib.root --hefhebComb Single211_${CMSSW_VERSION}_EE_AIR_SimHits_0/${calibVar}/HEFHEB_comb.root ${extraOpts} --compWeights ${outDir}/${var}/swcompweights.root --calib ${outDir}/${var}/calib_uncalib.root -v ${var};

      	done
    done
fi


#e.m. calibration
if [ "${step}" -eq "8" ]; then

    echo "********************************************"
    echo "e.m. final calibration"
    echo "********************************************"

    vars=("edep_rec" "edep_clus")
    baseDir="/store/cmst3/group/hgcal/CMSSW/ntuples/"
    ntuples=("Single22_CMSSW_6_2_0_SLHC23_patch1_RECO-PU0_SimHits" "Single22_CMSSW_6_2_0_SLHC23_patch1_pandoraRECO_SimHits" "Single22_CMSSW_6_2_0_SLHC23_patch1_mqRECO_SimHits")
    for var in ${vars[@]}; do
	for ntuple in ${ntuples[@]}; do
	    
	    outDir=${ntuple}/${var};
	    mkdir -p ${outDir};

	    #python test/analysis/runEMCalibration.py --vetoTrackInt -i ${baseDir}/${ntuple}.root -v ${var};

	    #mv ${ntuple}/*.* ${outDir};

	    python test/analysis/runEMCalibration.py -w ${outDir}/workspace.root  -v ${var};
	    python test/analysis/runEMCalibration.py -w ${outDir}/workspace.root -c ${outDir}/calib_uncalib.root -v ${var};
	done
    done
fi
