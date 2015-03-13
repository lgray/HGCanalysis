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
    echo "   picalib    perform the calibration of the hadronic sub-detectors"
    echo "   wrapup     wraps up the results for PF"
    echo ""
    exit -1
fi

energies=(2 3 5 8 10 20 40 50 75 100 125 175 250 400 500)
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
    pids=(211 22)
    #prods=(RECO-PU0)
    fileSplit=50
    for pid in ${pids[@]}; do
	for prod in ${prods[@]}; do
	    inputFiles=(`cmsLs /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/${prod} | awk '{print $5}'`);
	    nFiles=${#inputFiles[@]};
	    nJobs=$((nFiles/${fileSplit}));
	    echo "******* Starting with $Single${pid}_${CMSSW_VERSION}/${prod} ${nFiles} will be analyzed in ${nJobs} jobs"
	    for i in `seq 0 ${nJobs}`; do
		startFile=$((i*=${fileSplit}));
		cmsRun test/runHGCSimHitsAnalyzer_cfg.py Single${pid}_${CMSSW_VERSION}/${prod} ${startFile} ${fileSplit} &
	    done
	    nrunning="`ps -u $WHOAMI | grep -ir cmsRun | wc -l`"
	    while [ $nrunning -gt 0 ]; do
		echo "(...waiting for $nrunning jobs to finish...)"
		sleep 30;
		nrunning="`ps -u $WHOAMI | grep -ir cmsRun | wc -l`"
	    done
	    echo "******* All jobs finished, calling hadd"
	    hadd /tmp/$WHOAMI/Single${pid}_${CMSSW_VERSION}_${prod}_SimHits.root /tmp/$WHOAMI/Single${pid}_${CMSSW_VERSION}_${prod}_SimHits_*.root > /tmp/$WHOAMI/Single${pid}_${CMSSW_VERSION}_${prod}_SimHits.log; 
	    #rm /tmp/$WHOAMI/Single${pid}_${CMSSW_VERSION}_SimHits_*.root;
	    echo "******* Single${pid}_${CMSSW_VERSION}_${prod}_SimHits.root is ready for analysis ******"
	done
    done
fi

#e.m. calibration
if [ "${step}" == "emcalib" ]; then
    prods=(RECO-PU0-EE_HEF_AIR RECO-PU0-EE_AIR RECO-PU0)
    echo "********************************************"
    echo "e.m. calibration"
    echo "********************************************"
    vars=("edep_rec" "edep_sim")
    extraOpts=("" "--lambdaWeighting")
    for prod in ${prods[@]}; do
	
	sample=Single22_${CMSSW_VERSION}_${prod}_SimHits

	if [ ! -f ${sample}.root ]; then
	    continue
	fi
	
	baseOpts="--vetoTrackInt"
	
	echo "Launching calibration for ${sample}"
	for opt in "${extraOpts[@]}"; do
	    for var in ${vars[@]}; do
		echo "[${var}${opt}]"
		outDir=${sample}/${var}${opt};		

		if [[ ${prod} =~ .*EE_HEF_AIR.* && ${var} =~ .*rec.* ]]; then
		    python test/analysis/runEMCalibration.py ${opt} ${baseOpts} -i ${sample}.root -v ${var} --useSaturatedCalibModel;
		    #python test/analysis/runEMCalibration.py -w ${outDir}/workspace.root ${baseOpts} --useSaturatedCalibModel;
		else
		    python test/analysis/runEMCalibration.py ${opt} ${baseOpts} -i ${sample}.root -v ${var};
		    #python test/analysis/runEMCalibration.py -w ${outDir}/workspace.root ${baseOpts};
		fi

		mkdir -p ${outDir};
		mv ${sample}/*.* ${outDir};

		python test/analysis/runEMCalibration.py -w ${outDir}/workspace.root -c ${outDir}/calib_uncalib.root ${baseOpts};
	    done
	done
    done
fi

#Pion calibration
if [ "${step}" == "picalib" ]; then    
    #prods=(RECO-PU0-EE_HEF_AIR RECO-PU0-EE_AIR RECO-PU0)
    #prods=(RECO-PU0-EE_HEF_AIR)
    #prods=(RECO-PU0-EE_AIR)
    prods=(RECO-PU0)  
    #prods=(RECO-PU0-EE_AIR RECO-PU0) 
    echo "********************************************"
    echo "pion calibration"
    echo "********************************************"
    #vars=("edep_sim")
    vars=("edep_rec")
    #vars=("edep_rec" "edep_sim")
    for prod in ${prods[@]}; do
	
	sample=Single211_${CMSSW_VERSION}_${prod}_SimHits

	if [ ! -f ${sample}.root ]; then
	    continue
	fi

	echo "Launching calibration for ${sample}"
	for var in ${vars[@]}; do

	    emEE=Single22_${CMSSW_VERSION}_RECO-PU0_SimHits/${var}--lambdaWeighting/calib_uncalib.root
	    emHEF=Single22_${CMSSW_VERSION}_RECO-PU0-EE_AIR_SimHits/${var}--lambdaWeighting/calib_uncalib.root
	    emHEB=Single22_${CMSSW_VERSION}_RECO-PU0-EE_HEF_AIR_SimHits/${var}--lambdaWeighting/calib_uncalib.root

	    echo "[${var}]"
	    outDir=${sample}/${var};	

	    if [[ ${prod} =~ .*EE_HEF_AIR.* ]]; then
		baseOpts="--vetoTrackInt --noEE --noHEF";
	    elif [[ ${prod} =~ .*EE_AIR.* ]]; then
		baseOpts="--vetoTrackInt --noEE --hebResp=Single211_${CMSSW_VERSION}_RECO-PU0-EE_HEF_AIR_SimHits/${var}/HEB_response.root";
	    else
		baseOpts="--vetoTrackInt"
	    fi
	    
	    #python test/analysis/runPionCalibration.py -i ${sample}.root --emCalib EE:${emEE},HEF:${emHEF},HEB:${emHEB} -v ${var} ${baseOpts}
	    #mkdir -p ${outDir}
	    #mv ${sample}/*.* ${outDir};

	    if [[ ! ${prod} =~ .*AIR.* ]]; then

		#use separate combination
		#baseOpts="--hebResp=Single211_${CMSSW_VERSION}_RECO-PU0-EE_HEF_AIR_SimHits/${var}/HEB_response.root --hefResp Single211_${CMSSW_VERSION}_RECO-PU0-EE_AIR_SimHits/${var}/HEFHEB_response_corry.root"
		#python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts} --noResCalib
		#baseOpts="${baseOpts} --eeResp ${outDir}/EEHEFHEB_response_corry.root"
		#python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts} --noResCalib
		#baseOpts="${baseOpts} --banana ${outDir}/EEHEFHEB_response_corrx_corry.root"
		#python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts}
		#mkdir -p ${sample}/${var}-IndComb
		#mv ${outDir}/*.png ${sample}/${var}-IndComb
		#mv ${outDir}/calib*.root ${sample}/${var}-IndComb
		#mv ${outDir}/*response*.root ${sample}/${var}-IndComb
		
		#trivial combination
		baseOpts="--noComp"
		python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts}
		python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts} --calib ${outDir}/calib_uncalib.root;
		mkdir -p ${sample}/${var}-TrivialComb
		mv ${outDir}/*.png ${sample}/${var}-TrivialComb
		mv ${outDir}/calib*.root ${sample}/${var}-TrivialComb
		mv ${outDir}/*response*.root ${sample}/${var}-TrivialComb

		#combine Si detectors
		#baseOpts="--hebResp=Single211_${CMSSW_VERSION}_RECO-PU0-EE_HEF_AIR_SimHits/${var}/HEB_response.root --combineSi"
		#python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts}
		#baseOpts="--hebResp=Single211_${CMSSW_VERSION}_RECO-PU0-EE_HEF_AIR_SimHits/${var}/HEB_response.root  --eeResp ${outDir}/EEHEFHEB_response_corry.root --combineSi"
		#python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts} --calib ${outDir}/calib_uncalib.root;
	    else
		if [[ ${prod} =~ .*EE_HEF_AIR.* ]]; then
		    rm -v Single211_${CMSSW_VERSION}_RECO-PU0-EE_HEF_AIR_SimHits/${var}/HEB_response.root;
		else
		    rm -v Single211_${CMSSW_VERSION}_RECO-PU0-EE_AIR_SimHits/${var}/HEFHEB_response_corry.root;
		fi
		python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts} --noResCalib
		baseOpts="${baseOpts} --hebResp=Single211_${CMSSW_VERSION}_RECO-PU0-EE_HEF_AIR_SimHits/${var}/HEB_response.root --hefResp Single211_${CMSSW_VERSION}_RECO-PU0-EE_AIR_SimHits/${var}/HEFHEB_response_corry.root"
		python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts} --noResCalib

		if [[ ${prod} =~ .*EE_AIR.* ]]; then
		    baseOpts="${baseOpts} --banana Single211_${CMSSW_VERSION}_RECO-PU0-EE_AIR_SimHits/${var}/HEFHEB_response_corrx_corry.root"
		fi
		python test/analysis/runPionCalibration.py -w ${outDir}/workspace.root -v ${var} ${baseOpts}		
	    fi

	done
    done
fi


if [ "${step}" == "wrapup" ]; then

    echo "********************************************"
    echo "Wrapping up all information for PF          "
    echo "********************************************"
    python test/analysis/wrapUpCalibrationForPF.py;
fi


##
## the following is for test purpose only
##

#pion calibration
if [ "${step}" == "testpicalib" ]; then

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
if [ "${step}" == "testemcalib" ]; then

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
