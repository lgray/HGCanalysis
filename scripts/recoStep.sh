#!/bin/bash

#
# PARSE CONFIGURATION PARAMETERS
#
WORKDIR="/tmp/`whoami`/"
JOBNB=0
STOREDIR=${WORKDIR}
while getopts "hj:t:o:w:" opt; do
    case "$opt" in
    h)
        echo ""
        echo "recoStep.sh [OPTIONS]"
	echo "     -j      job number (file to re-digitize, default=$JOBNB)"
	echo "     -o      output directory (local or eos, default=$STOREDIR)"
	echo "     -t      process tag"
	echo "     -w      workdir (default=$WORKDIR)"
	echo "     -h      help"
        echo ""
	exit 0
        ;;
    o)  STOREDIR=$OPTARG
	;;
    j)  JOBNB=$OPTARG
	;;
    w)  WORKDIR=$OPTARG
	;;
    t)  TAG=$OPTARG
        ;;
    esac
done

#
# CONFIGURE JOB
#
BASEJOBNAME=Events_${JOBNB}_RECO
BASEJOBNAME=${BASEJOBNAME/","/"_"}
OUTFILE=${BASEJOBNAME}.root
PYFILE=${BASEJOBNAME}_cfg.py
LOGFILE=${BASEJOBNAME}.log
FILESFORTAG=(`cmsLs /store/cmst3/group/hgcal/CMSSW/${TAG} | awk '{print $5}'`)
FILEIN=${FILESFORTAG[${JOB}]}


#cmsDriver.py reco --filein=${FILEIN} -n -1 --conditions DES23_62_V1::All \
#    --python_filename ${WORKDIR}/${PYFILE} --fileout file:${WORKDIR}/${OUTFILE} \
#    -s RAW2DIGI,RECO --datatier GEN-SIM-DIGI-RECO --eventcontent FEVTDEBUGHLT \
#    --conditions auto:upgradePLS3 --beamspot Gauss --magField 38T_PostLS1 \
#    --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCalMuon \
#    --geometry "Extended2023HGCalMuon,Extended2023HGCalMuonReco" \
#    --no_exec 

cmsDriver.py step3 --mc --conditions auto:upgradePLS3 --magField 38T_PostLS1 -n -1 \
    --python_filename ${WORKDIR}/${PYFILE} \
    -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent FEVTDEBUGHLT \
    --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCalMuon \
    --geometry Extended2023HGCalMuon,Extended2023HGCalMuonReco \
    --filein ${FILEIN}  --fileout file:${WORKDIR}/${OUTFILE} \
    --no_exec

#
# RUN cmsRun and at the end move the output to the required directory
#
#echo "cmsRun ${CFG} ${TAG} ${JOBNB} 1 ${MINBIASTAG} ${PILEUP} ${WORKDIR}"
#cmsRun ${CFG} ${TAG} ${JOBNB} 1 ${MINBIASTAG} ${PILEUP} ${WORKDIR} > ${WORKDIR}/${LOGFILE} 2>&1

#move output
#if [[ $STOREDIR =~ .*/store/cmst3.* ]]; then
#    cmsMkdir ${STOREDIR}
#    cmsStage -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
#    rm ${WORKDIR}/${OUTFILE}
#elif [[ $STOREDIR =~ /afs/.* ]]; then
#    cmsMkdir ${STOREDIR}
#    cp -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
#    rm ${WORKDIR}/${OUTFILE}
#fi

