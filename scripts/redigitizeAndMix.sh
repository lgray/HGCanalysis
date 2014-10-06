#!/bin/bash


#
# PARSE CONFIGURATION PARAMETERS
#
PULSETAU=0
CFG=${CMSSW_BASE}/src/UserCode/HGCanalysis/test/rerunDigitization_cfg.py
WORKDIR="/tmp/`whoami`/"
JOBNB=1
STOREDIR=${WORKDIR}
while getopts "ho:t:m:p:j:" opt; do
    case "$opt" in
    h)
        echo ""
        echo "redigitizeAndMix.sh [OPTIONS]"
	echo "     -p      pulse shape tau"
	echo "     -j      job number (file to re-digitize)"
	echo "     -o      output directory (local or eos)"
	echo "     -t      process tag"
        echo "     -m      min. bias tag"
	echo "     -h      help"
        echo ""
	exit 0
        ;;
    p)  PULSETAU=$OPTARG
        ;;
    o)  STOREDIR=$OPTARG
	;;
    j)  JOBNB=$OPTARG
	;;
    w)  WORKDIR=$OPTARG
	;;
    t)  TAG=$OPTARG
        ;;
    m)  MINBIASTAG=$OPTARG
        ;;
    esac
done

#
# CONFIGURE JOB
#
BASEJOBNAME=HGCEvents_${TAG}_${PULSETAU}_${JOBNB}
BASEJOBNAME=${BASEJOBNAME/","/"_"}
OUTFILE=${BASEJOBNAME}.root
PYFILE=${BASEJOBNAME}_cfg.py
LOGFILE=${BASEJOBNAME}.log

if [ "$PILEUP" = "local" ]; then
    PILEUP="${PILEUP} --pileup_input ${PILEUPINPUT}"
fi

#
# RUN cmsRun and at the end move the output to the required directory
#
echo "cmsRun ${CFG} ${PULSETAU} ${TAG} ${JOBNB} 1 ${MINBIASTAG}"
cmsRun ${CFG} ${PULSETAU} ${TAG} ${JOBNB} 1 ${MINBIASTAG} ${WORKDIR} > ${WORKDIR}/${LOGFILE} 2>&1

#move output
if [[ $STOREDIR =~ .*/store/cmst3.* ]]; then
    cmsMkdir ${STOREDIR}
    cmsStage -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${OUTFILE}
elif [[ $STOREDIR =~ /afs/.* ]]; then
    cmsMkdir ${STOREDIR}
    cp -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${OUTFILE}
fi

