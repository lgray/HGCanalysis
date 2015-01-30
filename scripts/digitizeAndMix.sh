#!/bin/bash

#
# PARSE CONFIGURATION PARAMETERS
#
PILEUP=0
CFG=${CMSSW_BASE}/src/UserCode/HGCanalysis/test/digitizeAndMix_cfg.py
WORKDIR="/tmp/`whoami`/"
JOBNB=1
INPUTF=""
STOREDIR=${WORKDIR}
while getopts "ho:m:p:j:c:i:" opt; do
    case "$opt" in
    h)
        echo ""
        echo "digitizeAndMix.sh [OPTIONS]"
	echo "     -p      pileup"
	echo "     -j      job number (file to re-digitize)"
	echo "     -o      output directory (local or eos)"
	echo "     -i      input file"
        echo "     -c      cfg file"
        echo "     -m      min. bias tag"
	echo "     -h      help"
        echo ""
	exit 0
        ;;
    c)  CFG=$OPTARG
	;;
    i) INPUTF=$OPTARG
        ;;
    p)  PILEUP=$OPTARG
        ;;
    o)  STOREDIR=$OPTARG
	;;
    j)  JOBNB=$OPTARG
	;;
    w)  WORKDIR=$OPTARG
	;;
    m)  MINBIASTAG=$OPTARG
        ;;
    esac
done

#
# CONFIGURE JOB
#
DIGIOUTNAME=Events_PU${PILEUP}.root
BASEJOBNAME=Events_${JOBNB}_PU${PILEUP}
BASEJOBNAME=${BASEJOBNAME/","/"_"}
OUTFILE=${BASEJOBNAME}.root
PYFILE=${BASEJOBNAME}_cfg.py
LOGFILE=${BASEJOBNAME}.log

#
# RUN cmsRun and at the end move the output to the required directory
#
echo "cmsRun ${CFG} ${INPUTF} ${MINBIASTAG} ${PILEUP} ${WORKDIR}"
cmsRun ${CFG} ${INPUTF} ${MINBIASTAG} ${PILEUP} ${WORKDIR} > ${WORKDIR}/${LOGFILE} 2>&1

#move output
echo $STOREDIR
if [[ $STOREDIR =~ .*/store/cmst3.* ]]; then
    cmsMkdir ${STOREDIR}
    cmsStage -f ${WORKDIR}/${DIGIOUTNAME} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${DIGIOUTNAME}
elif [[ $STOREDIR =~ /afs/.* ]]; then
    cmsMkdir ${STOREDIR}
    cp -f ${WORKDIR}/${DIGIOUTNAME} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${DIGIOUTNAME}
fi

