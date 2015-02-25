#!/bin/bash

#
# PARSE CONFIGURATION PARAMETERS
#
WORKDIR="/tmp/`whoami`/"
JOBNB=0
STOREDIR=${WORKDIR}
CFG=''
while getopts "hj:t:o:w:c:" opt; do
    case "$opt" in
    h)
        echo ""
        echo "recoStep.sh [OPTIONS]"
	echo "     -j      job number (file to re-reco, default=$JOBNB)"
	echo "     -o      output directory (local or eos, default=$STOREDIR)"
	echo "     -t      process tag"
	echo "     -w      workdir (default=$WORKDIR)"
	echo "     -h      help"
        echo ""
	exit 0
        ;;
    c)  CFG=$OPTARG
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
BASEJOBNAME=Events_${JOBNB}
BASEJOBNAME=${BASEJOBNAME/","/"_"}
OUTFILE=${BASEJOBNAME}.root
#PYFILE=${BASEJOBNAME}_cfg.py
LOGFILE=${BASEJOBNAME}.log
#FILESFORTAG=(`cmsLs /store/cmst3/group/hgcal/CMSSW/${TAG} | awk '{print $5}'`)
#FILEIN=${FILESFORTAG[${JOB}]}

echo "cmsRun ${CFG} ${TAG} ${JOBNB} 1 ${WORKDIR}"
cmsRun ${CFG} ${TAG} ${JOBNB} 1 ${WORKDIR} > /dev/null 2>&1

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

