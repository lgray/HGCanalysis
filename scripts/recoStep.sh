#!/bin/bash

#
# PARSE CONFIGURATION PARAMETERS
#
WORKDIR="/tmp/`whoami`/"
JOBNB=0
STOREDIR=${WORKDIR}
PROCESSTHESE=1
CFG=''
while getopts "hj:t:o:w:c:f:n:" opt; do
    case "$opt" in
    h)
        echo ""
        echo "recoStep.sh [OPTIONS]"
	echo "     -j      job number will be used to skip (# events to process x job number) "
	echo "     -o      output directory (local or eos, default=$STOREDIR)"
	echo "     -n      number of events to process (1 by default)"
	echo "     -f      file to process in directory"
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
    f)  FILENUMBER=$OPTARG
	;;
    n)  PROCESSTHESE=$OPTARG
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
SKIPTHESE=$(((JOBNB-1)*PROCESSTHESE))
BASEJOBNAME=Events_${FILENUMBER}_${SKIPTHESE}_${PROCESSTHESE}
OUTFILE=${BASEJOBNAME}.root

echo "cmsRun ${CFG} ${TAG} ${FILENUMBER} 1 ${WORKDIR} ${SKIPTHESE} ${PROCESSTHESE}"
echo "Expecting to produce ${OUTFILE}"

cmsRun ${CFG} ${TAG} ${FILENUMBER} 1 ${WORKDIR} ${SKIPTHESE} ${PROCESSTHESE}

#move output
echo "Moving result to $STOREDIR"
if [[ $STOREDIR =~ .*/store/cmst3.* ]]; then
    cmsMkdir ${STOREDIR}
    cmsStage -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${OUTFILE}
elif [[ $STOREDIR =~ /afs/.* ]]; then
    cmsMkdir ${STOREDIR}
    cp -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${OUTFILE}
fi
