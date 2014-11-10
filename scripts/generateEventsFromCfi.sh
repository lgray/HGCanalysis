#!/bin/bash


#
# PARSE CONFIGURATION PARAMETERS
#
PID=211
ENERGY=25
NEVENTS=100
CFI=UserCode/HGCanalysis/python/particleGun_cfi.py
WORKDIR="/tmp/`whoami`/"
STOREDIR=${WORKDIR}
JOBNB=1
GEOMETRY="Extended2023HGCalMuon,Extended2023HGCalMuonReco"
PILEUP=""
TAG=""
PILEUPINPUT=root://eoscms//eos/cms/store/cmst3/group/hgcal/CMSSW/MinBias_CMSSW_6_2_X_SLHC_2014-09-10-0200/
PHYSLIST="QGSP_FTFP_BERT_EML"

while getopts "hp:e:n:c:o:w:j:g:ut:i:l:" opt; do
    case "$opt" in
    h)
        echo ""
        echo "generateEventsFromCfi.sh [OPTIONS]"
	echo "     -p      pdg ids to generate (csv)"
	echo "     -e      energy to generate"
	echo "     -n      number of events go generate"
	echo "     -c      cfi to use with cmsDriver command"
	echo "     -o      output directory (local or eos)"
	echo "     -w      local work directory (by default /tmp/user)"
	echo "     -l      GEANT4 physics list: QGSP_FTFP_BERT_EML (default) / FTFP_BERT_EML / FTFP_BERT_XS_EML / QBBC"
        echo "     -j      job number"
	echo "     -g      geometry"
	echo "             v4:            Extended2023HGCalV4Muon,Extended2023HGCalV4MuonReco"
	echo "             v5 (default) : ${GEOMETRY}"
	echo "     -u      Turn on Pileup"
        echo "     -i      pileup input file"
        echo "     -t      tag to name output file"
	echo "     -h      help"
        echo ""
	exit 0
        ;;
    l)  PHYSLIST=$OPTARG
	;;
    p)  PID=$OPTARG
        ;;
    e)  ENERGY=$OPTARG
        ;;
    n)  NEVENTS=$OPTARG
        ;;
    c)  CFI=$OPTARG
        ;;
    o)  STOREDIR=$OPTARG
	;;
    w)  WORKDIR=$OPTARG
	;;
    j)  JOBNB=$OPTARG
	;;
    g)  GEOMETRY=$OPTARG
	;;
    u)  PILEUP="--pileup AVE_140_BX_25ns"
        ;;
    t)  TAG=$OPTARG
        ;;
    i)  PILEUPINPUT=$OPTARG
        ;;
    esac
done

#
# CONFIGURE JOB
#
if [ "$TAG" = "" ]; then
   TAG=${PID}_${ENERGY}
fi
BASEJOBNAME=HGCEvents_${TAG}_${JOBNB}
BASEJOBNAME=${BASEJOBNAME/","/"_"}
OUTFILE=${BASEJOBNAME}.root
PYFILE=${BASEJOBNAME}_cfg.py
LOGFILE=${BASEJOBNAME}.log

if [ "$PILEUP" = "local" ]; then
    PILEUP="${PILEUP} --pileup_input ${PILEUPINPUT}"
fi

#run cmsDriver
#cmsDriver.py ${CFI} -n ${NEVENTS} \
#    --python_filename ${WORKDIR}/${PYFILE} --fileout file:${WORKDIR}/${OUTFILE} \
#    $PILEUP \
#    -s GEN,SIM --datatier GEN-SIM --eventcontent FEVTDEBUGHLT\
#    --conditions auto:upgradePLS3 --beamspot Gauss --magField 38T_PostLS1 \
#    --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCalMuon \
#    --geometry ${GEOMETRY} \
#    --no_exec 

cmsDriver.py ${CFI} -n ${NEVENTS} \
    --python_filename ${WORKDIR}/${PYFILE} --fileout file:${WORKDIR}/${OUTFILE} \
    $PILEUP \
    -s GEN,SIM,DIGI:pdigi_valid,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO --datatier GEN-SIM-DIGI-RECO --eventcontent FEVTDEBUGHLT\
    --conditions auto:upgradePLS3 --beamspot Gauss --magField 38T_PostLS1 \
    --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCalMuon \
    --geometry ${GEOMETRY} \
    --no_exec 

echo "cmsDriver.py ${CFI} -n ${NEVENTS} \
    --python_filename ${WORKDIR}/${PYFILE} --fileout file:${WORKDIR}/${OUTFILE} \
    $PILEUP \
    -s GEN,SIM,DIGI:pdigi_valid,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO --datatier GEN-SIM-DIGI-RECO --eventcontent FEVTDEBUGHLT\
    --conditions auto:upgradePLS3 --beamspot Gauss --magField 38T_PostLS1 \
    --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCalMuon \
    --geometry ${GEOMETRY} \
    --no_exec"

#customize with values to be generated
echo "process.g4SimHits.StackingAction.SaveFirstLevelSecondary = True" >> ${WORKDIR}/${PYFILE}
echo "process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(${JOBNB})" >> ${WORKDIR}/${PYFILE}
echo "process.source.firstEvent=cms.untracked.uint32($((NEVENTS*(JOBNB-1)+1)))" >> ${WORKDIR}/${PYFILE}
if [ -z ${PHYSLIST} ]; then
    echo "Changing default Geant4 physics list to $PHYSLIST"
    echo "process.g4SimHits.Physics.type=cms.string('SimG4Core/Physics/${PHYSLIST}')" >> ${WORKDIR}/${PYFILE}
fi
SUBSTRING="s/MinE = cms.double(0)/MinE = cms.double(${ENERGY})/"
SUBSTRING="${SUBSTRING};s/MaxE = cms.double(0)/MaxE = cms.double(${ENERGY})/"
SUBSTRING="${SUBSTRING};s/ParticleID = cms.vint32(0)/ParticleID = cms.vint32(${PID})/"
sed -i.bak "${SUBSTRING}" ${WORKDIR}/${PYFILE}
rm ${WORKDIR}/${PYFILE}.bak

#
# RUN cmsRun and at the end move the output to the required directory
#
cmsRun ${WORKDIR}/${PYFILE} > ${WORKDIR}/${LOGFILE} 2>&1

if [[ $STOREDIR =~ .*/store/cmst3.* ]]; then
    cmsMkdir ${STOREDIR}
    cmsStage -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${OUTFILE}
elif [[ $STOREDIR =~ /afs/.* ]]; then
    cmsMkdir ${STOREDIR}
    cp -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${OUTFILE}
fi


echo "Generated $NEVENTS events for pid=$PID with E=$ENERGY GeV"
echo "Local output @ `hostname` stored @ ${WORKDIR}/${OUTFILE} being moved to ${STOREDIR}" 
echo "cmsRun cfg file can be found in ${WORKDIR}/${PYFILE}"
echo "log file can be found in ${WORKDIR}/${LOGFILE}"
