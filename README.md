# HGCanalysis scripts

## Installation 

To add to your cmssw area

git clone git@github.com:PFCal-dev/HGCanalysis UserCode/HGCanalysis

To update in the remote

git push git@github.com:PFCal-dev/HGCanalysis

## GEN-SIM-RECO production

Use the generateEventsFromCfi.sh to steer the generation.
It will call cmsDriver.py and customize the configuration file.
The options can be inspected by calling:

generateEventsFromCfi.sh -h

A full production can be ran locally or submitted to the batch using 
the submitLocalHGCalProduction.py wrapper script. Two examples are given below:

### Particle gun 

For muons it's enough a single energy

#default geometry
python scripts/submitLocalHGCalProduction.py -q 1nd -n 10 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single13_${CMSSW_VERSION} -p 13 -n 100 -e 100";
#change geometry scenario
python scripts/submitLocalHGCalProduction.py -q 1nd -n 100 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single13_v4_${CMSSW_VERSION} -p 13 -n 100 -e 100 -g Extended2023HGCalV4Muon,Extended2023HGCalV4MuonReco";

For regression use flat gun

python scripts/submitLocalHGCalProduction.py -q 2nw -n 2500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/FlatPtYSingle22_${CMSSW_VERSION}/RECO_a -c UserCode/HGCanalysis/python/particlePtYGun_cfi.py -n 200 -p 22 -f";
python scripts/submitLocalHGCalProduction.py -q 1nw -n 2500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/FlatPtYSingle22_${CMSSW_VERSION}/RECO_b -c UserCode/HGCanalysis/python/particlePtYGun_cfi.py -n 200 -p 22 -f";

# DECEMBER JAMBOREE PRODUCTION

## SIM

energies=(10 20 40 50 75 100 125 175 250 400 500)
pids=(22 211 11)
pu=(140 0 200)
for p in ${pu[@]}; do
for pid in ${pids[@]}; do
    for en in ${energies[@]}; do
      python scripts/submitLocalHGCalProduction.py -q 1nd -n 500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/ -p ${pid} -n 25 -e ${en} -s -a ${p}";
     done
done
done

pids=(111 15)
for p in ${pu[@]}; do
for pid in ${pids[@]}; do
    for en in ${energies[@]}; do 
    	python scripts/submitLocalHGCalProduction.py -n 500 -q 2nd -s generateEventsFromCfi.sh -o "-c UserCode/HGCanalysis/python/jetGun_cfi.py -o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/ -p ${pid} -n 25 -e ${en} -s -a ${p}"; 
    done
done
done

for p in ${pu[@]}; do
    python scripts/submitLocalHGCalProduction.py -q 2nw -n 1000 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/VBFtoH125toTauTau_${CMSSW_VERSION} -c UserCode/HGCanalysis/python/VBFH125toTauTau_cfi.py -n 25 -p 25 -s -a ${p}";
    python scripts/submitLocalHGCalProduction.py -q 1nw -n 1000 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/QCDFlatPt15to3000_${CMSSW_VERSION} -c UserCode/HGCanalysis/python/QCDForPF_14TeV_cfi.py -n 25 -p 1 -s -a ${p}";
done



#test alternative physics lists for pions
phys=("QGSP_FTFP_BERT_EML" "FTFP_BERT_EML" "FTFP_BERT_XS_EML" "QBBC")
pids=(211)
energies=(30)
for pid in ${pids[@]}; do
    for en in ${energies[@]}; do
    	for p in ${phys[@]}; do
        python scripts/submitLocalHGCalProduction.py -q 1nd -n 10 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}_${p} -p ${pid} -n 400 -e ${en} -l ${p}";
        done
    done
done


#
# PRODUCTION WITH PILEUP
#
#simulation step
beamspots=(HLLHC_Fix HLLHCCrabKissing)
for b in ${beamspots[@]}; do
    python scripts/submitLocalHGCalProduction.py -q 2nd -n 500 -s generateEventsFromCfi.sh -o "-s -o /store/cmst3/group/hgcal/CMSSW/MinBias_${CMSSW_VERSION}/${b} -b ${b} -c UserCode/HGCanalysis/python/minBias_cfi.py -n 500";
    python scripts/submitLocalHGCalProduction.py -q 2nd -n 1000 -s generateEventsFromCfi.sh -o "-s -o /store/cmst3/group/hgcal/CMSSW/GluGluHtoGG_${CMSSW_VERSION}/${b}/SIM -b ${b} -c UserCode/HGCanalysis/python/hToGG_cfi.py -n 50 -p 22";
    python scripts/submitLocalHGCalProduction.py -q 1nd -n 1000 -s generateEventsFromCfi.sh -o "-s -o /store/cmst3/group/hgcal/CMSSW/VBFHtoInv_${CMSSW_VERSION}/${b}/SIM -b ${b} -c UserCode/HGCanalysis/python/VBFH125toInv_cfi.py -n 50 -p 0";
done

#SIM integrity check
samples=(VBFHtoInv GluGluHtoGG MinBias)
for b in ${beamspots[@]}; do
    for s in ${samples[@]}; do
    	python scripts/checkProductionIntegrity.py -d /store/cmst3/group/hgcal/CMSSW/${s}_${CMSSW_VERSION}/${b}/SIM &
    done
done

#digi step
pu=(0 140)
samples=(VBFHtoInv GluGluHtoGG)
for s in ${samples[@]}; do
    for p in ${pu[@]}; do
    	for b in ${beamspots[@]}; do
	    sigDir=${s}_${CMSSW_VERSION}/${b}
	    inputFiles=(`cmsLs /store/cmst3/group/hgcal/CMSSW/${sigDir}/SIM | awk '{print $5}'`)
            nFiles=${#inputFiles[@]};
	    echo "Submitting ${nFiles} jobs for events in ${sigDir}"
	    python scripts/submitLocalHGCalProduction.py -n ${nFiles} -q 2nw -s digitizeAndMix.sh -o "-o /store/cmst3/group/hgcal/CMSSW/${sigDir}/DIGI-PU${p} -c  ${CMSSW_BASE}/src/UserCode/HGCanalysis/test/digitizeAndMix_cfg.py -m MinBias_${CMSSW_VERSION}/${b}/SIM -p ${p} -i ${sigDir}/SIM";
        done
    done
done

#DIGI integrity check
for s in ${samples[@]}; do
    for p in ${pu[@]}; do
    	for b in ${beamspots[@]}; do
    	    python scripts/checkProductionIntegrity.py -d /store/cmst3/group/hgcal/CMSSW/${s}_${CMSSW_VERSION}/${b}/DIGI-PU${p} &	  
        done
    done
done
    
#RECO step

## Producing analysis ntuples

The ntuples are produced by plugins/HGCSimHitsAnalyzer.cc.  Change the code according to your needs.
To submit the production of the ntuples you can use the following script (it will printout the options)

cmsRun runHGCSimHitsAnalyzer_cfg.py

Submit several jobs to the batch and store the output in EOS
tags=("Single211_CMSSW_6_2_0_SLHC23_patch1" "Single22_CMSSW_6_2_0_SLHC23_patch1" "Single130-FixE_CMSSW_6_2_0_SLHC23_patch2")
#subdir="RECO-PU0"
subdir="pandoraRECO";
#subdir="mqRECO"
for tag in ${tags[@]}; do
    #python scripts/checkProductionIntegrity.py -d /store/cmst3/group/hgcal/CMSSW/${tag}/${subdir}; 
    inputFiles=(`cmsLs /store/cmst3/group/hgcal/CMSSW/${tag}/${subdir} | awk '{print $5}'`);
    nFiles=${#inputFiles[@]};	
    nJobs=$((nFiles/50));
    for i in `seq 0 ${nJobs}`; do
     	startFile=$((i*=50));
   	cmsRun test/runHGCSimHitsAnalyzer_cfg.py ${tag}/${subdir} ${startFile} 50 & 
	nrunning="`ps -u psilva | grep -ir cmsRun | wc -l`"
	while [ $nrunning -gt 5 ]; do
	      echo "$nrunning cmsRun processes launched, sleeping 30s"
	      sleep 30;
	      nrunning="`ps -u psilva | grep -ir cmsRun | wc -l`"
	done
    done
done

### JET ANALYZER
a=(pandoraRECO-PU140 pandoraRECO-PU0)
b=(/store/cmst3/group/hgcal/CMSSW/RelValQCD_Pt_80_120_14TeV /store/cmst3/group/hgcal/CMSSW/RelValQCDForPF_14TeV)
for i in ${a[@]}; do 
    for j in ${b[@]}; do 
    	cmsRun test/runHGCJetsAnalyzer_cfg.py ${j}/${i} "`basename ${j}`_${i}.root" &
    done
done

a=(`ls *.root`)
for i in ${a[@]}; do
   python test/analysis/showJetAnalysisPlots.py -i ${i} -o ~/public/html/HGCal/
done
