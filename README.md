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

### Particle gun (100 events per file x 100 jobs) 

energies=(5 10 20 30 50 75 100 150 250 500)
pids=(13 11 211)
for pid in ${pids[@]}; do
for en in ${energies[@]}; do
	#default geometry
        python scripts/submitLocalHGCalProduction.py -q 1nd -n 100 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION} -p ${pid} -n 100 -e ${en}";
	#change geometry scenario
	#python scripts/submitLocalHGCalProduction.py -q 1nd -n 100 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_v4_${CMSSW_VERSION} -p ${pid} -n 100 -e ${en} -g Extended2023HGCalV4Muon,Extended2023HGCalV4MuonReco";
done
done

### Minimum bias (1000 events per file x 500 jobs, should be ok for later mixing with particle gun)

python scripts/submitLocalHGCalProduction.py -q 2nd -n 500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/MinBias_${CMSSW_VERSION} -c UserCode/HGCanalysis/python/minBias_cfi.py -n 500";

python scripts/submitLocalHGCalProduction.py -q 2nd -n 500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/MinBias_v4_${CMSSW_VERSION} -c UserCode/HGCanalysis/python/minBias_cfi.py -n 500 -g Extended2023HGCalV4Muon,Extended2023HGCalV4MuonReco";

### Other processes

Can use the minimum bias example, just substitute the argument passed in the -c option to point to the new cfi snippet.

### Redigitization with pileup mixing (will run one job per file, randomizing the min.bias files at start)

tau=(0 10 20)
pids=(13)
for tau in ${taus[@]}; do
    for pid in ${pids[@]}; do
        #inputFiles=(`cmsLs /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}`)
	#nFiles=${#inputFiles[@]}
        nFiles=5 
	python scripts/submitLocalHGCalProduction.py -n ${nFiles} -q 1nd -s redigitizeAndMix.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/tau_${tau} -t Single${pid}_${CMSSW_VERSION} -m MinBias_${CMSSW_VERSION} -p ${tau}";
done
done
done
    

## Producing analysis ntuples

The ntuples are produced by plugins/HGCSimHitsAnalyzer.cc. 
Change the code according to your needs.
To submit the production of the ntuples you can use the following script



pfs, CHECK ME THIS POINT FORWARD!

Submit ntuples production
for i in `seq 0 50 1950`; do
	python scripts/submitLocalAnalysis_cfg.py -t MinBias_v14 -q 2nd -f ${i} -s 50;
done

#calibration studies
python test/analysis//testSimHits.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples -t SingleElectron_SLHC13_30um_SimHits
python test/analysis//testSimHits.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples -t SingleK0L_SLHC13_30um_SimHits
python test/analysis//testSimHits.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples -t SinglePion_SLHC13_30um_SimHits

#occupancy studies
pileup=(140 200)
noise=(0) # it is generated in the summary
sdType=(0 1)
for p in ${pileup[@]}; do
for n in ${noise[@]}; do 
for s in ${sdType[@]}; do
python test/analysis/runOccupancyAnalysis.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples/ -t MinBias_v16 -p ${p} -s ${s} -n ${n} &
done
done
done

for p in ${pileup[@]}; do
for n in ${noise[@]}; do
for s in ${sdType[@]}; do
python test/analysis/drawOccupancyAnalysisSummary.py -i MinBias_v16_occ_pu${p}_sd${s}.root -o pu${p}_sd${s} &
done
done
done