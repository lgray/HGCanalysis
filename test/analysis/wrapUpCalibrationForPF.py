import ROOT
import os
from datetime import datetime
 
i = datetime.now() 
calibHash=i.strftime('%Y-%m-%d')

CMSSW_BASE    = os.environ['CMSSW_BASE']
CMSSW_VERSION = os.environ['CMSSW_VERSION']
pids          = [22,211]
prods         = ['RECO-PU0-EE_HEF_AIR','RECO-PU0-EE_AIR','RECO-PU0']
calibs        = {}
for pid in pids:
	for prod in prods:

		#baseline e-scale and residual pi/e
		theDirs=['edep_rec']
		if pid==22 : theDirs.append('edep_rec--lambdaWeighting')
		for theDir in theDirs:
			fIn=ROOT.TFile.Open('%s/src/UserCode/HGCanalysis/Single%d_%s_%s_SimHits/%s/calib_uncalib.root'%(CMSSW_BASE,pid,CMSSW_VERSION,prod,theDir))
			key='e'
			if theDir.find('lambda')<0  and key=='e' : key='e_em'
			if pid==211                              : key  = 'pi'
			if prod.find('-EE_HEF_AIR')>=0           : key += '_BH'
			elif prod.find('-EE_AIR')>=0             : key += '_FH'
			else                                     : key += '_EE'
			calibs[key]=fIn.Get('simple_calib').Clone(key)
			fIn.Close()

		#pi/e
		if pid!=211: continue
		if prod.find('-EE_AIR')>=0 : continue
		piOverEfile, piovereGr,key = 'EEHEFHEB_response_corry.root','EEHEF_responseFunc','pioe_Si'
		if prod.find('-EE_HEF_AIR')>=0 : piOverEfile, piovereGr,key = 'HEB_response.root','HEB_responseFunc','pioe_Sci'
		fIn=ROOT.TFile.Open('%s/src/UserCode/HGCanalysis/Single%d_%s_%s_SimHits/edep_rec/%s'%(CMSSW_BASE,pid,CMSSW_VERSION,prod,piOverEfile))
		calibs[key]=fIn.Get(piovereGr).Clone(key)
		fIn.Close()

#dump to file
fOut=ROOT.TFile.Open('HGCRecHitCalib_%s.root'%calibHash,'RECREATE')
for key in calibs : calibs[key].Write()
fOut.Close()
print 'Energy scale functions written to : %s'%fOut.GetName()
