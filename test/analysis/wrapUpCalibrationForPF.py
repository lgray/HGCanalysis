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
		theDir='edep_rec'
		if pid==22 : 
			theDir='edep_rec--lambdaWeighting'
			key='e_em'
			if prod.find('-EE_HEF_AIR')>=0  : key += '_BH'
			elif prod.find('-EE_AIR')>=0    : key += '_FH'
			else                            : key += '_EE'

			fIn=ROOT.TFile.Open('%s/src/UserCode/HGCanalysis/Single%d_%s_%s_SimHits/%s/calib_uncalib.root'%(CMSSW_BASE,pid,CMSSW_VERSION,prod,theDir))
			func=fIn.Get('simple_calib')
			calibs[key]=()
			calibs[key]=calibs[key]+(func.GetExpFormula(),)
			for ip in xrange(0,func.GetNpar()):calibs[key]=calibs[key]+(func.GetParameter(ip),)
			fIn.Close()

		#pi/e
		else:
			if prod=='RECO-PU0' : theDir='edep_rec-IndComb'
			piOverEfile, piovereGr,key = 'EEHEFHEB_response_corrx_corry.root','EE_responseFunc','pioe_EE'
			if prod.find('-EE_AIR')>=0     : piOverEfile, piovereGr,key = 'HEFHEB_response_corry.root','HEF_responseFunc','pioe_HEF'
			if prod.find('-EE_HEF_AIR')>=0 : piOverEfile, piovereGr,key = 'HEB_response.root','HEB_responseFunc','pioe_HEB'
			fIn=ROOT.TFile.Open('%s/src/UserCode/HGCanalysis/Single%d_%s_%s_SimHits/%s/%s'%(CMSSW_BASE,pid,CMSSW_VERSION,prod,theDir,piOverEfile))
			func=fIn.Get(piovereGr).Clone(key)
			calibs[key]=()
			calibs[key]=calibs[key]+(func.GetExpFormula(),)
			for ip in xrange(0,func.GetNpar()):calibs[key]=calibs[key]+(func.GetParameter(ip),)
			if prod=='RECO-PU0':
				for i in xrange(0,3):
					func=fIn.Get('resEvolFunc_%d'%i)
					key='pioe_bananaParam_%d'%i
					calibs[key]=(func.GetExpFormula(),)
					for ip in xrange(0,func.GetNpar()):calibs[key]=calibs[key]+(func.GetParameter(ip),)
			fIn.Close()

#dump xml snippet 
print '******* XML snippet for Pandora ********'
print '<CMSGlobalHadronCompensation>'
print '\t<MipEnergyThreshold>10.0</MipEnergyThreshold>'
for key in calibs :
	print '\t<!-- %s function: '%key,calibs[key][0],' -->'
	print '\t<%s>'%key,
	for i in xrange(1,len(calibs[key])): 
		print calibs[key][i],
		if i<len(calibs[key])-1 : print ',',
	print '</%s>'%key
print '</CMSGlobalHadronCompensation>'

