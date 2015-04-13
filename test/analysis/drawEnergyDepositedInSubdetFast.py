import ROOT
fIn=ROOT.TFile.Open('RelValSingleMuPt100Extended_SimHits.root')
HGC=fIn.Get('analysis/HGC')

"""
Material overburden and e.m. scales to be applied
"""
def getWeightingScheme(mode):
    weightingScheme=None
    if mode=='lambda':
        weightingScheme={
            ('em_EE',0 ) :  0.229437, 
            ('EE',1 ):  0.0136, ('EE',2 ):  0.0461, ('EE',3 ):  0.0448, ('EE',4 ):  0.0241, ('EE',5 ):  0.0448, ('EE',6 ):  0.0241, ('EE',7 ):  0.0448, ('EE',8 ):  0.0241, ('EE',9 ):  0.0448, ('EE',10):  0.0241,
            ('EE',11):  0.0448, ('EE',12):  0.0347, ('EE',13):  0.0511, ('EE',14):  0.0347, ('EE',15):  0.0511, ('EE',16):  0.0347, ('EE',17):  0.0511, ('EE',18):  0.0347, ('EE',19):  0.0511, ('EE',20):  0.0347,
            ('EE',21):  0.0511, ('EE',22):  0.0488, ('EE',23):  0.0642, ('EE',24):  0.0488, ('EE',25):  0.0642, ('EE',26):  0.0488, ('EE',27):  0.0642, ('EE',28):  0.0488, ('EE',29):  0.0642, ('EE',30):  0.0488,
            ('em_HEF',0 ):  0.17804237, 
            ('HEF',1):  0.3377, ('HEF',2):  0.2727,
            ('em_HEB',0 ):  0.2406134,
            ('HEB',1):  0.4760
            }
    elif mode=='dedx':
        weightingScheme={
            ('em_EE',0 ) :  0.00179,
            ('EE',1 ):  2.372,  ('EE',2 ):  9.541,  ('EE',3 ):  8.816,  ('EE',4 ):  5.125,  ('EE',5 ):  8.816,  ('EE',6 ):  5.125,  ('EE',7 ):  8.816,  ('EE',8 ):  5.125,  ('EE',9 ):  8.816,  ('EE',10):  5.125,
            ('EE',11):  8.816,  ('EE',12):  7.445,  ('EE',13):  10.217, ('EE',14):  7.445,  ('EE',15):  10.217, ('EE',16):  7.445,  ('EE',17):  10.217, ('EE',18):  7.445,  ('EE',19):  10.217, ('EE',20):  10.539,
            ('EE',21):  10.217, ('EE',22):  10.539, ('EE',23):  13.148, ('EE',24):  10.539, ('EE',25):  13.148, ('EE',26):  10.539, ('EE',27):  13.148, ('EE',28):  10.539, ('EE',29):  13.148, ('EE',30):  10.539,
            ('em_HEF',0 ):  0.000918116,
            ('HEF',1):  65.001, ('HEF',2):  52.954,
            ('em_HEB',0 ):  0.00123365,
            ('HEB',1):  92.196
            }
    else:
        weightingScheme={
            ('em_EE',0 ) :  0.0121150,
            ('EE',1 ):  0.0798, ('EE',2 ):  0.9214, ('EE',3 ):  0.5960, ('EE',4 ):  0.5691, ('EE',5 ):  0.5960, ('EE',6 ):  0.5691, ('EE',7 ):  0.5960, ('EE',8 ):  0.5691, ('EE',9 ):  0.5960, ('EE',10):  0.5691,
            ('EE',11):  0.5960, ('EE',12):  0.8687, ('EE',13):  0.7920, ('EE',14):  0.8687, ('EE',15):  0.7920, ('EE',16):  0.8687, ('EE',17):  0.7920, ('EE',18):  0.8687, ('EE',19):  0.7920, ('EE',20):  0.8687,
            ('EE',21):  0.7920, ('EE',22):  1.2683, ('EE',23):  1.2019, ('EE',24):  1.2683, ('EE',25):  1.2019, ('EE',26):  1.2683, ('EE',27):  1.2019, ('EE',28):  1.2683, ('EE',29):  1.2019, ('EE',30):  1.2683,
            ('em_HEF',0 ):  0.01581528,
            ('HEF',1):  3.5803, ('HEF',2):  3.1029,
            ('em_HEB',0 ):  0.021909,
            ('HEB',1):  5.2279
            }
    return weightingScheme

#prepare histograms
#weightingScheme=getWeightingScheme('lambda')
#weightingScheme=getWeightingScheme('dedx')
weightingScheme=getWeightingScheme('x0')
histos={}
for key in ['ee','hef','heb']:
    histos[key]=ROOT.TH1F(key,';e.m. energy [GeV];Events (a.u.)',50,0,4)
    histos[key].SetDirectory(0)
    histos[key].Sumw2()

#fill histograms
for ev in xrange(0,HGC.GetEntriesFast()):
    HGC.GetEntry(ev)

    e_EE=0
    for i in xrange(0,30):
        ien=getattr(HGC,'edep_rec')[i]
        if ien<0.5: continue
        emWeight = weightingScheme[('em_EE',0)]
        weight   = weightingScheme[('EE',i+1)]
        e_EE    += ien*weight*emWeight
    if e_EE>0 : histos['ee'].Fill(e_EE)

    e_HEF=0
    for i in xrange(30,42):
        ien=getattr(HGC,'edep_rec')[i]
        if ien<0.5 : continue
        emWeight = weightingScheme[('em_HEF',0)]
        weight   = weightingScheme[('HEF',2)]
        if i==30 : weight = weightingScheme[('HEF',1)]
        e_HEF    += ien*weight*emWeight
    if e_HEF>0 : histos['hef'].Fill(e_HEF)

    e_HEB=0
    for i in xrange(42,54):
        ien=getattr(HGC,'edep_rec')[i]
        if ien<0.5 : continue
        emWeight = weightingScheme[('em_HEB',0)]
        weight   = weightingScheme[('HEB',1)]
        e_HEB    += ien*weight*emWeight
    if e_HEB>0 : histos['heb'].Fill(e_HEB)

fIn.Close()


from UserCode.HGCanalysis.PlotUtils import *
customROOTstyle()
ROOT.gROOT.SetBatch(False)
ROOT.gStyle.SetOptFit(1111)
c=ROOT.TCanvas('c','c',1500,500)
c.Divide(3,1)
c.cd(1)
histos['ee'].Scale(1./histos['ee'].Integral())
histos['ee'].Fit('landau')
histos['ee'].Draw('histsame')
MyPaveText('#bf{CMS} #it{simulation} EE')
c.cd(2)
histos['hef'].Scale(1./histos['hef'].Integral())
histos['hef'].Fit('landau')
histos['hef'].Draw('histsame')
MyPaveText('FH')
c.cd(3)
histos['heb'].Scale(1./histos['heb'].Integral())
histos['heb'].Fit('landau')
histos['heb'].Draw('histsame')
MyPaveText('BH')
raw_input()
