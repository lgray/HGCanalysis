import ROOT

#transverse momentum
pt=ROOT.TGraph()
pt.SetName("pt")
minPt,maxPt=1.,100.
nPts=1000
dPt=(maxPt-minPt)/nPts
for i in xrange(0,nPts): pt.SetPoint(i,minPt+i*dPt,1.)

#rapidity
rapidity=ROOT.TGraph()
rapidity.SetName("rapidity")
minY,maxY=1.5,3.0
nPts=1000
dY=(maxY-minY)/nPts
for i in xrange(0,nPts): rapidity.SetPoint(i,minY+i*dY,1.)

#save to file
fOut=ROOT.TFile.Open('flatptygun.root','RECREATE')
fOut.cd()
pt.Write()
rapidity.Write()
fOut.Close()
