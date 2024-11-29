import ROOT
import random
import math

#get histogram of already existing events
f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/bkgMC_test38_ALL.root")

hold = f.Get("ntuples/allGenPtEta") 
hold.Rebin(100)

hnew = ROOT.TH1F("hnew", "hnew", 100, 0.0, 100.0)
for xbin in range(hold.GetNbinsX()):
    hnew.SetBinContent(xbin, hold.GetBinContent(xbin)) 

#x, y arrays for the final TGraph
xlist = []
ylistMed = []
ylistMax = []

#maximum number of events to simulate
maxE = 50000000
fMaxPt = 65.0
fMinPt = 10.0
xpeak = 18
xmax = 10

for i in range(maxE+1):
    x = random.random()
    pt = ( min(-1/10.0 * math.log( x ), 1.0) ) * (fMaxPt - fMinPt) + fMinPt
    hnew.Fill(pt)

    #calculate the weights every million
    if i%1000000 == 0:
        peakWt = 6.3e11*3.1*38/((xpeak+0.5)**6 * hnew.GetBinContent(xpeak+1))
        maxWt = 6.3e11*3.1*38/((xmax+0.5)**6 * hnew.GetBinContent(xmax+1))
        print("i=%d, peakWt=%f, maxWt=%f"%(i, peakWt, maxWt)) 
        xlist.append(i)
        ylistMed.append(peakWt)
        ylistMax.append(maxWt) 

import array
xarr = array.array('f', xlist)
yarrMed = array.array('f', ylistMed)
yarrMax = array.array('f', ylistMax) 
#graph of peak weight vs number of new events generated
tgMed = ROOT.TGraph(len(xlist), xarr, yarrMed)
tgMed.SetName("tgMed")
tgMax = ROOT.TGraph(len(xlist), xarr, yarrMax)
tgMax.SetName("tgMax")
tgMax.GetXaxis().SetTitle("New events generated")
tgMax.GetYaxis().SetTitle("Estimated Weights")
tgMax.SetLineWidth(2)
tgMax.Draw() 
tgMed.SetLineWidth(2)
tgMed.SetLineColor(ROOT.kRed)
tgMed.Draw("same") 
leg = ROOT.TLegend()
leg.AddEntry(tgMax, "Highest Weight")
leg.AddEntry(tgMed, "Median Weight")
leg.Draw("same")

input("h")
