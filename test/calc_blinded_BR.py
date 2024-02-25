import ROOT
import sys
from array import array
sys.path.insert(1, "tm_analysis/analysis/python/utils")
import blinding

lumi = 38.48 #fb^-1 (all 2022)

selection = "nominal"
#selection = "loose"
#selection = "tight"

#Get 2mu2e acceptance!
if selection == "nominal":
    f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3837.root")
    #get this from the fit to data
    N2mu2e = 131.9
    Err2mu2e = 14.9
elif selection == "loose":
    f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3838.root")
    N2mu2e = 766.5
elif selection == "tight":
    f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3839.root")
    N2mu2e = 32.7

ptGen2mu2e1 = f2mu2e.Get("hpTGenAll")
#ptGen2mu2e.Rebin(5)
ptAcc2mu2e1 = f2mu2e.Get("hpTGenAccmmelel")
#ptAcc2mu2e = f2mu2e.Get("hpTGenAccmmlplp")
#bins should be 1 GeV wide!!!
#ptAcc2mu2e.Rebin(5)
#make pT bins for new histogram (at low pT, width is 1 GeV; then gradually widen as pT goes up)
newbins = [i*1.0 for i in range(5, 31)]
for i in range(32, 41, 2):
    newbins.append(1.0*i)
for i in range(45, 56, 5):
    newbins.append(1.0*i)
newbins.append(70.0)
newbins.append(100.0)
newbins.append(101.0)
print("new pT bins: " + str(newbins)) 
#last entry in newbins is the upper limit of the last bin, not a new bin
newnptbins = len(newbins)-1
ptGen2mu2e = ROOT.TH1F("ptGen2mu2e", "ptGen2mu2e", newnptbins, array('d', newbins)) 
ptAcc2mu2e = ROOT.TH1F("ptAcc2mu2e", "ptAcc2mu2e", newnptbins, array('d', newbins)) 
for i in range(ptAcc2mu2e1.GetNbinsX()):
    ptGen2mu2e.Fill(ptGen2mu2e1.GetBinCenter(i), ptGen2mu2e1.GetBinContent(i)) 
    ptAcc2mu2e.Fill(ptAcc2mu2e1.GetBinCenter(i), ptAcc2mu2e1.GetBinContent(i)) 
ptAcc2mu2e.Sumw2()
ptAcc2mu2e.Divide(ptGen2mu2e)
if selection == "nominal":
    xsfname = "xsec2022_nominal.root"
else:
    xsfname = "xsec2022.root"
xsecf = ROOT.TFile.Open(xsfname)
xsec = xsecf.Get("hXsecCor") 

#use the acceptance from the sigMC file, NOT from the xsec
acc2mu2e = ptAcc2mu2e

c1,c2 = blinding.get_blind_coeffs()

#denominator is a sum over all pT bins
denom = 0.0
for i in range(1, xsec.GetNbinsX()):
    denom += (c1*xsec.GetBinContent(i)*lumi*c2*acc2mu2e.GetBinContent(i))

#multiply by 1e-6 to make the number more meaningful
br = N2mu2e / (denom*10**-6)
#just statistical uncertainty for now!!
#unct = N2mu2e**0.5 / (denom*10**-6) 
unct = Err2mu2e / (denom*10**-6) 

print("%s Blinded BR: (%f +/- %f (stat))e-6"%(selection, br, unct)) 
