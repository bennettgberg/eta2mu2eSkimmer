import ROOT
import sys
from array import array
sys.path.insert(1, "tm_analysis/analysis/python/utils")
import blinding

#calc relative BR (B2mu2e/B2mu) instead of absolute? (uncertainties cancel then)
relative = False

year = 2022

if year == 2022:
    lumi = 38.48 #fb^-1 (all 2022)
elif year == 2023:
    lumi = 28.89 #fb^-1 (all 2023)

selection = "nominal"
#selection = "loose"
#selection = "tight"
#selection = "vloose"

#2mu sig model (nominal: DG, Cheb4)
refMod = "DG" 
refBkg = "Cheb4"

#refMod = "Voigtian"
#refBkg = "Cheb3"

new_wt = 3 #True

if year == 2022:
    #Get 2mu2e acceptance!
    if selection == "nominal" and not new_wt:
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3837.root")
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3860.root")
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3866.root")
        #get this from the fit to data
        #N2mu2e = 147.7
        #Err2mu2e = 14.0
        #DoubleGauss fit
        N2mu2e = 153.86
        Err2mu2e = 16.63
    elif selection == "nominal" and new_wt == 1:
    #new_wt: DG/Cheb4
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3876.root")
        #get this from the fit to data
        #DoubleGauss fit 155.262708 +/- 19.069812
        N2mu2e = 155.00 #26
        Err2mu2e = 18.57 #19.07
    elif selection == "nominal" and new_wt == 2:
        #155.088371 +/- 17.160677
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3879.root")
        N2mu2e = 155.09
        Err2mu2e = 17.16
    elif selection == "nominal" and new_wt == 3:
        #154.424684 +/- 24.088907 (stat.)
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3880.root")
        N2mu2e = 154.42
        Err2mu2e = 24.09
    elif selection == "loose":
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3838.root")
        N2mu2e = 766.5
    elif selection == "tight":
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3839.root")
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3871.root")
        #PU corex
        if new_wt:
            f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3875.root")
            #38.227636 +/- 12.546614
            N2mu2e = 38.23
            Err2mu2e = 12.55
        else:
            f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3872.root")
            #N2mu2e = 36.05
            #Err2mu2e = 6.83
            #DoubleGauss fit
            N2mu2e = 38.75
            Err2mu2e = 7.80
    elif selection == "vloose":
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3858.root")
        N2mu2e = 2260.1
        Err2mu2e = 72.9
elif year == 2023:
    f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2023/ultraskimmed/bparking_2023sigMCtest430.root")
    N2mu2e = 0
    Err2mu2e = 0

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
##if selection == "nominal" :
#if True:
#    if refMod == "Voigtian" and refBkg == "Cheb3":
#        xsfname = "xsec2022_nominal.root"
#    else:
#        xsfname = "xsec2022_%s%s.root"%(refMod if refMod != "Voigtian" else "", refBkg if refBkg != "Cheb3" else "")
#else:
#    xsfname = "xsec2022.root"
xsfname = "xsec2022_%s%s.root"%("Voigt" if refMod == "Voigtian" else refMod, refBkg)
xsecf = ROOT.TFile.Open(xsfname)
xsec = xsecf.Get("hXsecCor") 

#use the acceptance from the sigMC file, NOT from the xsec
acc2mu2e = ptAcc2mu2e

c1,c2 = blinding.get_blind_coeffs()

#denominator is a sum over all pT bins
denom = 0.0
for i in range(1, xsec.GetNbinsX()):
    if relative:
        denom += (c1*c2*acc2mu2e.GetBinContent(i)/acc2mu.GetBinContent(i))
    else:
        denom += (c1*xsec.GetBinContent(i)*lumi*c2*acc2mu2e.GetBinContent(i))

#multiply by 1e-6 to make the number more meaningful
br = N2mu2e / (denom*10**-6)
#just statistical uncertainty for now!!
#unct = N2mu2e**0.5 / (denom*10**-6) 
unct = Err2mu2e / (denom*10**-6) 

print("%s %s (2mu fit: %s, %s) Blinded BR: (%f +/- %f (stat))e-6"%(selection, "with new weights" if new_wt else "", refMod, refBkg, br, unct)) 
