import ROOT
import sys
from array import array
sys.path.insert(1, "tm_analysis/analysis/python/utils")
import blinding

year = 2022

if year == 2022:
    lumi = 38.48 #fb^-1 (all 2022)
elif year == 2023:
    lumi = 28.89 #fb^-1 (all 2023)

selection = "nominal"
#selection = "lowPt"
#selection = "nomuID"
#selection = "loose"
#selection = "tight"
#selection = "vloose"

#2mu sig model (nominal: DG, Cheb4)
refMod = "DG" 
#refMod = "Voigtian"

#refBkg = "Cheb2"
refBkg = "Cheb3"
#refBkg = "Cheb4"

new_wt = 1 #True

B2mu = 5.8e-6

if year == 2022:
    #Get 2mu2e acceptance!
    #if selection == "nominal" and not new_wt:
    if selection == "nomuID" and not new_wt:
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3837.root")
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3860.root")
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3866.root")
        #get this from the fit to data
        #N2mu2e = 147.7
        #Err2mu2e = 14.0
        #DoubleGauss fit
        N2mu2e = 153.86
        Err2mu2e = 16.63
    elif selection == "nomuID" and new_wt == 1:
    #new_wt: DG/Cheb4
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3876.root")
        #get this from the fit to data
        #DoubleGauss fit 155.262708 +/- 19.069812
        N2mu2e = 155.00 #26
        Err2mu2e = 18.57 #19.07
    elif selection == "nomuID" and new_wt == 2:
        #155.088371 +/- 17.160677
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3879.root")
        N2mu2e = 155.09
        Err2mu2e = 17.16
    elif selection == "nomuID" and new_wt == 3:
        #154.424684 +/- 24.088907 (stat.)
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3880.root")
        N2mu2e = 154.42
        Err2mu2e = 24.09
    elif selection == "nominal":
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3884.root")
        #N2mu2e = 135.17
        #Err2mu2e = 14.05
        #WITH triggercorrections: 136.664236 +/- 14.085243 (stat.)
        #Nominal frfrfrfrfrfrfr
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3888.root")
        #N2mu2e = 136.66
        #Err2mu2e = 14.09 
        #one trigger path only 135.832195 +/- 14.135855
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38110.root")
        #f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38110.root")
        #N2mu2e = 135.83
        #Err2mu2e = 14.14
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38130p474.root")
        #f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38130p474.root")
        #N2mu2e = 130.78
        #Err2mu2e = 13.62
        #124.148028 +/- 13.380399
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38137p478.root")
        #f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38137p478.root")
        #N2mu2e = 123.64
        #Err2mu2e = 13.39
        #N2mu2e = 128.36
        #Err2mu2e = 13.54
        #129.032638 +/- 13.570899
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38164p4733.root")
        #f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38137p478.root")
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38169p4735.root")
        #f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38169p4735.root")
        #N2mu2e = 129.03
        #Err2mu2e = 13.57
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38174p4740.root")
        f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38171p4737.root")
        N2mu2e = 124.23
        Err2mu2e = 13.38
    elif selection == "lowPt":
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38166p4734.root")
        f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38137p478.root")
        N2mu2e = 789.45
        Err2mu2e = 61.26
    elif selection == "loose":
        f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3838.root")
        N2mu2e = 766.5
    elif selection == "tight":
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3839.root")
        #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3871.root")
        #PU corex
        if new_wt:
            #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3875.root")
            #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38122.root")
            #f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38122.root")
            #f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38138p479.root")
            f2mu2e = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38147p4717.root")
            #f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38138p479.root")
            f2mu = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38139p4710.root")
            #38.227636 +/- 12.546614
            #N2mu2e = 38.23
            #Err2mu2e = 12.55
            #N2mu2e = 26.30
            #Err2mu2e = 5.75
            #25.546802 +/- 5.865937
            N2mu2e = 25.45
            Err2mu2e = 5.88
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
ptGen2mu2e1.SetName("ptGen2mu2e")
#ptGen2mu2e.Rebin(5)
#ptAcc2mu2e1 = f2mu2e.Get("hpTGenAccmmelel")
if selection == "lowPt":
    ptAcc2mu2e1 = f2mu2e.Get("hpTGenAccdRmmlplp")
else:
    ptAcc2mu2e1 = f2mu2e.Get("hpTGenAccdRmmelel")
ptGen2mu1 = f2mu.Get("hpTGenAll")
#ptAcc2mu1 = f2mu.Get("hpTGenAccmumu") 
ptAcc2mu1 = f2mu.Get("hpTGenAccdRmumu") 
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
ptGen2mu = ROOT.TH1F("ptGen2mu", "ptGen2mu", newnptbins, array('d', newbins)) 
ptAcc2mu = ROOT.TH1F("ptAcc2mu", "ptAcc2mu", newnptbins, array('d', newbins)) 
for i in range(ptAcc2mu2e1.GetNbinsX()):
    ptGen2mu2e.Fill(ptGen2mu2e1.GetBinCenter(i), ptGen2mu2e1.GetBinContent(i)) 
    ptAcc2mu2e.Fill(ptAcc2mu2e1.GetBinCenter(i), ptAcc2mu2e1.GetBinContent(i)) 
    ptGen2mu.Fill(ptGen2mu1.GetBinCenter(i), ptGen2mu1.GetBinContent(i)) 
    ptAcc2mu.Fill(ptAcc2mu1.GetBinCenter(i), ptAcc2mu1.GetBinContent(i)) 
ptAcc2mu2e.Sumw2()
ptAcc2mu2e.Divide(ptGen2mu2e)
ptAcc2mu.Sumw2()
ptAcc2mu.Divide(ptGen2mu)
##if selection == "nominal" :
#if True:
#    if refMod == "Voigtian" and refBkg == "Cheb3":
#        xsfname = "xsec2022_nominal.root"
#    else:
#        xsfname = "xsec2022_%s%s.root"%(refMod if refMod != "Voigtian" else "", refBkg if refBkg != "Cheb3" else "")
#else:
#    xsfname = "xsec2022.root"
if selection == "nomuID":
    xsfname = "xsec2022_%s%s.root"%("Voigt" if refMod == "Voigtian" else refMod, refBkg)
elif selection == "nominal" or selection == "lowPt":
    xsfname = "xsec2022_muID_%s%s.root"%("Voigt" if refMod == "Voigtian" else refMod, refBkg)
    #xsfname = "xsec2022.root"
elif selection == "tight":
    xsfname = "xsec2022_tightMuID_%s%s.root"%("Voigt" if refMod == "Voigtian" else refMod, refBkg)
xsecf = ROOT.TFile.Open(xsfname)
xsec = xsecf.Get("hXsecCor") 

#use the acceptance from the sigMC file, NOT from the xsec
acc2mu2e = ptAcc2mu2e
acc2mu = ptAcc2mu

c1,c2 = blinding.get_blind_coeffs()

#get N2mu from the raw yields
hN2mu = xsecf.Get("hRawYields")

#denominator is a sum over all pT bins
denom = 0.0
for i in range(1, xsec.GetNbinsX()):
    if acc2mu.GetBinContent(i) == 0:
        continue
    denom += (c1*c2*hN2mu.GetBinContent(i)*acc2mu2e.GetBinContent(i)/acc2mu.GetBinContent(i)) 

#multiply by 1e-6 to make the number more meaningful
br = N2mu2e / (denom*10**-6) * B2mu
#just statistical uncertainty for now!!
#unct = N2mu2e**0.5 / (denom*10**-6) 
unct = Err2mu2e / (denom*10**-6) * B2mu

print("%s %s (2mu fit: %s, %s) Blinded BR: (%f +/- %f (stat))e-6"%(selection, "with new weights" if new_wt else "", refMod, refBkg, br, unct)) 
