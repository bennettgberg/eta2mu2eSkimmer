import ROOT

sig = False
#Loosest cuts
#fNoCuts = ROOT.TFile.Open("bparking_bkgMCtest360.root")
#fNoCuts = ROOT.TFile.Open("bparking_bkgMCtest3811.root")
#fNoCuts = ROOT.TFile.Open("bparking_bkgMCtest3818.root")
if sig:
    fNoCuts = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3827.root")
else:
    fNoCuts = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3827.root")
hNoCuts = fNoCuts.Get("hMmmelel")
hNoCuts.SetName("hNoCuts")
if sig:
    hNoCuts.SetTitle("Signal Distribution Comparison")
else:
    hNoCuts.SetTitle("Resonant Background Distribution Comparison")

#Strictest cuts
#fIdCuts = ROOT.TFile.Open("bparking_bkgMCtest389.root")
#fIdCuts = ROOT.TFile.Open("bparking_bkgMCtest3819.root")
if sig:
    fIdCuts = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3826.root")
else:
    fIdCuts = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3826.root")
hIdCuts = fIdCuts.Get("hMmmelel")
hIdCuts.SetName("hIdCuts")
hIdCuts.SetTitle("hIdCuts")

#In-between cuts
#fBasicCuts = ROOT.TFile.Open("bparking_bkgMCtest3812.root")
if sig:
    #fBasicCuts = ROOT.TFile.Open("bparking_sigMCtest3828.root")
    fBasicCuts = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3830.root")
else:
    #fBasicCuts = ROOT.TFile.Open("bparking_bkgMCtest3820.root")
    #fBasicCuts = ROOT.TFile.Open("bparking_bkgMCtest3828.root")
    fBasicCuts = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3830.root")
hBasicCuts = fBasicCuts.Get("hMmmelel")
hBasicCuts.SetName("hBasicCuts")

hNoCuts.SetLineWidth(2)
hIdCuts.SetLineWidth(2)
hBasicCuts.SetLineWidth(2)

hNoCuts.SetLineColor(ROOT.kBlack)
hIdCuts.SetLineColor(ROOT.kGreen)
#hBasicCuts.SetLineColor(ROOT.kRed)
hBasicCuts.SetLineColor(ROOT.kMagenta)

if sig:
    rebin = 5
    binsize = ".005"
else:
    rebin = 10
    binsize = ".01"
hNoCuts.Rebin(rebin)
hIdCuts.Rebin(rebin)
hBasicCuts.Rebin(rebin)

hNoCuts.GetXaxis().SetTitle("m_{#mu#muee} (GeV)")
hNoCuts.GetYaxis().SetTitle("Events / %s GeV"%binsize)

c = ROOT.TCanvas()
c.cd()
hNoCuts.SetStats(ROOT.kFALSE)
if sig:
    hNoCuts.Draw("hist")
    hBasicCuts.Draw("hist same")
    hIdCuts.Draw("hist same")
else:
    hNoCuts.Draw()
    hBasicCuts.Draw("same")
    hIdCuts.Draw("same")

leg = ROOT.TLegend()
#leg.AddEntry(hNoCuts, "No selections")
leg.AddEntry(hNoCuts, "vProb > .1 ; nMissingHits <= 3")
#leg.AddEntry(hBasicCuts, "vProb > .5 ; nMissingHits <= 3")
leg.AddEntry(hBasicCuts, "vProb > .1 ; nMissingHits == 0")
leg.AddEntry(hIdCuts, "vProb > .5 ; nMissingHits == 0")
leg.Draw("same")

c.Modified()
c.Update()

int0 = hNoCuts.Integral()
int1 = hBasicCuts.Integral()
int2 = hIdCuts.Integral()

print("Integrals:")
print("vProb>.1,nMissing<=3: %f"%(int0))
#print("vProb>.5,nMissing<=3: %f"%(int1))
print("vProb>.1,nMissing==0: %f"%(int1))
print("vProb>.5,nMissing==0: %f"%(int2))
input("h")
