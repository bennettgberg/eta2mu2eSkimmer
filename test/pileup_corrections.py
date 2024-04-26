import ROOT

#use new weights instead of old?
new_wt = True

year = 2022

if year == 2022:
    fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3854_ALL.root")
    if not new_wt:
        fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3854.root")
        fbkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3854.root")
        fref = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3854.root")
    else:
        fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3877.root")
        fbkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3877.root")
        fref = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3877.root")
else:
    fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2023/ultraskimmed/bparking_2023datatest430_ALL.root")
    fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2023/ultraskimmed/bparking_2023sigMCtest430.root")
    #fbkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2023/ultraskimmed/bparking_2023bkgMCtest430.root")
    fref = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2023/ultraskimmed/bparking_2023mumuMCtest430.root")

hdata = fdata.Get("hnPV")
hdata.SetName("hnPVData")

hsig = fsig.Get("hnPV")
hsig.SetName("hnPVSigMC")

if year == 2022:
    hbkg = fbkg.Get("hnPV")
    hbkg.SetName("hnPVBkgMC")

href = fref.Get("hnPV")
href.SetName("hnPVRefMC")

hdata.SetMarkerStyle(20)
hdata.SetMarkerSize(0.6)
hdata.SetMarkerColor(ROOT.kBlack)
hdata.SetLineColor(ROOT.kBlack)
hdata.SetStats(ROOT.kFALSE)
hdata.GetXaxis().SetTitle("Number of Primary Vertices")
hdata.GetYaxis().SetTitle("Normalized Events") 
hdataNorm = hdata.DrawNormalized("PE")

hsig.SetLineWidth(2)
hsig.SetLineColor(ROOT.kOrange)
hsigNorm = hsig.DrawNormalized("hist same")
hsigNorm.Print()

if year == 2022:
    hbkg.SetLineWidth(2)
    hbkg.SetLineColor(ROOT.kMagenta)
    hbkgNorm = hbkg.DrawNormalized("hist same")
    hbkgNorm.Print()

href.SetLineWidth(2)
href.SetLineColor(ROOT.kCyan)
hrefNorm = href.DrawNormalized("hist same")
hrefNorm.Print()

leg = ROOT.TLegend()
leg.AddEntry(hdata, "Data")
leg.AddEntry(hsig, "Signal MC")
if year == 2022:
    leg.AddEntry(hbkg, "Resonant Background MC")
leg.AddEntry(href, "Reference Channel MC")
leg.Draw("same")

c2 = ROOT.TCanvas()
#the corrections are just the ratios of the histograms
hcorSig = hdataNorm.Clone()
hcorSig.SetName("hSigCorr%d"%year)
hcorSig.Divide(hsigNorm)
hcorSig.SetLineColor(ROOT.kOrange)
hcorSig.SetLineWidth(2)

if year == 2022:
    hcorBkg = hdataNorm.Clone()
    hcorBkg.SetName("hBkgCorr%d"%year)
    hcorBkg.Divide(hbkgNorm)
    hcorBkg.SetLineColor(ROOT.kMagenta)
    hcorBkg.SetLineWidth(2)

hcorRef = hdataNorm.Clone()
hcorRef.SetName("hRefCorr%d"%year)
hcorRef.Divide(hrefNorm)
hcorRef.SetLineColor(ROOT.kCyan)
hcorRef.SetLineWidth(2)
fname = "pileup_corrections%s_%d.root"%("_newWt" if new_wt else "", year)
fout = ROOT.TFile.Open(fname, "recreate")
hcorSig.Write()
if year == 2022:
    hcorBkg.Write()
hcorRef.Write()

c2.SetLogy()
hcorSig.GetYaxis().SetTitle("Correction")
hcorSig.SetStats(ROOT.kFALSE)
hcorSig.Draw("hist")
if year == 2022:
    hcorBkg.Draw("hist same")
hcorRef.Draw("hist same") 
leg2 = ROOT.TLegend()
leg2.AddEntry(hcorSig, "Signal MC")
if year == 2022:
    leg2.AddEntry(hcorBkg, "Background MC")
leg2.AddEntry(hcorRef, "Reference channel MC")
leg2.Draw("same")

input("end")
