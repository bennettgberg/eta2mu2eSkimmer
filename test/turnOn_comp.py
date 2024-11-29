import ROOT

#compare the trigger turn-ons for the nominal, up, and down varied toyMC trigger effs

rebin = 5

#use mumu MC instead of signal?
mumu = True
#what % variation to plot?
vary = 5
if mumu:
    distname = "hMmumu"
else:
    distname = "hMmmelel"
freal = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest38110.root"%("mumu" if mumu else "sig"))
hreal = freal.Get(distname)
hreal.SetName("hreal")
hreal.SetTitle("")
hreal.SetLineColor(ROOT.kBlack)
#hreal.SetLineWidth(2)
hreal.Rebin(rebin)
hreal.SetStats(ROOT.kFALSE)
hreal.GetXaxis().SetTitle("m_{2#mu2e} (GeV)")
hreal.GetYaxis().SetTitle("Event / 0.005 GeV") 
fnom = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest38114.root"%("mumu" if mumu else "sig"))
hnom = fnom.Get(distname)
hnom.SetName("hnom")
hnom.SetTitle("")
hnom.Rebin(rebin)
hnom.SetLineWidth(2)
hnom.SetLineColor(ROOT.kBlue)
hnom.SetStats(ROOT.kFALSE)
#fup = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3891.root")
if vary == 10:
    fup = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest38115.root"%("mumu" if mumu else "sig"))
elif vary == 20:
    fup = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest38119.root"%("mumu" if mumu else "sig"))
elif vary == 5:
    fup = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest38117.root"%("mumu" if mumu else "sig"))
hup = fup.Get(distname)
hup.SetName("hup")
hup.Rebin(rebin)
hup.SetLineWidth(2)
hup.SetLineColor(ROOT.kRed)
#fdn = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3892.root")
if vary == 10:
    fdn = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest38116.root"%("mumu" if mumu else "sig"))
elif vary == 20:
    fdn = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest38120.root"%("mumu" if mumu else "sig"))
elif vary == 5:
    fdn = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest38118.root"%("mumu" if mumu else "sig"))
hdn = fdn.Get(distname)
hdn.SetName("hdn")
hdn.Rebin(rebin)
hdn.SetLineWidth(2)
hdn.SetLineColor(ROOT.kGreen)
hdn.SetStats(ROOT.kFALSE)
hdn.SetTitle("")

hdn.GetXaxis().SetTitle("m_{2#mu2e} (GeV)")
hdn.GetYaxis().SetTitle("Event / 0.005 GeV") 

hdn.Draw("hist") 
hreal.Draw("PE same")
hnom.Draw("hist same")
hup.Draw("hist same")

leg = ROOT.TLegend()
leg.AddEntry(hreal, "Real MC")
leg.AddEntry(hnom, "Nominal")
leg.AddEntry(hup, "Turn-on threshold +%d%%"%(vary)) 
leg.AddEntry(hdn, "Turn-on threshold -%d%%"%(vary)) 
leg.Draw("same")

input("kj")
