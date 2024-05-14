import ROOT

#compare the trigger turn-ons for the nominal, up, and down varied toyMC trigger effs

rebin = 5

#use mumu MC instead of signal?
mumu = False
fnom = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3890.root")
hnom = fnom.Get("hMmmelel")
hnom.SetName("hnom")
hnom.SetTitle("")
hnom.Rebin(rebin)
hnom.SetLineWidth(2)
hnom.SetLineColor(ROOT.kBlack)
hnom.SetStats(ROOT.kFALSE)
hnom.GetXaxis().SetTitle("m_{2#mu2e} (GeV)")
hnom.GetYaxis().SetTitle("Event / 0.005 GeV") 
fup = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3891.root")
hup = fup.Get("hMmmelel")
hup.SetName("hup")
hup.Rebin(rebin)
hup.SetLineWidth(2)
hup.SetLineColor(ROOT.kRed)
fdn = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3892.root")
hdn = fdn.Get("hMmmelel")
hdn.SetName("hdn")
hdn.Rebin(rebin)
hdn.SetLineWidth(2)
hdn.SetLineColor(ROOT.kGreen)

hnom.Draw("hist")
hup.Draw("hist same")
hdn.Draw("hist same") 

leg = ROOT.TLegend()
leg.AddEntry(hnom, "Nominal")
leg.AddEntry(hup, "Turn-on threshold +10%%") 
leg.AddEntry(hdn, "Turn-on threshold -10%%") 
leg.Draw("same")

input("kj")
