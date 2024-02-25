import ROOT

fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3818_ALL.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3852.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3853.root")
fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3854.root")

hdata = fdata.Get("hMmumu")
hdata.SetName("hdata")
hmc = fmc.Get("hMmumu")
hmc.SetName("hmc")

hdata.SetMarkerStyle(20)
hdata.SetMarkerSize(0.5)
hdata.SetLineColor(ROOT.kBlack)
hdata.GetXaxis().SetTitle("m_{#mu#mu} (GeV)")
hdata.GetYaxis().SetTitle("Events / MeV")
hdata.SetTitle("Dimuon data with Overlaid #eta#rightarrow#mu#mu MC")
hdata.SetStats(0)
hdata.Draw("PE")

#add a constant 670K to each MC hist bin to align it with the data histogram
for i in range(hmc.GetNbinsX()):
    hmc.SetBinContent(i, hmc.GetBinContent(i) + 665000)

hmc.SetLineWidth(2)
hmc.Draw("same")

leg = ROOT.TLegend()
leg.AddEntry(hdata, "#mu#mu data")
leg.AddEntry(hmc, "#eta#rightarrow#mu#mu MC + 665K")
leg.Draw("hist same")

input("djhf")
