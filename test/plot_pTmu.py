import ROOT

inc_el = True

fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3876.root")
fref = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3875.root")

c = ROOT.TCanvas()

dist = "hpTMummelel"
hsig = fsig.Get(dist)
hsig.SetName("hsig")
hsig.SetTitle("")
hsig.Rebin(5)
hsig.SetLineColor(ROOT.kBlack)
hsig.SetLineWidth(2)
hsig.GetXaxis().SetTitle("Reconstructed p_{T} (GeV)")
hsig.GetYaxis().SetTitle("Weighted events / GeV")
hsig.SetStats(ROOT.kFALSE)

dist = "hpTMumumu"
href = fref.Get(dist)
href.SetName("href")
href.SetTitle("")
href.Rebin(5)
href.SetLineColor(ROOT.kRed)
href.SetLineWidth(2)

#hsig.Draw("hist")
#href.Draw("hist same")
hsig.DrawNormalized("hist")
href.DrawNormalized("hist same")

leg = ROOT.TLegend()
leg.AddEntry(hsig, "#eta #rightarrow #mu#muee MC simulation")
leg.AddEntry(href, "#eta #rightarrow #mu#mu MC simulation")
if inc_el:
    hel = fsig.Get("hpTElmmelel")
    hel.Rebin(5)
    hel.SetLineColor(ROOT.kGreen)
    hel.SetLineWidth(2)
    #hel.Draw("hist same")
    hel.DrawNormalized("hist same")
    leg.AddEntry(hel, "electrons in #eta#rightarrow#mu#mu ee MC") 

leg.Draw("same")

c.SetLogy()

input("kjkjkjkjkj")
