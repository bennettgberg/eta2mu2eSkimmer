import ROOT

eff = True
inc_el = True

#fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3878.root") 
fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3881.root") 

hpT0 = fsig.Get("hGenMupT0")
hpT1 = fsig.Get("hGenMupT1")
hpT0.Add(hpT1)

if inc_el:
    hpTElAll = fsig.Get("hGenElpTAll")
    hpTElAll.Rebin(100)
    hpTElRec = fsig.Get("hGenElpTRec")
    hpTElRec.Rebin(100)

hpT0.Rebin(100)

hpTRec = fsig.Get("hGenMupTRec")
hpTRec.Rebin(100)
hpTRec.Sumw2()

if eff:
    hpTRec.Divide(hpT0)
    if inc_el:  
        hpTElRec.Sumw2()
        hpTElRec.Divide(hpTElAll)

if eff:
    hpTRec.SetMarkerStyle(8)
    hpTRec.SetMarkerColor(ROOT.kRed)
hpT0.SetStats(ROOT.kFALSE)
hpTRec.SetStats(ROOT.kFALSE)
hpT0.SetTitle("")
hpTRec.SetTitle("") 
if inc_el and eff:
    hpTElRec.SetMarkerStyle(8)
    hpTElRec.SetMarkerColor(ROOT.kGreen)

if eff:
    if inc_el:
        hpTRec.GetXaxis().SetTitle("Lepton p_{T} (GeV)")
    else:
        hpTRec.GetXaxis().SetTitle("Muon p_{T} (GeV)")
    hpTRec.GetYaxis().SetTitle("Reconstruction efficiency")
else:
    hpT0.SetLineColor(ROOT.kBlack)
    hpTRec.SetLineColor(ROOT.kRed)
    if inc_el:
        hpT0.GetXaxis().SetTitle("Lepton p_{T} (GeV)")
    else:
        hpT0.GetXaxis().SetTitle("Muon p_{T} (GeV)")
    hpT0.GetYaxis().SetTitle("Events / GeV")
if inc_el or not eff:
    leg = ROOT.TLegend()
    if not eff:
        leg.AddEntry(hpT0, "All gen muons")
    leg.AddEntry(hpTRec, "Reconstructed gen-matched muons") 
    if inc_el:
        if not eff:
            hpTElAll.SetLineColor(ROOT.kBlue)
            hpTElRec.SetLineColor(ROOT.kGreen)
            leg.AddEntry(hpTElAll, "All gen electrons")
        leg.AddEntry(hpTElRec, "Reconstructed gen-matched electrons")

if eff:
    hpTRec.Draw("PE")
    if inc_el:
        hpTElRec.Draw("PE same")
else:
    hpT0.Draw("hist")
    hpTRec.Draw("hist same")
    if inc_el:
        hpTElAll.Draw("hist same")
        hpTElRec.Draw("hist same")
if inc_el or not eff:
    leg.Draw("same") 
input("h....") 
