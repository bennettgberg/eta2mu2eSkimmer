import ROOT

eff = True
inc_el = True
inc_mu = False #True

#rebin factor
rebin = 40 #100

#do a fit to the electron eff distribution?
do_fit = True

if not inc_el and not inc_mu:
    print("what the h*ck no muons or electrons bruv???")
    exit()

#fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3878.root") 
#fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3881.root") 
fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3884.root") 

if inc_mu:
    hpT0 = fsig.Get("hGenMupT0")
    hpT1 = fsig.Get("hGenMupT1")
    hpT0.Add(hpT1)
    hpT0.Rebin(rebin)

if inc_el:
    hpTElAll = fsig.Get("hGenElpTAll")
    hpTElAll.Rebin(rebin)
    hpTElRec = fsig.Get("hGenElpTRec")
    hpTElRec.Rebin(rebin)

if inc_mu:
    hpTRec = fsig.Get("hGenMupTRec")
    hpTRec.Rebin(rebin)
    hpTRec.Sumw2()

if eff:
    if inc_mu:
        hpTRec.Divide(hpT0)
    if inc_el:  
        hpTElRec.Sumw2()
        hpTElRec.Divide(hpTElAll)

if eff:
    if inc_mu:
        hpTRec.SetMarkerStyle(8)
        hpTRec.SetMarkerColor(ROOT.kRed)
if inc_mu:
    hpT0.SetStats(ROOT.kFALSE)
    hpTRec.SetStats(ROOT.kFALSE)
    hpT0.SetTitle("")
    hpTRec.SetTitle("") 
if inc_el and eff:
    hpTElRec.SetTitle("")
    hpTElRec.SetStats(ROOT.kFALSE)
    hpTElRec.SetMarkerStyle(8)
    hpTElRec.SetMarkerColor(ROOT.kGreen)

if eff:
    if inc_el and not inc_mu:
        hpTElRec.GetXaxis().SetTitle("Lepton p_{T} (GeV)")
        hpTElRec.GetYaxis().SetTitle("Reconstruction efficiency")
    elif inc_mu and inc_el:
        hpTRec.GetXaxis().SetTitle("gen lepton p_{T} (GeV)")
        hpTRec.GetYaxis().SetTitle("Reconstruction efficiency")
    else:
        hpTRec.GetXaxis().SetTitle("Muon p_{T} (GeV)")
        hpTRec.GetYaxis().SetTitle("Reconstruction efficiency")
else:
    hpT0.SetLineColor(ROOT.kBlack)
    hpTRec.SetLineColor(ROOT.kRed)
    if inc_el and inc_mu:
        hpT0.GetXaxis().SetTitle("gen lepton p_{T} (GeV)") 
    elif inc_el:
        hpTEl.GetXaxis().SetTitle("gen electron p_{T} (GeV)") 
    else:
        hpT0.GetXaxis().SetTitle("Muon p_{T} (GeV)")
    hpT0.GetYaxis().SetTitle("Events / GeV")
if inc_el or not eff:
    leg = ROOT.TLegend()
    if not eff:
        leg.AddEntry(hpT0, "All gen muons")
    if inc_mu:
        if eff:
            leg.AddEntry(hpTRec, "Muon reconstruction efficiency") 
        else:
            leg.AddEntry(hpTRec, "Reconstructed gen-matched muons") 
    if inc_el:
        if not eff:
            hpTElAll.SetLineColor(ROOT.kBlue)
            hpTElRec.SetLineColor(ROOT.kGreen)
            leg.AddEntry(hpTElAll, "All gen electrons")
            leg.AddEntry(hpTElRec, "Reconstructed gen-matched electrons")
        else:
            leg.AddEntry(hpTElRec, "Electron reconstruction efficiency")

if eff:
    if inc_mu:
        hpTRec.Draw("PE")
    if inc_el and not inc_mu:
        hpTElRec.Draw("PE")
    elif inc_el:
        hpTElRec.Draw("PE same")
else:
    if inc_mu:
        hpT0.Draw("hist")
        hpTRec.Draw("hist same")
    if inc_el:
        if inc_mu:
            hpTElAll.Draw("hist")
        else:
            hpTElAll.Draw("hist same")
        hpTElRec.Draw("hist same")
if (inc_el and inc_mu) or not eff:
    leg.Draw("same") 

if do_fit:
    #fit the electron eff distribution from 0 to 40 GeV
    #tf = ROOT.TF1("tf1", "([0] / ([1]*x + TMath::Exp(-(x - [2])/[3])))-[4]/x", 0.5, 40)
    tf = ROOT.TF1("tf1", "([0] / ([1]*x + TMath::Exp(-(x - [2])/[3])))-[4]/x+[5]*(x-[6])*(x-[6])", 0.5, 27)
    tf.SetParameter(0, 3.00)
    tf.SetParameter(1, 1.00)
    tf.SetParameter(2, 13.0)
    tf.SetParameter(3, 4.1)
    tf.SetParameter(4, 0.022)
    #tf.SetParameter(5, 0.8)
    #tf.SetParameter(6, 4.4)
    #tf.SetParameter(7, 0.8)
    ##tf = ROOT.TF1("tf1", "TMath::Landau(x, [0], [1], [2])", 10, 40)
    hpTElRec.Fit(tf)
    #value at 27 GeV, where the plateau begins
    v27 = tf.Eval(27.0)
    print("v27: " + str(v27)) 
    tline = ROOT.TLine(27, v27, 45, v27)
    tline.SetLineColor(ROOT.kBlue)
    tline.SetLineWidth(2)
    tline.Draw("same")
input("h....") 
