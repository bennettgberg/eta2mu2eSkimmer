#root [2] hpTGenAll->Draw()
#Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
#root [3] hpTGenReco->SetLineColor(kRed)
#root [4] hpTGenReco->Draw("same")
#root [5] hpTGenReco->Draw()
#Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
#root [6] hpTGenAll->Draw("same")
#root [7] hpTGenRecoe->SetLineColor(kGreen)
#root [8] hpTGenRecoe->Draw("same")
#root [9] hpTGenAll->Rebin(5)
#(TH1 *) 0x55fcb3e0c960
#root [10] hpTGenReco->Rebin(5)
#(TH1 *) 0x55fcb45460e0
#root [11] hpTGenRecoe->Rebin(5)
#(TH1 *) 0x55fcb47018a0
#root [12] hpTGenReco->Divide(hpTGenAll)
#(bool) true
#root [13] hpTGenRecoe->Divide(hpTGenAll)
#(bool) true
#root [14] hpTGenReco->Draw()
#Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
#root [15] hpTGenRecoe->Draw("same")
#root [16] .q
import ROOT
year = 2022
#specify tight to use tight ID result (WP80) instead of nominal (WP90)
tight = False
#what factor to rebin by
rebin = 5 #20
sig = True
mumu = False
eff = True 
#trig = True
gsf = True
#type of vertex to use
vtype = "mmelel"
#use pseudorapidity eta instead of pT?
eta = False
if mumu: vtype = "mumu"
#fname = "bparking_%sMCtest15.root"%("sig" if sig else "bkg") 
#fname = "bparking_%sMCtest23.root"%("sig" if sig else "bkg") 
#fname = "bparking_%sMCtest22.root"%("central") 
#fname = "bparking_%sMCtest337.root"%("sig" if sig else "bkg") 
#fname = "bparking_%sMCtest337.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest3311.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest371.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest373.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest374.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest377.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest387.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest3819.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest3826.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest3827.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#nominal
#fname = "bparking_%sMCtest3837.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
##fname = "bparking_%sMCtest3840.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#no elID
#fname = "bparking_%sMCtest3850.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#updated xsec weights, nd stuff
#fname = "bparking_%sMCtest3853.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#fname = "bparking_%sMCtest%d.root"%("sig" if sig else ("mumu" if mumu else "bkg"), 3865 if mumu else 3868) 
#fname = "rootfiles/bparking_%s%sMCtest%d.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), 38110) 
#fname = "rootfiles/bparking_%s%sMCtest%d.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), 38125) 
#fname = "rootfiles/bparking_%s%sMCtest%d.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), 470) 
#fname = "rootfiles/bparking_%s%sMCtest%s.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), "38126") 
#fname = "rootfiles/bparking_%s%sMCtest%s.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), "38128") 
#fname = "rootfiles/bparking_%s%sMCtest%s.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), "38130p474") 
#fname = "rootfiles/bparking_%s%sMCtest%s.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), "38133") 
#fname = "rootfiles/bparking_%s%sMCtest%s.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), "38134") 
#fname = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%s%sMCtest%s.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), "38137p478" if mumu else ("38164p4733" if sig else "38165") )
fname = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%s%sMCtest%s.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), "38192p4755" if mumu else ("38193p4756" if sig else "38194") )
#if year == 2023:
#    #fname = "bparking_2023%sMCtest400.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
#    fname = "bparking_2023%sMCtest430.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
##loose
#fname = "bparking_%sMCtest3838.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
##tight
#fname = "bparking_%sMCtest3839.root"%("sig" if sig else ("mumu" if mumu else "bkg")) 
print("Opening " + fname)
f = ROOT.TFile.Open(fname)
hpTGenAll = f.Get("h%sGenAll"%("Eta" if eta else "pT"))
if gsf:
    #hpTGenReco = f.Get("hpTGenRecoe") 
    hpTGenReco = f.Get("h%sGenReco%s"%("Eta" if eta else "pT", vtype))
    #for now, separate file for Acceptance (for memory resource reasons)
    if tight: 
        fAcc = ROOT.TFile.Open("bparking_%s%sMCtest%d.root"%(str(year) if year != 2022 else "", "sig" if sig else ("mumu" if mumu else "bkg"), (3865 if mumu else 3866) if year == 2022 else 430)) 
    else:
        #fAcc = ROOT.TFile.Open("rootfiles/bparking_%sMCtest%d.root"%("sig" if sig else ("mumu" if mumu else "bkg"), 3865 if mumu else 3871)) 
        fAcc = f
    #hpTGenAcc = f.Get("h%sGenAcc%s"%("Eta" if eta else "pT", vtype))
    hpTGenAcc = fAcc.Get("h%sGenAccdR%s"%("Eta" if eta else "pT", vtype))
else:
    hpTGenReco = f.Get("h%sGenReco"%("Eta" if eta else "pT")) 
    hpTGenAcc = f.Get("h%sGenAcc"%("Eta" if eta else "pT"))
hpTGenTrig = f.Get("h%sGenTrig"%("Eta" if eta else "pT")) 

hpTGenAll.Rebin(rebin)
hpTGenReco.Rebin(rebin)
hpTGenTrig.Rebin(rebin)
hpTGenAcc.Rebin(rebin)

hpTGenReco.SetLineColor(ROOT.kRed)
hpTGenTrig.SetLineColor(ROOT.kGreen+1)
hpTGenAcc.SetLineColor(ROOT.kBlue)
if eff:
    hpTGenReco.Sumw2()
    hpTGenTrig.Sumw2()
    hpTGenAcc.Sumw2()
hpTGenReco.SetLineWidth(2)
hpTGenTrig.SetLineWidth(2)
hpTGenAcc.SetLineWidth(2)

if eff:
    hpTGenTrig.Divide(hpTGenAll) 
    hpTGenReco.Divide(hpTGenAll)
    hpTGenAcc.Divide(hpTGenAll)
else:
    hpTGenAll.SetLineColor(ROOT.kBlack)
    hpTGenAll.SetLineWidth(2)
binsize = (hpTGenAll.GetXaxis().GetXmax() - hpTGenAll.GetXaxis().GetXmin()) / hpTGenAll.GetXaxis().GetNbins() 

if eff:
    hpTGenTrig.GetXaxis().SetTitle("#eta meson gen %s"%("pseudorapidity" if eta else "p_{T} (GeV)")) 
    hpTGenTrig.GetYaxis().SetTitle("Efficiency")
else:
    hpTGenAll.GetXaxis().SetTitle("#eta meson gen %s"%("pseudorapidity" if eta else "p_{T} (GeV)")) 
    if eta:
        hpTGenAll.GetYaxis().SetTitle("Events") 
    else:
        hpTGenAll.GetYaxis().SetTitle("Events / %.2f GeV"%binsize) 

if eff:
    hpTGenTrig.SetStats(ROOT.kFALSE)
    hpTGenTrig.SetTitle("")
    hpTGenTrig.Draw()
    hpTGenReco.Draw("same")
    hpTGenAcc.Draw("same")
else:
    hpTGenAll.SetStats(ROOT.kFALSE)
    hpTGenAll.SetTitle("")
    hpTGenAll.Draw("hist")
    hpTGenTrig.Draw("hist same")
    hpTGenReco.Draw("hist same")
    hpTGenAcc.Draw("hist same")
leg = ROOT.TLegend()
if not eff:
    leg.AddEntry(hpTGenAll, "All gen #eta mesons")
    leg.AddEntry(hpTGenTrig, "Passing trigger")
    leg.AddEntry(hpTGenReco, "Passing reconstruction")
    leg.AddEntry(hpTGenAcc, "Accepted")
else:
    leg.AddEntry(hpTGenTrig, "Trigger Efficiency")
    leg.AddEntry(hpTGenReco, "Reconstruction Efficiency")
    leg.AddEntry(hpTGenAcc, "Acceptance")
leg.Draw("same")

if eff and not eta:
    #calculate, print the overall weighted uncertainties!
    totW = 0.0
    trgW = 0.0
    recW = 0.0
    accW = 0.0
    #fitted parameter for the exponential in the xsec fit
    exp_param = 5.1121328
    for i in range(11, 66):
        totW += 1.0*1e11 / hpTGenAll.GetBinCenter(i)**exp_param
        trgW += 1.0*1e11*hpTGenTrig.GetBinContent(i) / hpTGenTrig.GetBinCenter(i)**exp_param
        recW += 1.0*1e11*hpTGenReco.GetBinContent(i) / hpTGenReco.GetBinCenter(i)**exp_param
        accW += 1.0*1e11*hpTGenAcc.GetBinContent(i) / hpTGenAcc.GetBinCenter(i)**exp_param
        #print("i=%d Running totW: %f Running thisW: %f"%(i, totW, thisW)) 
    print("Overall weighted efficiencies (from 10-65 GeV):") 
    print("Trigger: %f%%"%(trgW/totW*100.0)) 
    print("Reco: %f%%"%(recW/totW*100.0)) 
    print("Acceptance: %f%%"%(accW/totW*100.0)) 
#raw_input("h?")
input("h?")
