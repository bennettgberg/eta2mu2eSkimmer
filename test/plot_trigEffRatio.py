import ROOT

#plot the ratio of the efficiencies?
ratio = True 
#plot the efficiencies, or just the raw distributions?
eff = True
sig = True 
data = True  
if ratio and not data:
    print("Error: can't do ratio without data.")
    exit()
#-1 for EC only, 1 for barrel only, 0 for inclusive
etaregion = 1
#factor by which to rebin
rebin = 20
#lowest value of x to include in the hist
xlow = 1.0
#highest value of x to go up to with the hists
xhigh = 17.0
if xlow > 0:
    hL1TdataS = ROOT.TH1F("hL1TdataS", "", int((xhigh-xlow)/2), xlow, xhigh) 
    hL1TsigS = ROOT.TH1F("hL1TdataS", "", int((xhigh-xlow)/2), xlow, xhigh) 
    hHLTdataS = ROOT.TH1F("hHLTdataS", "", int((xhigh-xlow)/2), xlow, xhigh) 
    hHLTsigS = ROOT.TH1F("hHLTdataS", "", int((xhigh-xlow)/2), xlow, xhigh) 
    hL1TdataS.Sumw2()
    hL1TsigS.Sumw2()
    hHLTdataS.Sumw2()
    hHLTsigS.Sumw2()
if data:
    #fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest4746_ALL.root")
    fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest4747_ALL.root")
    if etaregion != -1:
        hL1Tdata = fdata.Get("hL1TC")
        hHLTdata = fdata.Get("hHLTC")
        if etaregion == 0:
            hL1Tdataf = fdata.Get("hL1TF")
            hL1Tdata.Add(hL1Tdataf)
            hHLTdataf = fdata.Get("hHLTF")
            hHLTdata.Add(hHLTdataf)
    else:
        hL1Tdata = fdata.Get("hL1TF")
        hHLTdata = fdata.Get("hHLTF")
    hL1Tdata.SetName("hL1Tdata")
    #the L1 path is prescaled by 1000 in data!!!
    hL1Tdata.Scale(1000)
    if xlow > 0:
        #fill the new histograms
        for ibin in range(hL1Tdata.GetNbinsX()):
            this_bin = hL1TdataS.FindBin(hL1Tdata.GetBinCenter(ibin))
            hL1TdataS.SetBinContent(this_bin, hL1TdataS.GetBinContent(this_bin) + hL1Tdata.GetBinContent(ibin))
            hL1TdataS.SetBinError(this_bin, (hL1TdataS.GetBinError(this_bin)**2 + hL1Tdata.GetBinError(ibin)**2)**0.5)
        for ibin in range(hHLTdata.GetNbinsX()):
            this_bin = hHLTdataS.FindBin(hHLTdata.GetBinCenter(ibin))
            hHLTdataS.SetBinContent(this_bin, hHLTdataS.GetBinContent(this_bin) + hHLTdata.GetBinContent(ibin))
            hHLTdataS.SetBinError(this_bin, (hHLTdataS.GetBinError(this_bin)**2 + hHLTdata.GetBinError(ibin)**2)**0.5)
        hL1Tdata = hL1TdataS
        hL1Tdata.SetName("hL1Tdata")
        hHLTdata.SetName("hHLTdata")
        hHLTdata = hHLTdataS
    elif rebin > 1:
        hL1Tdata.Rebin(rebin)
        hHLTdata.Rebin(rebin)
    hHLTdata.GetXaxis().SetTitle("sub-leading muon p_{T} (GeV)") 
    hHLTdata.SetStats(ROOT.kFALSE)
    hHLTdata.SetLineWidth(2)
    hHLTdata.SetLineColor(ROOT.kBlack)
    if etaregion == 0:
        hHLTdata.SetTitle("Inclusive pseudorapidity")
    elif etaregion == -1:
        hHLTdata.SetTitle("|#eta| > 1.2")
    elif etaregion == 1:
        hHLTdata.SetTitle("|#eta| < 1.2")
    if eff:
        hHLTdata.GetYaxis().SetTitle("Trigger efficiency")
        hHLTdata.Sumw2()
        hHLTdata.Divide(hL1Tdata) 
if sig:
    #fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38184p4746.root")
    #fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38185p4747.root")
    fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest4748.root")
    if etaregion != -1:
        hL1Tsig = fsig.Get("hL1TC")
        hHLTsig = fsig.Get("hHLTC")
        if etaregion == 0:
            hL1Tsigf = fsig.Get("hL1TF")
            hL1Tsig.Add(hL1Tsigf)
            hHLTsigf = fsig.Get("hHLTF")
            hHLTsig.Add(hHLTsigf)
    else:
        hL1Tsig = fsig.Get("hL1TF")
        hHLTsig = fsig.Get("hHLTF")
    hL1Tsig.SetName("hL1Tsig")
    if xlow > 0:
        for ibin in range(hL1Tsig.GetNbinsX()):
            this_bin = hL1TsigS.FindBin(hL1Tsig.GetBinCenter(ibin))
            hL1TsigS.SetBinContent(this_bin, hL1TsigS.GetBinContent(this_bin) + hL1Tsig.GetBinContent(ibin))
            hL1TsigS.SetBinError(this_bin, (hL1TsigS.GetBinError(this_bin)**2 + hL1Tsig.GetBinError(ibin)**2)**0.5)
        for ibin in range(hHLTsig.GetNbinsX()):
            this_bin = hHLTsigS.FindBin(hHLTsig.GetBinCenter(ibin))
            hHLTsigS.SetBinContent(this_bin, hHLTsigS.GetBinContent(this_bin) + hHLTsig.GetBinContent(ibin))
            hHLTsigS.SetBinError(this_bin, (hHLTsigS.GetBinError(this_bin)**2 + hHLTsig.GetBinError(ibin)**2)**0.5)
        hL1Tsig = hL1TsigS
        hHLTsig = hHLTsigS
        hL1Tsig.SetName("hL1Tsig")
        hHLTsig.SetName("hHLTsig")
    elif rebin > 1:
        hL1Tsig.Rebin(rebin)
        hHLTsig.Rebin(rebin)
    hHLTsig.GetXaxis().SetTitle("sub-leading muon p_{T} (GeV)") 
    hHLTsig.SetStats(ROOT.kFALSE)
    hHLTsig.SetLineWidth(2)
    hHLTsig.SetLineColor(ROOT.kRed)
    if etaregion == 0:
        hHLTsig.SetTitle("Inclusive pseudorapidity")
    elif etaregion == -1:
        hHLTsig.SetTitle("|#eta| > 1.2")
    elif etaregion == 1:
        hHLTsig.SetTitle("|#eta| < 1.2")
    if eff:
        hHLTsig.GetYaxis().SetTitle("Trigger efficiency")
        hHLTsig.Sumw2()
        hHLTsig.Divide(hL1Tsig)
    
if ratio:
    c = ROOT.TCanvas()
    c.Divide(1,2)
    top_pad = c.cd(1)
    top_pad.SetPad(0, 0.3, 1, 1)  # Set top pad to cover top 70% of the canvas
    top_pad.SetBottomMargin(0)
if sig:
    hHLTsig.Draw() 
    if data:
        hHLTdata.Draw("same")
        leg = ROOT.TLegend()
        leg.AddEntry(hHLTsig, "Reference MC")
        leg.AddEntry(hHLTdata, "Data")
        leg.Draw("same")
elif data:
    hHLTdata.Draw()
if ratio:
    rat = hHLTdata.Clone("Ratio")
    bottom_pad = c.cd(2)
    bottom_pad.SetPad(0, 0, 1, 0.3)  # Set bottom pad to cover bottom 30% of the canvas
    bottom_pad.SetTopMargin(0)       # Remove top margin for bottom pad
    bottom_pad.SetBottomMargin(0.3)
    rat.GetYaxis().SetTitle("Data / MC ratio")
    rat.SetLineWidth(2)
    rat.SetTitle("")
    rat.GetXaxis().SetTitleSize(0.1)
    rat.GetYaxis().SetTitleSize(0.08)
    rat.GetXaxis().SetLabelSize(0.08)
    rat.GetYaxis().SetLabelSize(0.08)
    rat.GetYaxis().SetTitleOffset(0.5)
    rat.SetMinimum(0.75)
    rat.SetMaximum(1.00)
    rat.Draw()
input("Write out the new trigger corrections? (ratio mode only)")

if ratio:
    #write out the new trigger corrections
    hname = "trigger_corrections%s"%("C" if etaregion == 1 else "F") 
    outfn = "%s.root"%(hname)
    outf = ROOT.TFile.Open(outfn, "recreate")
    rat.SetName(hname)
    rat.Write()
    outf.Close()
    print("%s written."%(outfn))
