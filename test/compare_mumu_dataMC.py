import ROOT

#include the pre/post-EE distinction?
#incEE = True

#compare pt instead of invariant mass?
compt = False
#subtract the sidebands instead of showing everything added together?
subtract = False

#use signal MC instead of mumu
sig = False

#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3818_ALL.root")
#additional cuts, nd stuff
#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3860_ALL.root")
#loose muID req'd
#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3884_ALL.root")
#one trigger path only instead of 6
#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38105_ALL.root")
#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38107_ALL.root")
#4_3 path only
#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38107_ALL.root")
#cutting on mu pT 3,4 GeV
#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38126_ALL.root")
#NOMINAL
#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38130_ALL.root")
fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest4733_ALL.root")
#TIGHTTTTTTTTTTTT
#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38139_ALL.root")
#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38130_ALL.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3854.root")
#additional cuts, nd stuff (no PU corrections)
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3867.root")
#with PU corrections
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3866.root")
#new evt weights (with DG/Cheb4 2mu fits) -- no PU corrections
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3877.root")
#yes PU
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3878.root")
#loose muID req'd, pileup correction, trigger eff correction
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3885.root")
#back to no trigger eff coreq, but yes loose muID
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3886.root")
#now with mumu-specific trig eff corrections
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3887.root")
#trigger corrections only above 5 GeV. -- NOMINAL frfrfr
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3888.root")
#above 4 GeV only instead of 5
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3889.root")
###toy MC files, for uncty
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3892.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3894.root")
#one trigger path only instead of 6
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38106.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38107.root")
#4_3 path only
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38110.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38110.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38126.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38130.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38136.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38137p478.root")
if sig:
    #tightttttttt
    #fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38139p4710.root")
    #nominal
    #fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38145p4715.root")
    fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38164p4733.root")
else:
    #tighttt
    #fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38139p4710.root")
    #nominal
    fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38145p4715.root")
    #toy trig 0
    #fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38149p4719.root")
    #toy trig -5
    #fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38150p4720.root")
    #toy trig +5
    #fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38151p4721.root")
    #toy trig -8
    #fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38152p4722.root")
    #toy trig -12
    #fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38153p4723.root")
#if incEE:
#    fEE = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest477.root")
#toy trigger varied +10%
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38115.root")
#toy trigger varied -10%
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38116.root")
#toy trigger varied +5%
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38117.root")
#toy trigger varied +8%
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38121.root")
#toy trigger varied -8%
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38122.root")
#miny = 665000
#miny = 460000
#miny = 340000
#miny = 337000
miny = 321000
#miny = 600

c1 = ROOT.TCanvas()
if compt:
    if sig:
        hdata = fdata.Get("hpTSigmmelel")
        #sidebands!
        hsides = fdata.Get("hpTSidemmelel")
        #normalize the sideband to correspond to the number of background events fitted in data
        nSide = hsides.GetSumOfWeights()
        nSideData = 19.0
        hsides.Scale(nSideData / nSide)
    else:
        hdata = fdata.Get("hpTSigmumu")
        #sidebands!
        hsides = fdata.Get("hpTSidemumu")
    frbkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38141.root")
    hrbkg = frbkg.Get("hpTModmu2")
    if subtract:
        hdata.Add(hsides, -1)
        hdata.Add(hrbkg, -1)
    else:
        htot = hsides.Clone("htot")
        htot.Sumw2()
        if sig:
            htot.Add(hrbkg)
        htot.SetStats(ROOT.kFALSE)
else:
    hdata = fdata.Get("hMmumu")
hdata.SetName("hdata")
if compt:
    if sig:
        hmc = fmc.Get("hpTSigmmelel")
    else:
        hmc = fmc.Get("hpTModmu2")
    print("hmc: ")
    hmc.Print()
else:
    hmc = fmc.Get("hMmumu")
    #hmc = fmc.Get("hMModmu3")
hmc.SetName("hmc")
#if incEE:
#    if compt:
#        hEE = fEE.Get("hpTSigmumu")
#    else:
#        #hEE = fEE.Get("hMmumu")
#        hEE = fEE.Get("hMModmu2")
#    #hEE.Add(hmc, 28.25/38.01)
#    hEE.Add(hmc)
#    hmc = hEE

hdata.SetMarkerStyle(20)
hdata.SetMarkerSize(0.5)
hdata.SetLineColor(ROOT.kBlack)
if compt:
    if sig:
        rebin = 50
        hmc.Rebin(rebin)
        hsides.Rebin(rebin)
        if not subtract:
            hrbkg.Rebin(rebin)
            #how much to scale the signal MC by (to account for the blinding factor)
            scale = 0.35 
            hmc.Scale(scale)
            htot.Add(hmc)
            htot.GetXaxis().SetTitle("#mu#muee p_{T} (GeV)")
        else:
            hdata.GetXaxis().SetTitle("#mu#muee p_{T} (GeV)")
    else:
        rebin = 2 #5
        hmc.Rebin(rebin)
        hdata.GetXaxis().SetTitle("#mu#mu p_{T} (GeV)")
        if not subtract:
            hsides.Rebin(rebin)
            htot.Add(hmc)
    hdata.Rebin(rebin)
    if subtract:
        hdata.GetYaxis().SetTitle("Normalized Events")
        #hdata.GetYaxis().SetTitle("Events / GeV")
    else:
        htot.GetYaxis().SetTitle("Events / 10 GeV")
else:
    hdata.GetXaxis().SetTitle("m_{#mu#mu} (GeV)")
    hdata.GetYaxis().SetTitle("Events / MeV")
if sig:
    if compt and not subtract:
        htot.SetTitle("#mu#muee data with Overlaid MC")
    else:
        hdata.SetTitle("#mu#muee data with Overlaid  #eta#rightarrow#mu#muee MC")
else:
    hdata.SetTitle("Dimuon data with Overlaid  #eta#rightarrow#mu#mu MC")
hdata.SetStats(0)
#hdata.Draw("PE")

##try fitting the background to a 3rd order polynomial
#hsub = hdata.Clone()
#hsub.SetName("hsub")
#for i in range(hmc.GetNbinsX()):
#    hsub.SetBinContent(i, hsub.GetBinContent(i) - hmc.GetBinContent(i)) 
##tf3 = ROOT.TF1("tf3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", hsub.GetXaxis().GetXmin(), hsub.GetXaxis().GetXmax()) 
#tf4 = ROOT.TF1("tf4", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", hsub.GetXaxis().GetXmin(), hsub.GetXaxis().GetXmax()) 
##fit3 = hsub.Fit(tf3, "S", "", 0.45, 0.65)
#fit4 = hsub.Fit(tf4, "S", "", 0.45, 0.65)
#hsub.Draw()
#input("continue?")
#
##params = fit3.floatParsFinal()
##params = fit3.GetParams()
#params = fit4.GetParams()
if not compt:
    #add a constant 670K to each MC hist bin to align it with the data histogram
    for i in range(hmc.GetNbinsX()):
        xval = hmc.GetBinCenter(i)
        #bkgi = params[0] + params[1]*xval + params[2]*xval*xval + params[3]*xval*xval*xval
        #bkgi = params[0] + params[1]*xval + params[2]*xval*xval + params[3]*xval*xval*xval + params[4]*xval*xval*xval*xval
        hmc.SetBinContent(i, hmc.GetBinContent(i) + miny)
        #hmc.SetBinContent(i, hmc.GetBinContent(i) + bkgi)

hmc.SetLineWidth(2)

leg = ROOT.TLegend()
if compt:
    if sig:
        if subtract:
            leg.AddEntry(hdata, "Background-subtracted #mu#muee data")
            leg.AddEntry(hmc, "#eta#rightarrow#mu#muee MC")
        else:
            leg.AddEntry(hdata, "#mu#muee data") 
            leg.AddEntry(hsides, "Background taken from the sidebands", 'l')
            leg.AddEntry(hrbkg, "#eta#rightarrow#mu#mu#gamma background MC", 'l')
            leg.AddEntry(hmc, "#eta#rightarrow#mu#muee MC #times #bar{B} #times %.2f"%(scale), 'l') 
            leg.AddEntry(htot, "Sum S+B", 'l') 
            htot.SetLineColor(ROOT.kRed)
            htot.SetLineWidth(2)
            hsides.SetLineStyle(2)
            hsides.SetLineColor(ROOT.kBlue)
            hrbkg.SetLineStyle(2)
            hrbkg.SetLineColor(ROOT.kMagenta)
            hmc.SetLineStyle(2)
            hmc.SetLineColor(ROOT.kOrange)
    else:
        if subtract:
            leg.AddEntry(hdata, "Sideband-subtracted #mu#mu data")
            leg.AddEntry(hmc, "#eta#rightarrow#mu#mu MC")
        else:
            leg.AddEntry(hdata, "#mu#mu data") 
            leg.AddEntry(hsides, "Background taken from the sidebands", 'l')
            leg.AddEntry(hmc, "#eta#rightarrow#mu#mu MC", 'l') 
            leg.AddEntry(htot, "Sum S+B", 'l') 
            htot.SetLineWidth(2)
            htot.SetLineColor(ROOT.kRed)
            hsides.SetLineStyle(2)
            hsides.SetLineColor(ROOT.kBlue)
            hmc.SetLineStyle(2)
            hmc.SetLineColor(ROOT.kOrange)
    if subtract:
        hdata.DrawNormalized("PE")
        hmc.DrawNormalized("same")
    else:
        htot.Draw("hist E")
        hdata.Draw("PE same")
        hsides.Draw("hist same")
        hrbkg.Draw("hist same")
        hmc.Draw("hist same") 
        htot.GetXaxis().SetRangeUser(0, 60)
    #hdata.Draw("PE")
    #hmc.Draw("same")
else:
    leg.AddEntry(hdata, "#mu#mu data")
    hdata.Draw("PE")
    hmc.Draw("same")
    leg.AddEntry(hmc, "#eta#rightarrow#mu#mu MC + %dK"%(int(miny/1000)))
#leg.AddEntry(hmc, "#eta#rightarrow#mu#mu MC + Bkg fit")
leg.Draw("hist same")

input("djhf")

if not compt:
    #now make a new pad underneath that shows the ratio of the two histograms
    c = ROOT.TCanvas()
    hpad = ROOT.TPad("hpad", "Pad for histograms", 0, 0.3, 1, 1.0)
    rpad = ROOT.TPad("rpad", "Pad for ratio", 0, 0.0, 1, 0.3)
    hpad.SetBottomMargin(0)
    rpad.SetTopMargin(0)
    rpad.SetBottomMargin(0.3)
    hpad.Draw()
    rpad.Draw()
    hpad.cd()
    hdata.Draw()
    hmc.Draw("same")
    leg.Draw("same")

    rpad.cd()
    hratio = hmc.Clone("hratio")
    hratio.Divide(hdata)
    hratio.SetLineWidth(2)
    hratio.SetLineColor(ROOT.kBlue)
    hratio.SetMarkerSize(1)
    hratio.GetYaxis().SetTitle("(MC+%dk)/Data"%(int(miny/1000)))
    hratio.Draw("PE")
    hratio.GetXaxis().SetTitle("m_{#mu#mu} (GeV)") 
    hratio.GetXaxis().SetLabelSize(0.1) 
    hratio.GetXaxis().SetTitleSize(0.1)
    hratio.SetTitle("Ratio") 
    hratio.SetStats(ROOT.kFALSE)
    hratio.GetYaxis().SetLabelSize(0.075)
    hratio.GetYaxis().SetTitleSize(0.075)
    hratio.GetYaxis().SetTitleOffset(0.5)

    #Now zoom in to just 0.52 to 0.57
    hratio.GetXaxis().SetRangeUser(0.52, 0.57)
    hdata.GetXaxis().SetRangeUser(0.52, 0.57)
    hdata.GetYaxis().SetRangeUser(miny-2000, miny+153000) 

    input("khlkjhl") 
else:
    
    c1.SaveAs("AN_Figures/DataMC_pt_%s.C"%("mumu" if not sig else "sig") )

