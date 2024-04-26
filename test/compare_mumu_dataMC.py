import ROOT

#fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3818_ALL.root")
#additional cuts, nd stuff
fdata = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3860_ALL.root")
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3854.root")
#additional cuts, nd stuff (no PU corrections)
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3867.root")
#with PU corrections
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3866.root")
#new evt weights (with DG/Cheb4 2mu fits) -- no PU corrections
#fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3877.root")
#yes PU
fmc = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3878.root")
#miny = 665000
miny = 460000

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
#add a constant 670K to each MC hist bin to align it with the data histogram
for i in range(hmc.GetNbinsX()):
    xval = hmc.GetBinCenter(i)
    #bkgi = params[0] + params[1]*xval + params[2]*xval*xval + params[3]*xval*xval*xval
    #bkgi = params[0] + params[1]*xval + params[2]*xval*xval + params[3]*xval*xval*xval + params[4]*xval*xval*xval*xval
    hmc.SetBinContent(i, hmc.GetBinContent(i) + miny)
    #hmc.SetBinContent(i, hmc.GetBinContent(i) + bkgi)

hdata.Draw("PE")
hmc.SetLineWidth(2)
hmc.Draw("same")

leg = ROOT.TLegend()
leg.AddEntry(hdata, "#mu#mu data")
leg.AddEntry(hmc, "#eta#rightarrow#mu#mu MC + %dK"%(int(miny/1000)))
#leg.AddEntry(hmc, "#eta#rightarrow#mu#mu MC + Bkg fit")
leg.Draw("hist same")

input("djhf")
