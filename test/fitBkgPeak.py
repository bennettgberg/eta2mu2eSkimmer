import ROOT, sys

distname = "hMlplp"

bkgfile = "bparking_bkgMCtest31.root"
f = ROOT.TFile.Open(bkgfile)
h = f.Get(distname)
h.Rebin(5)

nbins = h.GetNbinsX()
xmin = h.GetXaxis().GetXmin()
xmax = h.GetXaxis().GetXmax()

h.Draw()

c = ROOT.TCanvas("cnew","cnew")
c.cd()

rrv = ROOT.RooRealVar("m_{2#mu2e}", "", .48, .7) #xmin, xmax) 
rrv.setRange('peak', 0.51, 0.63)
rrv.setRange("full", xmin, xmax)

sys.path.append("tm_analysis/analysis/python/")
import utils.fit_function_library as library

#CrystalBall function?
#bkgModel = library.get_fit_function('CB', rrv)
#bkgModel.set_params( mcb=library.Param(.555, .545, .575), acb=library.Param(-4.5, -6, -0.1), ncb=library.Param(30, 15, 40), scb=library.Param(.0195, .01, .05) )
#Voigtian?
#bkgModel = library.get_fit_function('Voigtian', rrv)
#bkgModel.set_params( mv=library.Param(.555, .545, .575), wv=library.Param(.005, .001, .05), sv=library.Param(.005, .001, .05) )
#DoubleGaussian??
#bkgModel = library.get_fit_function('DoubleGauss', rrv)
#bkgModel.set_params( mg=library.Param(.555, .545, .575), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
#Landau function??
#bkgModel = library.get_fit_function('Landau', rrv)
#bkgModel.set_params( ml=library.Param(.555, .545, .575), sl=library.Param(.005, .001, .01) )
#CrystalBall + Gaussian??
bkgModel = library.get_fit_function('CB_Gauss', rrv)
bkgModel.set_params( mcb=library.Param(.555, .545, .575), acb=library.Param(-4.5, -6, -0.1), ncb=library.Param(30, 15, 40), scb=library.Param(.0195, .005, .05), sg=library.Param(0.01, 0.001, 0.1), CB_frac=library.Param(0.8, 0.0, 1.0) )


bkgMC = ROOT.RooDataHist("bkgMC", "bkgMC", rrv, ROOT.RooFit.Import(h)) 
nbkg = ROOT.RooRealVar("nbkg", "nbkg", 50, 1, 10000)

model = ROOT.RooExtendPdf("model", "extended pdf", bkgModel(), nbkg)
model.fitTo(bkgMC, ROOT.RooFit.Save())
frame = rrv.frame()
frame.SetTitle("")
frame.GetXaxis().SetTitle("m_{2#mu2e} [GeV]") 
binsize = (xmax - xmin) / nbins
frame.GetYaxis().SetTitle("Events / (%f GeV)"%(binsize))
#data.plotOn(frame, Name="Data", DrawOption="PEZ")
bkgMC.plotOn(frame, ROOT.RooFit.Name("bkgMC"), ROOT.RooFit.DrawOption("PEZ"))
model.plotOn(frame, ROOT.RooFit.Name("Bkg"), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineColor(ROOT.kBlue))
frame.Draw("AC")
chi2_SB = frame.chiSquare() #* ndf
print("***chi2/ndf of the fit: %f***"%(chi2_SB)) 
frame.SetName("")

c.Modified()
c.Update()
import utils.CMSStyle as cmsstyle
cmsstyle.setCMSLumiStyle(c, 0, era='2022')
c.Update()
input("press enter to continue bruv")
