import ROOT, sys

weighted = True #False

#include uncertainties on the event weights too (instead of just statistical)?
incWtUnct = False

#how many electrons to require elID: 0, 1, or 2 (both)?
req_elID = 2

#distname = "hMlplp"
if weighted:
    distname = "hMmmelel"
    if incWtUnct:
        distnameUp = "hMUpmmelel"
        distnameDn = "hMDnmmelel"
else:
    distname = "hMNoWtmmelel"

#sigfile = "bparking_sigMCtest30.root"
#sigfile = "bparking_sigMCtest31.root"
#sigfile = "bparking_sigMCtest374.root"
if req_elID == 2:
    #nMiss==0, vProb>.5
    #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3826.root"
    #nMiss<=3, vProb>.1
    sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3827.root"
elif req_elID == 1:
    sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3829.root"
else:
    print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    exit()
f = ROOT.TFile.Open(sigfile)
h = f.Get(distname)
rebin = 5
h.Rebin(rebin)

nbins = h.GetNbinsX()
xmin = h.GetXaxis().GetXmin()
xmax = h.GetXaxis().GetXmax()

if incWtUnct:
    hUp = f.Get(distnameUp)
    hDn = f.Get(distnameDn)
    hUp.Rebin(rebin)
    hDn.Rebin(rebin)
    for j in range(nbins):
        binCon = h.GetBinContent(j)
        statErr = h.GetBinError(j)
        #wtErrUp = hUp.GetBinContent(j) + hUp.GetBinError(j) - binCon
        wtErrUp = hUp.GetBinContent(j) - binCon
        #wtErrDn = binCon - (hDn.GetBinContent(j) - hDn.GetBinError(j))
        wtErrDn = binCon - hDn.GetBinContent(j)
        #print("j=%d, binCenter= %f, binCon=%f, wtErrUp: %f, wtErrDn: %f"%(j, h.GetBinCenter(j), binCon, wtErrUp, wtErrDn)) 
        wtErr = max(wtErrUp, wtErrDn)
        #totErr = wtErr
        totErr = ((wtErr)**2 + (statErr)**2)**0.5
        h.SetBinError(j, totErr)

h.Draw()

c = ROOT.TCanvas("cnew","cnew")
c.cd()

fitmin = .505 #.50 #.505
fitmax = .585 #.60 #.585
#rrv = ROOT.RooRealVar("m_{2#mu2e}", "", .48, .7) #xmin, xmax) 
#rrv = ROOT.RooRealVar("m_{2#mu2e}", "", .48, .62) #xmin, xmax) 
rrv = ROOT.RooRealVar("m_{2#mu2e}", "", fitmin, fitmax) #xmin, xmax) 
#rrv.setRange('peak', 0.51, 0.63)
rrv.setRange('peak', fitmin, fitmax)
rrv.setRange("full", xmin, xmax)

sys.path.append("tm_analysis/analysis/python/")
import utils.fit_function_library as library

##CrystalBall function?
#sigModel = library.get_fit_function('CB', rrv)
#nparam = 5
#sigModel.set_params( mcb=library.Param(.548, .543, .553), acb=library.Param(-1.5, -100, -0.01), ncb=library.Param(30, .1, 1000), scb=library.Param(.0195, .005, .05) )
##sigModel.set_params( mcb=library.Param(.548, .543, .553), acb=library.Param(-4.5, -6, -0.1), ncb=library.Param(40, .1, 400), scb=library.Param(.0195, .01, .05) )
##Breit-Wigner?
#sigModel = library.get_fit_function('BreitWigner', rrv)
#nparam = 3
#sigModel.set_params( mb=library.Param(.548, .543, .553), wb=library.Param(.005, .0005, .05) )
##Voigtian?
#sigModel = library.get_fit_function('Voigtian', rrv)
#nparam = 4
#sigModel.set_params( mv=library.Param(.548, .543, .553), wv=library.Param(.005, .0005, .05), sv=library.Param(.005, .0000001, .05) )
##DoubleGaussian??
#sigModel = library.get_fit_function('DoubleGauss', rrv)
#nparam = 5
#sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
##Landau function??
#sigModel = library.get_fit_function('Landau', rrv)
#nparam = 3
#sigModel.set_params( ml=library.Param(.548, .543, .553), sl=library.Param(.005, .001, .01) )
##CrystalBall + Gaussian??
#sigModel = library.get_fit_function('CB_Gauss', rrv)
#nparam = 7
#sigModel.set_params( mcb=library.Param(.548, .543, .553), acb=library.Param(-4.5, -6, -0.1), ncb=library.Param(30, 5, 1000), scb=library.Param(.0195, .005, .05), sg=library.Param(0.01, 0.001, 0.1), CB_frac=library.Param(0.5, 0.05, 0.95) )
##DoubleCB?
#sigModel = library.get_fit_function('DoubleCB', rrv)
#nparam = 7
#sigModel.set_params( mcb=library.Param(.548, .543, .553), acb1=library.Param(-.508, -6, -0.1), acb2=library.Param(-.470, -6, -0.1), ncb=library.Param(30, .1, 1000), scb1=library.Param(.01, .005, .05), scb2=library.Param(.02, .005, .05) )
##DoubleCB + Gauss
#sigModel = library.get_fit_function('DoubleCB_Gauss', rrv)
#nparam = 8
#sigModel.set_params( mcb=library.Param(.548, .547, .553), acb2=library.Param(0.265, 0, 5.0), acb1=library.Param(-0.257, -6, -0.1), scb1=library.Param(.0145, .001, .05), scb2=library.Param(.0101, .001, .09), sg=library.Param(.0144, .001, .05), sig1frac=library.Param(.440, .01, .99) )
#TripleGaussian??
sigModel = library.get_fit_function('TripleGauss', rrv)
nparam = 7
sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.03, .001, .5), sg2=library.Param(.007, .001, .5), sg3=library.Param(.001, .0001, .01), sig1frac=library.Param(0.8, 0.1, 0.9), sig2frac=library.Param(0.2, 0.05, 0.95) )


sigMC = ROOT.RooDataHist("sigMC", "sigMC", rrv, ROOT.RooFit.Import(h)) 
nsig = ROOT.RooRealVar("nsig", "nsig", 50, 1, 10000)

model = ROOT.RooExtendPdf("model", "extended pdf", sigModel(), nsig)
model.fitTo(sigMC, ROOT.RooFit.Save())
frame = rrv.frame()
frame.SetTitle("")
frame.GetXaxis().SetTitle("m_{2#mu2e} [GeV]") 
binsize = (xmax - xmin) / nbins
if weighted:
    frame.GetYaxis().SetTitle("Weighted Events / (%.3f GeV)"%(binsize))
else:
    frame.GetYaxis().SetTitle("Unweighted Events / (%.3f GeV)"%(binsize))
#data.plotOn(frame, Name="Data", DrawOption="PEZ")
sigMC.plotOn(frame, ROOT.RooFit.Name("sigMC"), ROOT.RooFit.DrawOption("PEZ"))
model.plotOn(frame, ROOT.RooFit.Name("Bkg"), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineColor(ROOT.kBlue))
frame.Draw("AC")
#now also overlay the unweighted histogram
#hMNoWt = f.Get("hMNoWtlplp")
#hMNoWt.Rebin(5)
#hMNoWt.Draw("hist same")
##########################################
ndf = int( (fitmax-fitmin)/binsize ) - nparam #nbins - 4 #7
chi2_SB = frame.chiSquare()
print("***chi2 of the fit: %f***"%(chi2_SB*ndf)) 
print("ndf: %d --> chi2/ndf=%f"%(ndf, chi2_SB)) 

text = "#chi^{2}/ndf = %.2f/%d = %.2f"%(chi2_SB*ndf, ndf, chi2_SB)

leg = ROOT.TLegend()
leg.SetHeader(text)
leg.Draw("same")

frame.SetName("")

c.Modified()
c.Update()
import utils.CMSStyle as cmsstyle
cmsstyle.setCMSLumiStyle(c, 0, era='2022')
c.Update()
input("press enter to continue bruv")
