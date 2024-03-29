import ROOT, sys

#include uncertainties on the event weights too (instead of just statistical)?
incWtUnct = False # True

#how many electrons to require elID: 0, 1, or 2 (both)?
req_elID = 2

#cut out .04 < M_ee < .09, or nah?
cut_mee = False

#distname = "hMlplp"
distname = "hMmmelel"
if incWtUnct:
    distnameUp = "hMUpmmelel"
    distnameDn = "hMDnmmelel"

if req_elID == 2:
    #bkgfile = "bparking_bkgMCtest31.root"
    #bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3826.root"
    #bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3827.root"
    #WP90 only!
    #bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3831.root"
    #elpt > 2 and 2022+2023
    #bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3837.root"
    ##one cand only and new xsec weight
    ##bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3844.root"
    #cut out .04 < M_ee < .09
    bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3851.root"
elif req_elID == 1:
    #bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3829.root"
    #bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3834.root"
    bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3838.root"
elif req_elID == 3:
    bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3839.root"
elif req_elID == 0:
    if cut_mee:
        #with .04 < M_ee < .09 cut out!
        bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3849.root"
    else:
        #normal
        #bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3850.root"
        #updated event weights
        bkgfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3854.root"
f = ROOT.TFile.Open(bkgfile)
h = f.Get(distname)

rebin = 5 #6 #4

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

#temporary fix for 2022 lumi only
if req_elID == 2 or req_elID == 1:
    h.Scale(38.48/(38.48+28.89)) 

h.Draw()

c = ROOT.TCanvas("cnew","cnew")
c.cd()

fitmin = .50
fitmax = .60
#rrv = ROOT.RooRealVar("m_{2#mu2e}", "", .45, .8) #xmin, xmax) 
rrv = ROOT.RooRealVar("m_{2#mu2e}", "", fitmin, fitmax) #xmin, xmax) 
rrv.setRange('peak', 0.51, 0.59)
rrv.setRange("full", xmin, xmax)

sys.path.append("tm_analysis/analysis/python/")
import utils.fit_function_library as library

##for very small res bkg, simple SingleGauss works good enough!
#bkgModel = library.get_fit_function('SingleGauss', rrv)
#nparam = 3
#bkgModel.set_params( mg=library.Param(.548, .545, .565), sg=library.Param(.01, .0001, .1) )
#CrystalBall function?
#bkgModel = library.get_fit_function('CB', rrv)
#bkgModel.set_params( mcb=library.Param(.555, .545, .575), acb=library.Param(-4.5, -6, -0.1), ncb=library.Param(30, 15, 40), scb=library.Param(.0195, .01, .05) )
#Breit-Wigner?
bkgModel = library.get_fit_function('BreitWigner', rrv)
nparam = 3
bkgModel.set_params( mb=library.Param(.555, .545, .575), wb=library.Param(.005, .001, .05) )
##Voigtian?
#bkgModel = library.get_fit_function('Voigtian', rrv)
#nparam = 4
#bkgModel.set_params( mv=library.Param(.555, .545, .575), wv=library.Param(.005, .001, .05), sv=library.Param(.005, .000001, .05) )
#DoubleGaussian??
#bkgModel = library.get_fit_function('DoubleGauss', rrv)
#bkgModel.set_params( mg=library.Param(.555, .545, .575), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
#Landau function??
#bkgModel = library.get_fit_function('Landau', rrv)
#bkgModel.set_params( ml=library.Param(.555, .545, .575), sl=library.Param(.005, .001, .01) )
##CrystalBall + Gaussian??
#bkgModel = library.get_fit_function('CB_Gauss', rrv)
#bkgModel.set_params( mcb=library.Param(.555, .545, .575), acb=library.Param(-4.5, -6, -0.1), ncb=library.Param(30, 15, 40), scb=library.Param(.0195, .005, .05), sg=library.Param(0.01, 0.001, 0.1), CB_frac=library.Param(0.8, 0.0, 1.0) )

##cheating and setting small bin weight to 0 to see if makes chi2 reasonable?-- yes it does....
#h.SetBinContent(59, 0)
#h.SetBinError(59, 0)

bkgMC = ROOT.RooDataHist("bkgMC", "bkgMC", rrv, ROOT.RooFit.Import(h)) 
nbkg = ROOT.RooRealVar("nbkg", "nbkg", 50, 1, 10000)

model = ROOT.RooExtendPdf("model", "extended pdf", bkgModel(), nbkg)
model.fitTo(bkgMC, ROOT.RooFit.Save())
frame = rrv.frame()
frame.SetTitle("")
frame.GetXaxis().SetTitle("m_{2#mu2e} [GeV]") 
binsize = (xmax - xmin) / nbins
frame.GetYaxis().SetTitle("Events / (%.3f GeV)"%(binsize))
#data.plotOn(frame, Name="Data", DrawOption="PEZ")
bkgMC.plotOn(frame, ROOT.RooFit.Name("bkgMC"), ROOT.RooFit.DrawOption("PEZ"))
model.plotOn(frame, ROOT.RooFit.Name("Bkg"), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineColor(ROOT.kBlue))
frame.Draw("AC")
ndf = int( (fitmax-fitmin)/binsize ) - nparam #nbins - 4 #7
chi2_SB = frame.chiSquare()
print("***chi2 of the fit: %f***"%(chi2_SB*ndf)) 
print("ndf: %d --> chi2/ndf=%f"%(ndf, chi2_SB)) 

text = "#chi^{2}/ndf = %.2f/%d = %.2f"%(chi2_SB*ndf, ndf, chi2_SB)

leg = ROOT.TLegend()
leg.SetHeader(text)
leg.Draw("same")

##make TPaveText box to show fit results!
#pav = ROOT.TPaveText(.5005, 3, .54, 8)
#pav.AddText("Fit results:")
#pav.AddText("#mu: %.4f +/- %.5f"%(model.pars['mg'].getVal(), model.pars['mg'].getError()))
#pav.AddText("#sigma: %.4f +/- %.5f"%(model.pars['sg'].getVal(), model.pars['sg'].getError()))
#pav.Draw("same")

frame.SetName("")

c.Modified()
c.Update()
import utils.CMSStyle as cmsstyle
cmsstyle.setCMSLumiStyle(c, 0, era='2022')
c.Update()
input("press enter to continue bruv")
