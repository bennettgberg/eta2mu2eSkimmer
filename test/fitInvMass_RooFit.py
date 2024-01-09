#get the data from infile, fit to const bkg + gaussian.
#from ROOT import TFile, TF1, TH1F, TCanvas, RooRealVar, RooDataHist, RooFit
import ROOT
import sys
#sys.path.append("../../tm_analysis/analysis/python/utils")
#sys.path.append("/uscms_data/d3/bgreenbe/CMSSW_10_2_9/src/tm_analysis/analysis/python")
sys.path.append("tm_analysis/analysis/python/")
#import fitter
import utils.fitter as fitter
#import utils.fitter_2mu2e as fitter
from ctypes import c_double

#how many electrons to require elID: 0, 1, or 2 (both)?
req_elID = 2 #1

#infile = "hipTpair/data_DoubleMuon_Run2_skimmed.root"
#infile = "multiquad/data_DoubleMuon_Run2_skimmed.root"
#infile = "multimmg/data_DoubleMuon_Run2_skimmed.root"
#infile = "oppocharge/data_DoubleMuon_Run2_skimmed.root"
#infile = "data_DoubleMuon_Run2_ntuple_skimmed.root"
#infile = "/eos/uscms/store/user/bgreenbe/BParking2022/bparking_test6_ALL.root"
#infile = "/eos/uscms/store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest7_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest7_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest32_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest321_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest322_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest323_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/bparking_datatest335_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest33_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest323_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest335_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3310_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest36_ALL.root"
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest369_ALL.root" #elID req'd
#infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest374_ALL.root" #elID req'd
if req_elID == 2:
    infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3819_ALL.root" #elID req'd
elif req_elID == 1:
    infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3829_ALL.root" #elID req'd on ONE ele
#distname = "Mmmee"
#distname = "Mmmg"
#distname = "hMe"
#distname = "hMlpe"
#distname = "hMlplp"
#distname = "hMmmlplp"
distname = "hMmmelel"
#accidentally used the wrong name!
#distname = "Ptmmee"
compare_MC = True

f = ROOT.TFile.Open(infile)
f.Print()

#t = f.Get("Events") 
h = f.Get(distname)
h.Print()
#h = f.Get("hMlpe")

h.Rebin(5) #15
h.Print()

#nbins = 40 #30 #25 
#xmin = .25 #.4 #.8
#xmax = .9 #1.2
nbins = h.GetNbinsX()
xmin = h.GetXaxis().GetXmin()
xmax = h.GetXaxis().GetXmax()

#initiate histogram
#h = ROOT.TH1F("h", "mu mu e e", nbins, xmin, xmax)
#fill histogram fast
#t.Draw(distname + ">>h")
h.Draw()

c = ROOT.TCanvas("cnew","cnew")
c.cd()
#cgfit = TF1("cgfit", "([0] + [4]*x + gaus(1))", xmin, xmax);
#cgfit = TF1("cgfit", "([0] + [6]*x + [7]*x*x + crystalball(1))", xmin, xmax);
#cgfit = TF1("cgfit", "([0] + gaus(1))", xmin, xmax);
#cgfit.Print()


#signal shape is best described by a crystalBall.
#rrv = ROOT.RooRealVar("m_{2#mu2e}","m_{2#mu2e} [GeV]",xmin,xmax)
#rrv = ROOT.RooRealVar("m_{2#mu2e}","",xmin,xmax)
rrv = ROOT.RooRealVar("m_{2#mu2e}","",.45,xmax)
#rrv = ROOT.RooRealVar("m_{2#mu2e}","m_{2#mu2e} [GeV]",0.5,0.9)
#rrv.setRange("peak", 0.53, 0.57)
#rrv.setRange("peak", 0.51, 0.63)
rrv.setRange("peak", 0.51, 0.60)
#rrv.setRange("full", xmin, xmax)
rrv.setRange("full", .45, xmax)
#try a different bkg model!!
#myfitter = fitter.fitter_4mu(mass=rrv, bkg_model='Threshold', sig_model='CB')
#myfitter = fitter.fitter_4mu(mass=rrv, bkg_model='Cheb2', sig_model='CB')
#myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb2', sig_model='CB')
#myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb2', sig_model='Voigtian')
##Voigtian with const combinatorial bkg!!
#myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb0', sig_model='Voigtian')
##Bkg model is Gaussian (resonant bkg--already fitted!) + constant !
#myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='ConstGauss', sig_model='Voigtian')
if req_elID == 2:
    #myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb0', resBkg_model='SingleGauss', sig_model='Voigtian')
    #myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb0', resBkg_model='BreitWigner', sig_model='Voigtian')
    #myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb0', resBkg_model='BreitWigner', sig_model='TripleGauss')
    myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb0', resBkg_model='BreitWigner', sig_model='BreitWigner')
elif req_elID == 1:
    #myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb1', resBkg_model='Voigtian', sig_model='Voigtian')
    myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb1', resBkg_model='BreitWigner', sig_model='Voigtian')
##Double-Gaussian with const combinatorial bkg!!
#myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb0', sig_model='DoubleGauss')
#myfitter = fitter.fitter_4mu(mass=rrv, bkg_model='Cheb2', sig_model='Gauss2')
#now set the parameters appropriately.
#from collections import namedtuple
#import fit_function_library as library
import utils.fit_function_library as library
#try to do like this: 'Param', ['val', 'min', 'max']
#myfitter.set_bkg_params( alpha=library.Param(1, 0.5, 5), x0=library.ConstParam(.2122) ) #for Threshold
#myfitter.set_bkg_params( a1=library.Param(0.5, -1, 1), a2=library.Param(0.2, -1, 1) ) #for Cheb2
#myfitter.set_bkg_params( a1=library.Param(0.5, -1, 1) ) #for Cheb1
if req_elID == 2:
    myfitter.set_bkg_params( ) #a1=library.Param(0.5, .1, 1) ) 
elif req_elID == 1:
    myfitter.set_bkg_params( a1=library.Param(0.5, .1, 1) ) 
#myfitter.set_resBkg_params( mg=library.Param(.545, .544, .549), sg=library.Param(2.64e-2, 2.63e-2, 2.65e-2) )
if req_elID == 2:
    #myfitter.set_resBkg_params( mg=library.ConstParam(.555, .554, .556), sg=library.Param(1.82e-2, 1.80e-2, 1.84e-2) )
    myfitter.set_resBkg_params( mb=library.ConstParam(.552), wb=library.ConstParam(1.87e-2) )
    myfitter.set_resBkg_norm(21.3)
elif req_elID == 1:
    #myfitter.set_resBkg_params( mv=library.Param(.557, .556, .558), sv=library.Param(1.01e-3, 1.00e-3, 1.02e-3), wv=library.Param(2.97e-2, 2.95e-2, 2.99e-2) )
    myfitter.set_resBkg_params( mv=library.ConstParam(.557), sv=library.ConstParam(1.01e-3), wv=library.ConstParam(2.99e-2) )
    myfitter.set_resBkg_norm(199.4)
#CrystalBall function for signal model
#myfitter.set_sig_params( mcb=library.Param(.549, .5, .6), acb=library.Param(-4.5, -6, -0.5), ncb=library.Param(21, 15, 25), scb=library.Param(.0195, .01, .02) )
#Voigtian function
#myfitter.set_sig_params( mv=library.Param(.548, .543, .553), wv=library.Param(.005, .0001, .05), sv=library.Param(.005, .0001, .05) )
##set constant those params found from the fit to signal MC
#myfitter.set_sig_params( mv=library.ConstParam(.5475), wv=library.ConstParam(.0141) ) #, sv=library.ConstParam(.00306)))
#including weight uncertainties in the error bars
#myfitter.set_sig_params( mv=library.Param(.5475, .5474, .5476), wv=library.Param(.0141, .0140, .0142), sv=library.Param(.00306, .00305, .00307))
#not including weight unctys in the error bars
#myfitter.set_sig_params( mv=library.Param(.5482, .5481, .5483), wv=library.Param(.0172, .01, .023), sv=library.Param(.0005, .0001, .005))
if req_elID == 2:
    #myfitter.set_sig_params( mv=library.Param(.5472, .5470, .548), wv=library.Param(.0172, .01, .023), sv=library.Param(.00151, .001, .01))
    #Voigtian
    #myfitter.set_sig_params( mv=library.ConstParam(.5472), wv=library.ConstParam(.0172), sv=library.ConstParam(.00151) )
    #BW
    myfitter.set_sig_params( mb=library.Param(.5473, .5465, .5481), wb=library.Param(.01609, .01415, .01803) )
    ##TripleGauss function
    #myfitter.set_sig_params( mg=library.ConstParam(0.5475), sg1=library.ConstParam(0.3875), sg2=library.ConstParam(0.01227), sg3=library.ConstParam(6.1425e-3), sig1frac=library.ConstParam(0.2638), sig2frac=library.ConstParam(0.1917) )
elif req_elID == 1:
    #myfitter.set_sig_params( mv=library.Param(.5479, .5478, .5480), wv=library.Param(.0151, .01, .025), sv=library.Param(.00050, .00049, .00051))
    myfitter.set_sig_params( mv=library.ConstParam(.5479), wv=library.ConstParam(.0151), sv=library.Param(.00050, .00005, .05) )
#DoubleGauss function
#myfitter.set_sig_params( mg=library.Param(0.547612, 0.5475, 0.5477), sg1=library.Param(0.0234698, 0.0234, 0.0235), sg2=library.Param(0.00645498, 0.00645, 0.00646), sig1frac=library.Param(0.3826, 0.382, 0.383) )
#myfitter.set_sig_params( mg=library.Param(0.547612, 0.5475, 0.5477), sg1=library.Param(0.0234698, 0.0234, 0.0235), sg2=library.Param(0.00645498, 0.0001, 0.05), sig1frac=library.Param(0.3826, 0.382, 0.383) )

#new signal model: 2xGauss (NOT DoubleGauss)
#myfitter.set_sig_params( mg1=library.Param(.549, .545, .552), sg1=library.Param(.01, .005, .03), mg2=library.Param(.555, .545, .58), sg2=library.Param(.01, .005, .03), sig1frac=library.Param(0.5, 0.1, 5.0) )

data = ROOT.RooDataHist("data", "data", rrv, ROOT.RooFit.Import(h))

#do the fit
fitres = myfitter.model.fitTo(data, ROOT.RooFit.Save())
#print("fitres: " + str(fitres)) 


#plot the fit
#frame = rrv.frame()
frame = rrv.frame(ROOT.RooFit.Name(f"pullframe"), ROOT.RooFit.Title(" "))
frame.SetTitle("")
frame.GetXaxis().SetTitle("m_{2#mu2e} [GeV]")
binsize = (xmax - xmin) / nbins
frame.GetYaxis().SetTitle("Events / (%.3f GeV)"%(binsize))
#data.plotOn(frame, Name="Data", DrawOption="PEZ")
data.plotOn(frame, ROOT.RooFit.Name("Data"), ROOT.RooFit.DrawOption("PEZ"))

myfitter.model.plotOn(frame, ROOT.RooFit.Name("Bkg"), ROOT.RooFit.Components('bkg'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(ROOT.kGreen-1))
if req_elID == 2:
    myfitter.model.plotOn(frame, ROOT.RooFit.Name("ResBkg"), ROOT.RooFit.Components('resbkg'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(ROOT.kCyan))
elif req_elID == 1:
    myfitter.model.plotOn(frame, ROOT.RooFit.Name("ResBkg"), ROOT.RooFit.Components('resbkg'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(ROOT.kCyan))

#myfitter.model.plotOn(frame, ROOT.RooFit.Name("Sig"), ROOT.RooFit.Components('sig'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(4), ROOT.RooFit.LineColor(ROOT.kRed+1))
#myfitter.model.plotOn(frame, ROOT.RooFit.Name("Sig"), ROOT.RooFit.Components('CB'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(4), ROOT.RooFit.LineColor(ROOT.kRed+1))
myfitter.model.plotOn(frame, ROOT.RooFit.Name("Sig"), ROOT.RooFit.Components('sig'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(4), ROOT.RooFit.LineColor(ROOT.kRed+1))
#myfitter.model.plotOn(frame, ROOT.RooFit.Name("Sig"), ROOT.RooFit.Components('sig'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(ROOT.kRed-1))

myfitter.model.plotOn(frame, ROOT.RooFit.Name("Tot"), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineColor(ROOT.kBlue))
####
#minX = .45
#maxX = .8
##trying to add pull frame underneath! #####
##pullframe = rrv.frame(ROOT.RooFit.Name(f"pullframe"), ROOT.RooFit.Title(" "))
#print("Bkg exists??? : " + str(frame.findObject("Bkg"))) 
##pullframe.GetYaxis().SetRangeUser(-3,7)
#pullh = frame.pullHist("Data", "Bkg")
#print("pullh: " + str(pullh)) 
#####
##myfitter.model.plotOn(frame, ROOT.RooFit.Name("Bkg"), ROOT.RooFit.Components({myfitter.bkg}), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle('--'), ROOT.RooFit.LineColor(ROOT.kGreen-1))
##myfitter.model.plotOn(frame, ROOT.RooFit.Name("Sig"), ROOT.RooFit.Components({myfitter.sig}), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle('-.'), ROOT.RooFit.LineColor(ROOT.kRed+1))
#frame.Print()
frame.Draw("AC")

#ndf = 40 - 4 - 3
#ndf = nbins - 4 - 3
chi2_SB = frame.chiSquare() #* ndf
print("***chi2/ndf of the fit: %f***"%chi2_SB) 
#??
frame.SetName("")

c.Modified()
c.Update()

#ndata = h.Integral(h.FindBin(0.51), h.FindBin(0.63))
ndata = h.Integral(h.FindBin(0.51), h.FindBin(0.60))
#ndata = h.Integral(h.FindBin(0.90), h.FindBin(1.0))
argset = ROOT.RooArgSet(rrv)
sig_int = myfitter.sig.func.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
bkg_int = myfitter.bkg.func.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
resBkg_int = myfitter.resBkg.func.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
tot_int = myfitter.model.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
print("sigvals: %f, %f; bkgvals: %f, %f"%(sig_int.getVal(), myfitter.nsig.getVal(), bkg_int.getVal(), myfitter.nbkg.getVal())) 
nsig = sig_int.getVal() * myfitter.nsig.getVal()
nbkg = bkg_int.getVal() * myfitter.nbkg.getVal()
nrbkg = resBkg_int.getVal() * myfitter.nrbkg.getVal()
ntot = tot_int.getVal() * (myfitter.nsig.getVal() + myfitter.nbkg.getVal())

bkghi = bkg_int.getVal() * myfitter.nbkg.getAsymErrorHi() 
bkglo = bkg_int.getVal() * myfitter.nbkg.getAsymErrorLo() 

print("bkg errors: %f, %f "%(bkghi, bkglo)) 

#leg = ROOT.TLegend(0.4, 0.65, 0.8, 0.85)
#leg = ROOT.TLegend(0.15, 0.65, 0.45, 0.85)
leg = ROOT.TLegend(0.40, 0.50, 0.85, 0.85)
#leg.SetHeader("N: m_{2#mu2e} #in [0.51, 0.63] GeV")
leg.SetHeader(".51 < m_{2#mu2e} < 0.63 GeV")
#leg.SetHeader("N: m_{2#mu2e} #in [0.9, 1.0] GeV")
leg.SetLineWidth(0)
leg.AddEntry("Sig", f"Signal (CB): N = {nsig:.1f}", "l")
leg.AddEntry("Bkg", f"Comb. Background: N = {nbkg:.1f}", "l")
leg.AddEntry("ResBkg", f"Res. Background: N = {nrbkg:.1f}", "l")
leg.AddEntry("Tot", f"Sum: N = {ntot:.1f}", "l")
leg.AddEntry("Data", f"Data, N = {ndata:n}", "lep")
leg.Draw()
####################################
###CrystalBall starting params####
#cgfit.SetParameter(0, 0.3) #const
#cgfit.SetParameter(1, 8)  #norm
#cgfit.SetParameter(2, 0.548) #mean
##fix eta mass or nah??
##cgfit.FixParameter(2, .5479)
#cgfit.SetParameter(3, 0.018) #sigma
#
#cgfit.SetParameter(4, -2) #alph
#cgfit.SetParameter(5, 2) #n
#cgfit.SetParameter(6, .01) #slope
#cgfit.SetParameter(7, .01) #quadratic coef

#initiate histogram
#h = TH1F("h", "mu mu gamma", nbins, xmin, xmax)

#xax = h.GetXaxis()
#xax.Print()
#xax.SetTitle("4-lepton invariant mass (GeV)")
##xax.SetTitle("dimu-gamma invariant mass (GeV)")
#yax = h.GetYaxis()
#binsize = (xmax - xmin) / nbins
#yax.SetTitle("Events / %f GeV"%binsize)
#
#
#draw = True # False
#fitresult = h.Fit("cgfit", "LSB")
#
#c = TCanvas("c", "c")
#c.cd()
#h.Draw()
#
#params = fitresult.GetParams()
#const = params[0]
#slope = params[6]
##slope = 0
#mean = params[2]
#sigma = params[3]
#
##area under the bkg only
#hi = .63 # mean+3*sigma
#lo = .51 #mean-3*sigma
#bkgInt = (0.5*slope*(hi**2 - lo**2) + const*(hi - lo) + 1.0/3 * params[7]*(hi**3 - lo**3) ) / binsize
#
#sigInt = cgfit.Integral(lo, hi) / binsize - bkgInt
#
#print("About %f signal events and %f +/- %f background events under the peak."%(sigInt, bkgInt, sigmaB))
#

#use official CMS Style.
import utils.CMSStyle as cmsstyle
#cmsstyle.setCMSLumiStyle(c, 0, era='Run2')
cmsstyle.setCMSLumiStyle(c, 0, era='2022')

from math import log
##l = log(1 + sigInt/bkgInt) 
##print("l: " + str(l))
##l2 = l*(sigInt + bkgInt)
##l3 = l2 - sigInt
##print("l3: " + str(l3)) 
#
##approximate significance (neglecting sigmaB)
#significance = (2*((sigInt + bkgInt)*log(1 + sigInt/bkgInt) - sigInt))**0.5

#significance taking sigmaB into account
sigInt = nsig
bkgInt = nbkg
sigmaB = max(bkglo, bkghi)
#significance = (2*((sigInt + bkgInt)*log((sigInt+bkgInt)*(bkgInt+sigmaB**2)/(bkgInt**2 + (sigInt+bkgInt)*sigmaB**2)) - bkgInt**2/sigmaB**2 *log(1 + sigmaB**2*sigInt/(bkgInt*(bkgInt+sigmaB**2)))))**0.5 
#
#print("Significance is about %f sigmas."%(significance))

mc_sig_int = -1
mc_bkg_int = -1
#draw MC first (if drawing it)
if compare_MC:
    #fMC = ROOT.TFile.Open("EtaTo2Mu2E_2018_0_skimmedMCtest.root")
    #fMC = ROOT.TFile.Open("EtaTo2Mu2E_2018_0_ntuple_skimmedsignalMC.root")
    #tMC = fMC.Get("Events")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest141.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest17.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest19.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest231.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest32.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest321.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest322.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest323.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest324.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest335.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest339.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest3310.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest3313.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest339.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest351.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest352.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest36.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest369.root")
    #fMC = ROOT.TFile.Open("bparking_sigMCtest374.root")
    #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3819.root")
    if req_elID == 2:
        #nMiss==0, vProb>.5
        #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3826.root")
        #nMiss<=3, vProb>.1
        fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3827.root")
    elif req_elID == 1:
        fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3829.root")
    #tMC = fMC.Get("hMe")
    tMC = fMC.Get(distname)

    #MC_scales = [2.0, 1.0, 0.5]
    #MC_scales = [1.0]
    MC_scales = [0.3]
    hMC = [None for mcs in MC_scales]
    #draw the blinded MC scaled by a few different values
    for i,scale in enumerate(MC_scales):
        #hMC[i] = ROOT.TH1F("hMC"+str(i), "#mu#muee", nbins, xmin, xmax)
        hMC[i] = tMC.Clone()
        hMC[i].SetLineWidth(2)
        color = 7+i
        #10 is just white \throwingUpEmoji
        if color >= 9: color += 2
        #hMC[i].SetLineColor(color)
        hMC[i].SetLineColor(ROOT.kOrange)
        #switch to a temporary canvas so don't draw bs on the good canvas??
        #ctemp = TCanvas("ctemp", "ctemp")
        #ctemp.cd()
        drawopt = "hist same"
        #if i > 0: drawopt += " same"
        #tMC.Draw(distname + ">>hMC"+str(i), str(scale)+"*Weight*(Weight>0)", drawopt)
        rebin = hMC[i].GetNbinsX() / nbins
        print("rebin: " + str(rebin))
        irebin = int(rebin)
        hMC[i].Rebin(irebin)
        hMC[i].Scale(scale)
        hMC[i].Draw("hist same")
        #c.cd()
        #hMC.Draw("hist same")
        if i == 0:
            mc_sig_int = hMC[i].Integral(hMC[i].FindBin(0.51), hMC[i].FindBin(0.63))
#if compare_MC:
#    for i,scale in enumerate(MC_scales):
        leg.AddEntry(hMC[i],"#eta#rightarrow2#mu2e MC Sig * #bar{B} * %.1f: N = %.1f"%(scale, mc_sig_int), "l")
    c.Update()

    #fBkg = ROOT.TFile.Open("EtaToMuMuGamma_2018_0_ntuple_skimmedBkgMC.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest17.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest19.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest231.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest32.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest321.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest322.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest323.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest324.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest335.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest339.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest3310.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest3313.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest339.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest352.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest36.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest369.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest374.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest382.root")
    #fBkg = ROOT.TFile.Open("bparking_bkgMCtest3819.root")
    #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3819.root")
    if req_elID == 2:
        #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3826.root")
        #vProb>.1, nMiss<=3
        fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3827.root")
    elif req_elID == 1:
        fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3829.root")
    #tBkg = fBkg.Get("Events") 
    #hBkg = ROOT.TH1F("hBkg", "#mu#mu#gamma", nbins, xmin, xmax)
    #hBkg = fBkg.Get("hMe")
    hBkg = fBkg.Get(distname)
    rebin = hBkg.GetNbinsX() / nbins
    print("bkg rebin: " + str(rebin))
    irebin = int(rebin)
    hBkg.Rebin(irebin)
    hBkg.SetLineWidth(2)
    hBkg.SetLineColor(ROOT.kViolet)
    #tBkg.Draw(distname + ">>hBkg", "Weight*(Weight>0)", "hist same")
    hBkg.Draw("hist same")
    #mc_bkg_int = hBkg.Integral(hBkg.FindBin(0.51), hBkg.FindBin(0.63))
    mc_bkg_err_db = c_double(-1.0)
    mc_bkg_int = hBkg.IntegralAndError(hBkg.FindBin(0.51), hBkg.FindBin(0.63), mc_bkg_err_db)
    print("bkg int:")
    print(mc_bkg_int)
    mc_bkg_err = float(mc_bkg_err_db.value) 
    leg.AddEntry(hBkg, "#eta#rightarrow#mu#mu#gamma MC Bkg: N = %.1f +/- %.1f"%(mc_bkg_int, mc_bkg_err), "l")


if mc_bkg_int > 0: 
    sigInt -= mc_bkg_int
    bkgInt = (mc_bkg_int**2 + bkgInt**2)**0.5
significance = (2*((sigInt + bkgInt)*log((sigInt+bkgInt)*(bkgInt+sigmaB**2)/(bkgInt**2 + (sigInt+bkgInt)*sigmaB**2)) - bkgInt**2/sigmaB**2 *log(1 + sigmaB**2*sigInt/(bkgInt*(bkgInt+sigmaB**2)))))**0.5 
print("Significance is about %f sigmas."%(significance))

c.Update()

#minX = .45
#maxX = .8
##trying to add pull frame underneath! #####
#pullframe = rrv.frame(ROOT.RooFit.Name(f"pullframe"), ROOT.RooFit.Title(" "))
##pullframe.GetYaxis().SetRangeUser(-3,7)
#pullh = frame.pullHist("Data", "Bkg")
##print("pullh: " + str(pullh)) 
#residh = frame.residHist("Data", "Bkg")
#binning = rrv.getBinning()
#for i in range(1, h.GetNbinsX()+1):
#    rrv.setRange("range_for_bin", binning.binLow(i), binning.binHigh(i))
#    normset = ROOT.RooFit.NormSet(rrv)
#    bkgPdfIntegral = myfitter.bkg().createIntegral(ROOT.RooArgSet(rrv), "range_for_bin")
#    #print("bkgPdfIntegral: " + str(bkgPdfIntegral)) 
#    pset = ROOT.RooArgList()
#    pset.add(bkgPdfIntegral) 
#    bkgYield = ROOT.RooProduct("bkgYield", "bkgYield", pset)
#    fitError = bkgYield.getPropagatedError(fitres)
#    old_sigma_i = residh.GetPointY(i)/pullh.GetPointY(i)
#    new_sigma_i = (old_sigma_i**2 - fitError**2)**(1/2)
#    residh.SetPointY(i, residh.GetPointY(i)/new_sigma_i)
#    residh.SetPointEYhigh(i, residh.GetErrorYhigh(i)/new_sigma_i)
#    residh.SetPointEYlow(i, residh.GetErrorYlow(i)/new_sigma_i)
#pullframe.addPlotable(residh, 'P')
##pullframe.GetYaxis().SetRangeUser(-5,8)
#pullframe.GetXaxis().SetTitle("m_{#mu#muee} (GeV)")
#pullframe.GetYaxis().SetTitle("Pull") 
###pullframe.addPlotable(pullh, 'P')
#line = ROOT.TLine(minX, 0, maxX, 0)
#line.SetLineColor(ROOT.kGray)
#pullframe.addObject(line)
#c1 = cmsstyle.get_fit_canvas(frame, pullframe)
#pad1 = c1.FindObject("pad1")
#pad1.cd()
#pullframe.Draw("same")
#for i,scale in enumerate(MC_scales):
#    hMC[i].Draw("hist same")
#hBkg.Draw("hist same")
#leg.Draw("same")
#c1.Update()
#############################3

if compare_MC:
    print("MC integrals:\nSignal: %f\nBackground: %f"%(mc_sig_int, mc_bkg_int)) 

input("h")

##now do the fit again, but this time without the signal model.
## so make the sig_model a Pol1, but then make it 0.
#newfitter = fitter.fitter_4mu(mass=rrv, bkg_model='Threshold', sig_model='Pol1')
#newfitter.set_bkg_params( alpha=library.Param(1, 0.5, 5), x0=library.ConstParam(.2122) )
#newfitter.set_sig_params( p0=library.ConstParam(0), p1=library.ConstParam(0) )
#newfitter.sig.args['p0'].setConstant(True)
#newfitter.sig.args['p1'].setConstant(True)
#
#newfitres = newfitter.model.fitTo(data, ROOT.RooFit.Save())
##number of degrees of freedom for the fit with NO signal model
#newndf = 40 - 3
#chi2_B = frame.chiSquare() * newndf
#print("***chi2 of the fit with NO signal model: %f***"%chi2_B) 
##difference in number of degrees of freedom betwixt the fit with CB and the fit with no signal model
#dNdf = 4
#dChi2 = (chi2_B - chi2_SB) / dNdf
#print("dChi2 of the two fits: %f"%(dChi2)) 
