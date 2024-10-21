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

year = 2022 #2223

#how many electrons to require elID: 0, 1, 2 (both), or 3 (not actually 3 but just means WP80 required on both)?
req_elID = 2

#require muon ID or nah?
req_muID = True

#True if want to use the selection where events with .04 < M_ee < .09 GeV are cut
cut_mee = False

#include pileup corrections or nah
do_pileup = True

#include trigger efficiency corrections or nah?
do_trigCor = True

#use new event weights from DG/Cheb4 2mu fits?
new_wt = 1 #True

#use lowPt electrons instead of regular ones?
useLowPt = False

#include pre/post-EE distinction?
#incEE = True

#use modified invar mass distribution (for electron efficiency calculation) ? 
modnum = -1
if len(sys.argv) > 1:
    modnum = int(sys.argv[1]) 

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
    #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3819_ALL.root" #elID req'd
    #vProb>.1, nMiss<=3
    #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3818_ALL.root" #elID req'd
    #WP90 only
    #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3831_ALL.root"
    if year == 2223:
        #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3837_allRun3.root"
        infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3891_allRun3.root"
    if cut_mee:
        #cut .04 < M_ee < .09
        infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3851_ALL.root"
    else:
        if req_muID:
            #medium ID
            #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3883_ALL.root"
            #loose ID
            #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3884_ALL.root"
            #nominal!!
            #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38107_ALL.root"
            #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38126_ALL.root"
            #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38130_ALL.root"
            #reprocessed, nd stuff
            #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest4715_ALL.root"
            if useLowPt:
                infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest4734_ALL.root"
            else:
                infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest4733_ALL.root"
        else:
            #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3837_ALL.root"
            #including cuts on muon pt, eta, etc
            #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3860_ALL.root"
            infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3882_ALL.root"
    #electron pT > 2 AND 2023 data included
elif req_elID == 1:
    if year != 2022:
        print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa 2022 only aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
        sys.exit()
    #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3829_ALL.root" #elID req'd on ONE ele
    infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3838_ALL.root" #elID req'd on ONE ele
elif req_elID == 3:
    #tight (WP80) elID req'd on both electrons
    #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3839_ALL.root"
    #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3871_ALL.root"
    #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38122_ALL.root"
    infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest38138_ALL.root"
elif req_elID == 0:
    if cut_mee:
        #cut out .04 < M_ee < .09
        infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3849_ALL.root"
    else:
        #normal
        #infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3850_ALL.root"
        infile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3858_ALL.root"
    
#distname = "Mmmee"
#distname = "Mmmg"
#distname = "hMe"
#distname = "hMlpe"
#distname = "hMlplp"
#distname = "hMmmlplp"
if modnum > -1:
    distname = "hMMod"+str(modnum)
else:
    if useLowPt:
        distname = "hMmmlplp"
    else:
        distname = "hMmmelel"
    #distname = "hMModmu2"
    #distname = "hMModmu3"
    #distname = "hMMod6"
#accidentally used the wrong name!
#distname = "Ptmmee"
compare_sigMC = False
compare_bkgMC = True

f = ROOT.TFile.Open(infile)
f.Print()

#t = f.Get("Events") 
#h = f.Get(distname)
if useLowPt:
    h = f.Get("hMmmlplp")
else:
    h = f.Get("hMmmelel")
#if year == 2223:
#    #get the 2023 data
#    h.SetName("hM2223")
#    f23 = ROOT.TFile.Open(infile23)
#    h23 = f23.Get(distname)
#    h.Add(h23)

h.Print()
#h = f.Get("hMlpe")

#factor by which to rebin the histograms
rebin = 5 #4 #6

h.Rebin(rebin) #5) #15
h.Print()

#nbins = 40 #30 #25 
#xmin = .25 #.4 #.8
#xmax = .9 #1.2
nbins = h.GetNbinsX()
xmin = h.GetXaxis().GetXmin()
xmax = h.GetXaxis().GetXmax()

#where to start and stop the fit
fitmin = .45
fitmax = .8

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


peakmin = 0.51
peakmax = 0.60

#signal shape is best described by a crystalBall.
#rrv = ROOT.RooRealVar("m_{2#mu2e}","m_{2#mu2e} [GeV]",xmin,xmax)
#rrv = ROOT.RooRealVar("m_{2#mu2e}","",xmin,xmax)
rrv = ROOT.RooRealVar("m_{2#mu2e}","",.45,xmax)
#rrv = ROOT.RooRealVar("m_{2#mu2e}","m_{2#mu2e} [GeV]",0.5,0.9)
#rrv.setRange("peak", 0.53, 0.57)
#rrv.setRange("peak", 0.51, 0.63)
rrv.setRange("peak", peakmin, peakmax)
#rrv.setRange("full", xmin, xmax)
rrv.setRange("full", fitmin, fitmax)
#try a different bkg model!!
#myfitter = fitter.fitter_4mu(mass=rrv, bkg_model='Threshold', sig_model='CB')
#myfitter = fitter.fitter_4mu(mass=rrv, bkg_model='Cheb2', sig_model='CB')
#myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb2', sig_model='CB')
#myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb2', sig_model='Voigtian')
##Voigtian with const combinatorial bkg!!
#myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb0', sig_model='Voigtian')
##Bkg model is Gaussian (resonant bkg--already fitted!) + constant !
#myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='ConstGauss', sig_model='Voigtian')

#set this to -1 for a smaller fit range (.52-.575 GeV) or +1 for a big fit range (.50-.60 GeV)
fitsize = 0
nresparam = 1 #only normalization allowed to float for resonant background
if req_elID == 2 or req_elID == 3:
    #sigMod = 'BreitWigner'
    #nsigparam = 3
    #sigMod = 'Voigtian'
    #nsigparam = 4
    sigMod = ''
    nsigparam = 0
    #sigMod = 'DoubleGauss' #nominal!!!
    #nsigparam = 5 
    #CB_Gauss
    #sigMod = 'CB_Gauss'
    #nsigparam = 7
    #comBkgMod = 'Cheb0'
    #ncomparam = 1
    comBkgMod = 'Threshold2mu2e' #new nominal!
    ncomparam = 2
    #comBkgMod = 'ThreshExp2mu2e'
    #ncomparam = 3
    #comBkgMod = 'Cheb1'
    #ncomparam = 2
    #comBkgMod = 'Cheb2'
    #ncomparam = 3
    resBkgMod = 'BreitWigner' #nominal
    #resBkgMod = 'SingleGauss' #alternative
    #resBkgMod = 'Voigtian' #other alternative
    nparam = nsigparam + ncomparam + nresparam
    myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model=comBkgMod, resBkg_model=resBkgMod, sig_model=sigMod)
elif req_elID == 1:
    #sigMod = 'BreitWigner'
    #nparam = 3 + 1 + 2
    sigMod = 'DoubleGauss'
    nparam = 5 + 1 + 2
    myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb1', resBkg_model='BreitWigner', sig_model=sigMod)
elif req_elID == 0:
    sigMod = 'TripleGauss'
    nparam = 7 + 1 + 2
    myfitter = fitter.fitter_2mu2e(mass=rrv, bkg_model='Cheb1', resBkg_model='BreitWigner', sig_model=sigMod)
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
if req_elID == 2 or req_elID == 3:
    if comBkgMod == 'Cheb0':
        myfitter.set_bkg_params( ) #a1=library.Param(0.5, .1, 1) ) 
    elif comBkgMod == 'Threshold2mu2e':
        if sigMod == '':
            myfitter.set_bkg_params( alpha=library.Param(0.01, 0.0, 10.0) )
        else:
            myfitter.set_bkg_params( alpha=library.Param(1.0, 0.5, 5.0) )
    elif comBkgMod == 'ThreshExp':
        myfitter.set_bkg_params( alpha=library.Param(1.0, 0.5, 5.0), C=library.Param(7.0, -10.0, 10.0) )
    elif comBkgMod == 'Cheb1':
        myfitter.set_bkg_params( a1=library.Param(0.5, .1, 1) ) 
    elif comBkgMod == 'Cheb2':
        myfitter.set_bkg_params( a1=library.Param(0.5, .1, 1), a2=library.Param(0.5, .1, 1) ) 
elif req_elID == 1:
    myfitter.set_bkg_params( a1=library.Param(0.5, .1, 1) ) 

#read resonant background parameters in from the txt file!
bparamfname = "fit_results/bkgMCParams_" + resBkgMod
if useLowPt:
    bparamfname += "_lowPt"
if new_wt:
    bparamfname += "_newWt%s"%(str(new_wt) if new_wt > 1 else "")
if modnum > -1:
    bparamfname += "_mod"+str(modnum)
if req_muID:
    bparamfname += "_reqMuID"
if req_elID != 2:
    bparamfname += "_req%d"%req_elID
if not do_pileup:
    bparamfname += "_noPU"
if do_trigCor:
    bparamfname += "_trigCor"
if year != 2022:
    bparamfname += "_2023"
if rebin != 5:
    bparamfname += "_%dMeVbins"%rebin
#if fitsize != 0:
#    if fitsize < 0:
#        bparamfname += "_smallRange"
#    else:
#        bparamfname += "_bigRange"
bparamfname += ".txt"
import os
if os.path.exists(bparamfname):
    bparamf = open(bparamfname, "r")
    nbkgMC = 10
    comm = "myfitter.set_resBkg_params(" 
    for i,line in enumerate(bparamf):
        words = line.split()
        pname = words[0].strip()
        pval = float(words[1])
        perr = float(words[2])
        if pname == "nbkg":
            nbkgMC = pval
        else:
            if i != 0:
                comm += ", "
            comm += "%s=library.ConstParam(%f)"%(pname, pval) 
    comm += " )"
    exec(comm)
    myfitter.set_resBkg_norm(nbkgMC, 0.2) 
else:
    print("aaaaaaaaaaaaaaaaaaaaaaaaaaa %s doesn't exist aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"%bparamfname)

    #myfitter.set_resBkg_params( mg=library.Param(.545, .544, .549), sg=library.Param(2.64e-2, 2.63e-2, 2.65e-2) )
    if req_elID == 2:
        #myfitter.set_resBkg_params( mg=library.ConstParam(.555, .554, .556), sg=library.Param(1.82e-2, 1.80e-2, 1.84e-2) )
        #myfitter.set_resBkg_params( mb=library.ConstParam(.552), wb=library.ConstParam(1.87e-2) )
        if rebin == 5:
            if cut_mee:
                #cut .04 < M_ee < .09 (kills almost all the resBkg)
                myfitter.set_resBkg_params( mb=library.ConstParam(.5750), wb=library.ConstParam(8.469e-3) )
                myfitter.set_resBkg_norm(1.061, 0.2)
            else:
                #normal
                #myfitter.set_resBkg_params( mb=library.ConstParam(.5537), wb=library.ConstParam(1.807e-2) )
                #myfitter.set_resBkg_norm(17.63, 0.2)
                #adding muon pt, eta cuts, updated event weight, pileup corex, etc
                #myfitter.set_resBkg_params( mb=library.ConstParam(.5549), wb=library.ConstParam(1.652e-2) )
                #myfitter.set_resBkg_norm(5.24, 0.2)
                if do_pileup:
                    #3866
                    if resBkgMod == 'BreitWigner':
                        myfitter.set_resBkg_params( mb=library.ConstParam(.5548), wb=library.ConstParam(1.626e-2) )
                        myfitter.set_resBkg_norm(13.35, 0.2)
                    elif resBkgMod == 'SingleGauss':
                        myfitter.set_resBkg_params( mg=library.ConstParam(.5580), sg=library.ConstParam(1.339e-02) )
                        myfitter.set_resBkg_norm(13.35, 0.2)
                    elif resBkgMod == 'Voigtian':
                        myfitter.set_resBkg_params( mv=library.ConstParam(.5588), wv=library.ConstParam(4.543e-03), sv=library.ConstParam(1.045e-02) )
                        myfitter.set_resBkg_norm(13.35, 0.2)
                else:
                    #3867
                    myfitter.set_resBkg_params( mb=library.ConstParam(.5549), wb=library.ConstParam(1.774e-2) )
                    myfitter.set_resBkg_norm(9.04, 0.2)
        elif rebin == 4:
            myfitter.set_resBkg_params( mb=library.ConstParam(.5539), wb=library.ConstParam(1.433e-2) ) 
            myfitter.set_resBkg_norm(13.37, 0.2)
        elif rebin == 6:
            myfitter.set_resBkg_params( mb=library.ConstParam(.5558), wb=library.ConstParam(1.872e-2) ) 
            myfitter.set_resBkg_norm(14.00, 0.2)
        #myfitter.set_resBkg_norm(19.82)
    elif req_elID == 1:
        #myfitter.set_resBkg_params( mv=library.ConstParam(.557), sv=library.ConstParam(1.01e-3), wv=library.ConstParam(2.99e-2) )
        #myfitter.set_resBkg_norm(199.4)
        myfitter.set_resBkg_params( mb=library.ConstParam(.5595), wb=library.ConstParam(2.346e-2) )
        #myfitter.set_resBkg_norm(199.3)
        myfitter.set_resBkg_norm(163.2, 0.2)
    elif req_elID == 3:
        #myfitter.set_resBkg_params( mb=library.ConstParam(.5540), wb=library.ConstParam(1.309e-2) )
        #myfitter.set_resBkg_norm(6.0, 0.2)
        #3871
        #myfitter.set_resBkg_params( mb=library.ConstParam(.5541), wb=library.ConstParam(1.310e-2) )
        #myfitter.set_resBkg_norm(5.91, 0.2)
        #3872--with PU corex
        myfitter.set_resBkg_params( mb=library.ConstParam(.5545), wb=library.ConstParam(1.426e-2) )
        myfitter.set_resBkg_norm(6.09, 0.2)
    elif req_elID == 0:
        ##normal
        #myfitter.set_resBkg_params( mb=library.ConstParam(0.5606), wb=library.ConstParam(1.912e-2) )
        #myfitter.set_resBkg_norm(1002.4, 0.2)
        ## .04 < M_ee < .09 cut
        #myfitter.set_resBkg_params( mb=library.ConstParam(0.5750), wb=library.ConstParam(4.897e-2) )
        #myfitter.set_resBkg_norm(131.4, 0.2)
        myfitter.set_resBkg_params( mb=library.ConstParam(0.5605), wb=library.ConstParam(1.850e-2) )
        myfitter.set_resBkg_norm(640.8, 0.2)
    
if sigMod != '':
    #read signal parameters in from the txt file!
    sparamfname = "fit_results/sigMCParams_" + sigMod
    if useLowPt:
        sparamfname += "_lowPt"
    if new_wt:
        sparamfname += "_newWt%s"%(str(new_wt) if new_wt > 1 else "")
    if "Modmu" in distname:
        sparamfname += "_modmu%d"%(int(distname[-1])) 
    elif "hMMod" in distname:
        modnum = int(distname[-1]) 
        sparamfname += "_mod"+str(modnum)
    if req_muID != 2:
        sparamfname += "_reqMuID"
    if req_elID != 2:
        sparamfname += "_req%d"%req_elID
    if not do_pileup:
        sparamfname += "_noPU"
    if do_trigCor:
        sparamfname += "_trigCor"
    if year != 2022:
        sparamfname += "_2023"
    if rebin != 5:
        sparamfname += "_%dMeVbins"%rebin
    if fitsize != 0:
        if fitsize < 0:
            sparamfname += "_smallRange"
        else:
            sparamfname += "_bigRange"
    sparamfname += ".txt"
    if os.path.exists(sparamfname):
        sparamf = open(sparamfname, "r")
        comm = "myfitter.set_sig_params(" 
        for i,line in enumerate(sparamf):
            words = line.split()
            pname = words[0].strip()
            if pname == "nsig": continue
            pval = float(words[1])
            perr = float(words[2])
            if i != 0:
                comm += ", "
            if rebin == 6:
                comm += "%s=library.Param(%f, %f, %f)"%(pname, pval, pval-2*perr, pval+2*perr) 
            else:
                if sigMod == "CB_Gauss":
                    comm += "%s=library.Param(%f, %f, %f)"%(pname, pval, pval-2*perr, pval+2*perr) 
                else:
                    comm += "%s=library.Param(%f, %f, %f)"%(pname, pval, pval-perr, pval+perr) 
        comm += " )"
        exec(comm)
    else:
        print("aaaaaaaaaaaaaaaaaaaaaaaaaaa %s doesn't exist aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"%sparamfname)
        exit()

data = ROOT.RooDataHist("data", "data", rrv, ROOT.RooFit.Import(h))

#do the fit
fitres = myfitter.model.fitTo(data, ROOT.RooFit.Save())
print("fitres: " + str(fitres)) 


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
myfitter.model.plotOn(frame, ROOT.RooFit.Name("ResBkg"), ROOT.RooFit.Components('resbkg'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(ROOT.kCyan))

#myfitter.model.plotOn(frame, ROOT.RooFit.Name("Sig"), ROOT.RooFit.Components('sig'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(4), ROOT.RooFit.LineColor(ROOT.kRed+1))
#myfitter.model.plotOn(frame, ROOT.RooFit.Name("Sig"), ROOT.RooFit.Components('CB'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(4), ROOT.RooFit.LineColor(ROOT.kRed+1))
if sigMod != '': 
    if sigMod == 'CB_Gauss':
        myfitter.model.plotOn(frame, ROOT.RooFit.Name("Sig"), ROOT.RooFit.Components('CB_Gauss'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(4), ROOT.RooFit.LineColor(ROOT.kRed+1))
    else:
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
nfitbins = int((fitmax - fitmin)/binsize)
ndf = nfitbins - nparam
chi2_SB = frame.chiSquare() #* ndf
chi2 = chi2_SB*ndf
print("***chi2/ndf of the fit: %f***"%chi2_SB) 
print("ndf = %d (nfitbins=%d, nparam=%d) --> total chi2 = %f"%(ndf, nfitbins, nparam, chi2))
#??
frame.SetName("")

c.Modified()
c.Update()

#ndata = h.Integral(h.FindBin(0.51), h.FindBin(0.63))
ndata = h.Integral(h.FindBin(peakmin), h.FindBin(peakmax))
#ndata = h.Integral(h.FindBin(0.90), h.FindBin(1.0))
argset = ROOT.RooArgSet(rrv)
if sigMod == '':
    sig_int = 0
    sig_int_full = 0
else:
    sig_int = myfitter.sig.func.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
    sig_int_full = myfitter.sig.func.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("full"))
bkg_int = myfitter.bkg.func.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
resBkg_int = myfitter.resBkg.func.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
tot_int = myfitter.model.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
if sigMod != '':
    print("sigvals: %f, %f; bkgvals: %f, %f"%(sig_int.getVal(), myfitter.nsig.getVal(), bkg_int.getVal(), myfitter.nbkg.getVal())) 
    nsig = sig_int.getVal() * myfitter.nsig.getVal()
    nsigErr = sig_int.getVal()*myfitter.nsig.getError()
    nsig_full = sig_int_full.getVal() * myfitter.nsig.getVal()
    nsigErr_full = sig_int_full.getVal()*myfitter.nsig.getError() 
    print("***Fitted signal events: %f +/- %f (stat.)***"%(nsig_full, nsigErr_full))
    print("   (Fitted signal events just under the peak: %f +/- %f)"%(nsig, nsigErr)) 
    propErr = myfitter.nsig.getPropagatedError(fitres)
    print("    Propagated error: %f ; nsig.getError(): %f"%(propErr, myfitter.nsig.getError()))  
    print("       sig_int : %f +/- %f"%(sig_int.getVal(), sig_int.getPropagatedError(fitres))) 
    print("       myfitter.nsig.getError() = %f"%(myfitter.nsig.getError())) 
    #mom2 = myfitter.nsig.moment(rrv, 2, ROOT.kFALSE, ROOT.kFALSE)
    #mom2 = myfitter.nsig.moment(rrv, 2, ROOT.kTRUE, ROOT.kFALSE)
    mom2 = myfitter.sig.func.moment(rrv, 2, ROOT.kTRUE , ROOT.kTRUE )
    #sig2_int = mom2.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("full")) 
    sig2_int = mom2.createIntegral(argset, ROOT.RooFit.Range("full")) 
    errVal = sig2_int.getVal()
    print("          2nd moment integral = %f"%(errVal)) 
else:
    nsig = 0
nbkg = bkg_int.getVal() * myfitter.nbkg.getVal()
nrbkg = resBkg_int.getVal() * myfitter.nrbkg.getVal()
if sigMod != '':
    ntot = tot_int.getVal() * (myfitter.nsig.getVal() + myfitter.nbkg.getVal() + myfitter.nrbkg.getVal() )
else:
    ntot = tot_int.getVal() * ( myfitter.nbkg.getVal() + myfitter.nrbkg.getVal() )

bkghi = bkg_int.getVal() * myfitter.nbkg.getAsymErrorHi() 
bkglo = bkg_int.getVal() * myfitter.nbkg.getAsymErrorLo() 

print("bkg errors: %f, %f "%(bkghi, bkglo)) 

#leg = ROOT.TLegend(0.4, 0.65, 0.8, 0.85)
#leg = ROOT.TLegend(0.15, 0.65, 0.45, 0.85)
leg = ROOT.TLegend(0.40, 0.50, 0.85, 0.85)
#leg.SetHeader("N: m_{2#mu2e} #in [0.51, 0.63] GeV")
leg.SetHeader("%.2f < m_{2#mu2e} < %.2f GeV"%(peakmin, peakmax))
#leg.SetHeader("N: m_{2#mu2e} #in [0.9, 1.0] GeV")
leg.SetLineWidth(0)
if sigMod != '':
    leg.AddEntry("Sig", f"Signal (CB): N = {nsig:.1f}", "l")
leg.AddEntry("Bkg", f"Comb. Background: N = {nbkg:.1f}", "l")
leg.AddEntry("ResBkg", f"Res. Background: N = {nrbkg:.1f}", "l")
leg.AddEntry("Tot", f"Sum: N = {ntot:.1f}", "l")
leg.AddEntry("Data", f"Data, N = {ndata:n}", "lep")
leg.Draw()
if req_elID == 2:
    if year == 2022:
        if rebin == 5:
            if req_muID:
                pav = ROOT.TPaveText(.46, 25, .54, 30)
                #pav = ROOT.TPaveText(.455, 22, .54, 30)
            else:
                pav = ROOT.TPaveText(.46, 30, .54, 35)
        elif rebin == 4:
            pav = ROOT.TPaveText(.46, 25, .54, 30)
        elif rebin == 6:
            pav = ROOT.TPaveText(.46, 35, .54, 40)
    elif year == 2223:
        pav = ROOT.TPaveText(.46, 40, .54, 48)
elif req_elID == 3:
    pav = ROOT.TPaveText(.46, 10, .54, 12)
elif req_elID == 1:
    pav = ROOT.TPaveText(.46, 150, .54, 180)
elif req_elID == 0:
    pav = ROOT.TPaveText(.46, 450, .54, 520)
pav.AddText("#chi^{2}/ndf = %.3f / %d "%(chi2, ndf)) 
pav.AddText("= %.3f"%(chi2_SB)) 
pav.Draw("same")
####################################
#use official CMS Style.
import utils.CMSStyle as cmsstyle
#cmsstyle.setCMSLumiStyle(c, 0, era='Run2')
if year == 2022:
    cmsstyle.setCMSLumiStyle(c, 0, era='2022')
elif year == 2223:
    cmsstyle.setCMSLumiStyle(c, 0, era='2022-2023')

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
if compare_sigMC:
    if req_elID == 2:
        #nMiss==0, vProb>.5
        #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3826.root")
        #nMiss<=3, vProb>.1
        #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3827.root")
        #WP90 only
        #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3831.root")
        if cut_mee:
            #cut .04 < M_ee < .09
            fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3851.root")
        else:
            #elpt > 2 and 2022+2023 (Nominal)
            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3837.root")
            #pileup corrections added, etc (2022 only)
            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3860.root")
            #updated weights
            if do_pileup:
                if new_wt == 1:
                    if req_muID:
                        if do_trigCor:
                            #nominal frfrfrfrfr
                            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3888.root")
                            #only one trigger path instead of 6
                            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38110.root")
                            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38130.root")
                            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38136.root")
                            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38145p4715.root")
                            if useLowPt:
                                fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38166p4734.root")
                            else:
                                #cut transition electrons, do additional corrections, nd stuff
                                #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38164p4733.root")
                                fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38174p4740.root")
                        else:
                            #medium ID
                            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3883.root")
                            #loose ID
                            fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3884.root")
                    else:
                        #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3876.root")
                        #updated PU corrections
                        #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3878.root")
                        fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3882.root")
                elif new_wt == 2:
                    fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3879.root")
                elif new_wt == 3:
                    fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3880.root")
                else:
                    fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3866.root")
            else:
                fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3867.root")
    elif req_elID == 1:
        #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3829.root")
        fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3838.root")
    elif req_elID == 3:
        #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3839.root")
        #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3871.root")
        #with PU corex
        if new_wt:
            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3875.root")
            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38122.root")
            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38138p479.root")
            fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38147p4717.root")
        else:
            fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3872.root")
    elif req_elID == 0:
        if cut_mee:
            #cut .04 < M_ee < .09
            fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3849.root")
        else:
            #normal
            #fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3850.root")
            fMC = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3858.root")
    #tMC = fMC.Get("hMe")
    tMC = fMC.Get(distname)
    #tMC = fMC.Get("hMModmu2")
    #tMC = fMC.Get("hMDnmmelel")

    if year == 2223:
        tMC.SetName("hM2223")
        print("tMC before adding 2023 MC:")
        tMC.Print()
        fMC23 = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_2023sigMCtest390.root")
        tMC23 = fMC23.Get(distname)
        print("tMC23:")
        tMC23.Print()
        tMC.Add(tMC23)
        print("new tMC:")
        tMC.Print()

    #MC_scales = [2.0, 1.0, 0.5]
    #MC_scales = [1.0]
    MC_scales = [0.40] #4]
    if req_elID == 0:
        MC_scales = [0.7]
    hMC = [None for mcs in MC_scales]
    #draw the blinded MC scaled by a few different values
    for i,scale in enumerate(MC_scales):
        #hMC[i] = ROOT.TH1F("hMC"+str(i), "#mu#muee", nbins, xmin, xmax)
        hMC[i] = tMC.Clone()
        #if incEE:
        #    hMC[i].SetName("hpostEE")
        #    #fee = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest472.root")
        #    #fee = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest474.root")
        #    fee = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest477.root")
        #    hee = fee.Get(distname)
        #    #hee = fee.Get("hMModmu2")
        #    #hee = fee.Get("hMModmu11")
        #    #hee.Add(hMC[i], 28.25/38.01) 
        #    hee.Add(hMC[i]) 
        #    hMC[i] = hee
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
        #rebin = hMC[i].GetNbinsX() / nbins
        #print("rebin: " + str(rebin))
        #irebin = int(rebin)
        hMC[i].Rebin(rebin) #irebin)
        hMC[i].Scale(scale)
        ##temporary fix: set MC lumi to just 2022 instead of 2022+23
        #hMC[i].Scale(scale*38.48/(38.48+28.89))
        hMC[i].Draw("hist same")
        #c.cd()
        #hMC.Draw("hist same")
        if i == 0:
            #mc_sig_int = hMC[i].Integral(hMC[i].FindBin(0.51), hMC[i].FindBin(0.63))
            mc_sig_int = hMC[i].Integral(hMC[i].FindBin(0.51), hMC[i].FindBin(0.60))
#if compare_MC:
#    for i,scale in enumerate(MC_scales):
        leg.AddEntry(hMC[i],"#eta#rightarrow2#mu2e MC Sig * #bar{B} * %.2f: N = %.1f"%(scale, mc_sig_int), "l")
    c.Update()

if compare_bkgMC:
    #fBkg = ROOT.TFile.Open("EtaToMuMuGamma_2018_0_ntuple_skimmedBkgMC.root")
    if req_elID == 2:
        #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3826.root")
        #vProb>.1, nMiss<=3
        #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3827.root")
        #WP90 only
        #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3831.root")
        if cut_mee:
            #cut .04 < M_ee < .09
            fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3851.root")
        else:
            #elpt > 2 and 2022+2023 (Nominal)
            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3837.root")
            #nominal, 2022 only, pileup corex, etc
            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3860.root")
            #updated weights
            if do_pileup:
                if new_wt == 1:
                    if req_muID:
                        if do_trigCor:
                            #nominal frfrfrfrfrfr
                            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3888.root")
                            #only one trigger path instead of 6
                            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38110.root")
                            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38129.root")
                            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38146.root")
                            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38145.root")
                            if useLowPt:
                                fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38167.root")
                            else:
                                #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38165.root")
                                fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38175.root")
                        else:
                            #medium ID
                            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3883.root")
                            #loose ID
                            fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3884.root")
                    else:
                        #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3876.root")
                        #updated PU corrections
                        #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3878.root")
                        fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3882.root")
                elif new_wt == 2:
                    fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3879.root")
                elif new_wt == 3:
                    fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3880.root")
                else:
                    fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3866.root")
            else:
                fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3867.root")
    elif req_elID == 1:
        #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3829.root")
        fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3838.root")
    elif req_elID == 3:
        #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3839.root")
        #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3871.root")
        #with PU corex
        if new_wt:
            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3875.root")
            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38122.root")
            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38140.root")
            fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest38148.root")
        else:
            fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3872.root")
    elif req_elID == 0:
        if cut_mee:
            #cut .04 < M_ee < .09
            fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3849.root")
        else:
            #normal
            #fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3850.root")
            fBkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3858.root")
    print("fBkg: ") 
    fBkg.Print()
    #tBkg = fBkg.Get("Events") 
    #hBkg = ROOT.TH1F("hBkg", "#mu#mu#gamma", nbins, xmin, xmax)
    #hBkg = fBkg.Get("hMe")
    hBkg = fBkg.Get(distname)
    #hBkg = fBkg.Get("hMModmu2")
    ##temporary fix: scale bkg to get back to just 2022 lumi instead of 2022-23
    #if (req_elID == 2 or req_elID == 1) and year != 2223 and not cut_mee:
    #    hBkg.Scale(38.48/(38.48+28.89)) 
    if "38145" in fBkg:
        hBkg.Scale((28.25+9.76)/28.25) 
    #rebin = hBkg.GetNbinsX() / nbins
    #print("bkg rebin: " + str(rebin))
    #irebin = int(rebin)
    hBkg.Rebin(rebin) #irebin)
    hBkg.SetLineWidth(2)
    hBkg.SetLineColor(ROOT.kViolet)
    #tBkg.Draw(distname + ">>hBkg", "Weight*(Weight>0)", "hist same")
    hBkg.Draw("hist same")
    #mc_bkg_int = hBkg.Integral(hBkg.FindBin(0.51), hBkg.FindBin(0.63))
    mc_bkg_err_db = c_double(-1.0)
    #mc_bkg_int = hBkg.IntegralAndError(hBkg.FindBin(0.51), hBkg.FindBin(0.63), mc_bkg_err_db)
    mc_bkg_int = hBkg.IntegralAndError(hBkg.FindBin(0.51), hBkg.FindBin(0.60), mc_bkg_err_db)
    print("bkg int:")
    print(mc_bkg_int)
    mc_bkg_err = float(mc_bkg_err_db.value) 
    leg.AddEntry(hBkg, "#eta#rightarrow#mu#mu#gamma MC Bkg: N = %.1f +/- %.1f"%(mc_bkg_int, mc_bkg_err), "l")


#if mc_bkg_int > 0: 
#    sigInt -= mc_bkg_int
#    bkgInt = (mc_bkg_int**2 + bkgInt**2)**0.5
#if sigMod != '' and bkgInt > 0:
#    significance = (2*((sigInt + bkgInt)*log((sigInt+bkgInt)*(bkgInt+sigmaB**2)/(bkgInt**2 + (sigInt+bkgInt)*sigmaB**2)) - bkgInt**2/sigmaB**2 *log(1 + sigmaB**2*sigInt/(bkgInt*(bkgInt+sigmaB**2)))))**0.5 
#    print("Significance is about %f sigmas."%(significance))

c.Update()

minX = .45
maxX = .8
##trying to add pull frame underneath! #####
pullframe = rrv.frame(ROOT.RooFit.Name(f"pullframe"), ROOT.RooFit.Title(" "))
pullh = frame.pullHist("Data", "Tot") #Bkg")
#print("pullh: " + str(pullh)) 
residh = frame.residHist("Data", "Tot") #Bkg")
binning = rrv.getBinning()
for i in range(1, h.GetNbinsX()+1):
    rrv.setRange("range_for_bin", binning.binLow(i), binning.binHigh(i))
    normset = ROOT.RooFit.NormSet(rrv)
    bkgPdfIntegral = myfitter.bkg().createIntegral(ROOT.RooArgSet(rrv), "range_for_bin")
    #print("bkgPdfIntegral: " + str(bkgPdfIntegral)) 
    pset = ROOT.RooArgList()
    pset.add(bkgPdfIntegral) 
    bkgYield = ROOT.RooProduct("bkgYield", "bkgYield", pset)
    fitError = bkgYield.getPropagatedError(fitres)
    if pullh.GetPointY(i) == 0:
        old_sigma_i = 1.0
    else:
        old_sigma_i = residh.GetPointY(i)/pullh.GetPointY(i)
    new_sigma_i = (old_sigma_i**2 - fitError**2)**(1/2)
    residh.SetPointY(i, residh.GetPointY(i)/new_sigma_i)
    residh.SetPointEYhigh(i, residh.GetErrorYhigh(i)/new_sigma_i)
    residh.SetPointEYlow(i, residh.GetErrorYlow(i)/new_sigma_i)
pullframe.addPlotable(residh, 'P')
pullframe.GetXaxis().SetTitle("m_{#mu#muee} (GeV)")
pullframe.GetYaxis().SetTitle("Pull") 
line = ROOT.TLine(minX, 0, maxX, 0)
line.SetLineColor(ROOT.kGray)
pullframe.addObject(line)
c1 = cmsstyle.get_fit_canvas(frame, pullframe)
pad1 = c1.FindObject("pad1")
pad1.cd()
if compare_sigMC:
    for i,scale in enumerate(MC_scales):
        hMC[i].Draw("hist same")
if compare_bkgMC:
    hBkg.Draw("hist same")
leg.Draw("same")
pav.Draw("same")
c1.Update()
############################3

if compare_sigMC and compare_bkgMC:
    print("MC integrals:\nSignal: %f\nBackground: %f"%(mc_sig_int, mc_bkg_int)) 

input("h")

#save the canvas!
canfname = "twoMu2EDataFit"
if useLowPt:
    canfname += "_lowPt"
if new_wt:
    canfname += "_newWt%s"%(str(new_wt) if new_wt > 1 else "")
if "Modmu" in distname:
    canfname += "_modmu%d"%(int(distname[-1])) 
elif "hMMod" in distname:
    canfname += "_mod%d"%(int(distname[-1])) 
if req_muID:
    canfname += "_muID"
if req_elID == 3:
    canfname += "_tightID"
elif req_elID == 0:
    canfname += "_NoElID"
if not do_pileup:
    canfname += "_NoPU"
if do_trigCor:
    canfname += "_trigCor"
if fitsize < 0:
    canfname += "_smallRange"
elif fitsize > 0:
    canfname += "_bigRange"
if rebin != 5:
    canfname += "_%dMeVbins"%rebin
canfname += "_%s_%s%sBkg"%(sigMod, comBkgMod, resBkgMod)
#c.SaveAs("AN_Figures/" + canfname + ".pdf") 
c1.SaveAs("AN_Figures/" + canfname + ".pdf") 

print("opening outf.")
#now write the results to a txt file so can look at it later
outf = open("fit_results/" + canfname + ".txt", "w") 
params = fitres.floatParsFinal()
print("params: " + str(params)) 
for i in range(len(params)):
    par = params.at(i).getVal()
    err = params.at(i).getError()
    nam = params.at(i).GetName()
    outf.write("%s\t%f\t%f\n"%(nam, par, err)) 
outf.write("***chi2/ndf of the fit: %f***\n"%chi2_SB) 
outf.write("ndf = %d (nfitbins=%d, nparam=%d) --> total chi2 = %f\n"%(ndf, nfitbins, nparam, chi2))
if sigMod != '':
    outf.write("***Fitted signal events: %f +/- %f (stat.)***\n"%(nsig_full, nsigErr_full))
outf.close()

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
