import ROOT, sys

weighted = True

year = 2022

#fit mumu peak instead of signal?
mumu = False

#include uncertainties on the event weights too (instead of just statistical)?
incWtUnct = False

#how many electrons to require elID: 0, 1, or 2 (both)? or 3 (meaning tightID WP80 req'd on both)
req_elID = 2 

#require muon ID or nah?
req_muID = True

#include pileup corrections, or nah?
do_pileup = True

#include trigger corrections, or nah?
do_trigCor = True

#use LowPt electrons instead of regular??
useLowPt = False #True

#include distinction b/t pre/post-EE?
#incEE = True

##use only 1 trigger path instead of 6?
##onetrig = True

#use the toy trigger results instead of the real ones? 0: no, 1: nominal toy, 2: down-varied 10% toy, 3: up-varied 10% toy
# 4: down-varied 5% toy, 5: up-varied 5% toy, 6: w param down-varied by 10%, 7: w param up-varied by 10%
#  8: p50 param up-varied 2%, 9: p50 param down-varied 2%, 10: p50 down 7%, 11: p50 up 7%
# -1: real MC, but trigger correction varied up by stat. uncty; -2: trigger eff correction varied down by stat. uncty
toyTrig = 0

#use new weights found with DG/Cheb4 2mu fits?
new_wt = 1 #True

#use the modified muons instead of electrons?
modmu = False
#use modified invar mass distribution (for electron efficiency calculation) ?
modnum = -1
if len(sys.argv) > 1:
    modnum = int(sys.argv[1]) 

#distname = "hMlplp"
if weighted:
    if modnum > -1:
        distname = "hMMod%s%d"%("mu" if modmu else "", modnum)
    elif mumu:
        distname = "hMmumu"
    elif useLowPt:
        distname = "hMmmlplp"
    else:
        distname = "hMmmelel"
        if incWtUnct:
            distnameUp = "hMUpmmelel"
            distnameDn = "hMDnmmelel"
else:
    distname = "hMNoWtmmelel"

if year == 2022:
    #sigfile = "bparking_sigMCtest30.root"
    #sigfile = "bparking_sigMCtest31.root"
    #sigfile = "bparking_sigMCtest374.root"
    if mumu:
        if toyTrig == 0:
            #nominal mumu
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3888.root"
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38130.root"
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38136.root"
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38137p478.root"
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38145p4715.root"
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38171p4737.root"
            #if incEE:
            #    #sigEEfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest472.root"
            #    #sigEEfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38131.root"
            #    sigEEfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest477.root"
        elif toyTrig == 1:
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3890.root"
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38111.root"
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38113.root"
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38114.root"
        elif toyTrig == 2:
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3892.root"
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38116.root"
        elif toyTrig == 3:
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3891.root"
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38115.root"
        elif toyTrig == 4:
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3894.root"
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38118.root"
        elif toyTrig == 5:
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3893.root"
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38117.root"
        elif toyTrig == 6:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3896.root"
        elif toyTrig == 7:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3895.root"
        elif toyTrig == 8:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3898.root"
        elif toyTrig == 9:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3897.root"
        elif toyTrig == 10:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38100.root"
        elif toyTrig == 11:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3899.root"
        elif toyTrig == -1:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38101.root"
        elif toyTrig == -2:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest38102.root"
    elif req_elID == 2:
        #nMiss==0, vProb>.5
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3826.root"
        #nMiss<=3, vProb>.1
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3827.root"
        #also WP90 only
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3831.root"
        ##electron pt > 2 GeV
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3837.root"
        ##singleVertex only
        ##sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3840.root"
        ##sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3844.root"
        ##cutting out .04 < M_ee < .09
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3851.root"
        #muon pt,eta cuts too, and pileup corrections
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3860.root"
        #updated weights, nd stuff
        if do_pileup:
            if new_wt == 1:
                if req_muID:
                    if do_trigCor:
                        if toyTrig == -1:
                            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38106.root"
                        else:
                            if toyTrig == 0:
                                #nominal frfrfrfrfr
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3888.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38103.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38110.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38126.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38130.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38136.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38144p4714.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38145p4715.root"
                                if useLowPt:
                                    sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38166p4734.root"
                                else:
                                    #cut electrons from transition b/t ECAL barrel and EC
                                    #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38164p4733.root"
                                    #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38172p4738.root"
                                    #with SFs
                                    sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38174p4740.root"
                                #if incEE:
                                #    #sigEEfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest472.root"
                                #    #sigEEfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38131.root"
                                #    sigEEfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest477.root"
                                #updated electron Reco eff uncty calculation
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38123.root"
                            elif toyTrig == 1:
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3890.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38111.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38112.root"
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38113.root"
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38114.root"
                            elif toyTrig == 2:
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3892.root"
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38116.root"
                            elif toyTrig == 3:
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3891.root"
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38115.root"
                            elif toyTrig == 4:
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3894.root"
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38118.root"
                            elif toyTrig == 5:
                                #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3893.root"
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38117.root"
                            elif toyTrig == 6:
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3896.root"
                            elif toyTrig == 7:
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3895.root"
                            elif toyTrig == 8:
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3898.root"
                            elif toyTrig == 9:
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3897.root"
                            elif toyTrig == 10:
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38100.root"
                            elif toyTrig == 11:
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3899.root"
                            elif toyTrig == -1:
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38101.root"
                            elif toyTrig == -2:
                                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38102.root"
                    else:
                        #medium ID
                        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3883.root"
                        #loose ID
                        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3884.root"
                        #TM trig path only
                        sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38107.root"
                else:
                    #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3876.root"
                    #updated new wt pileup corrections -- will it make any diff?? prolly not...
                    #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3878.root"
                    sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3882.root"
            elif new_wt == 2:
                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3879.root"
            elif new_wt == 3:
                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3880.root"
            else:
                sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3866.root"
        else:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3867.root"
    elif req_elID == 1:
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3829.root"
        #WP90 only!
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3834.root"
        sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3838.root"
    elif req_elID == 3:
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3839.root"
        #new weights, cuts, nd stuff
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3871.root"
        #PU corex
        if new_wt:
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3875.root"
            #one trigger path only, medium muon ID
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38122.root"
            #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38138p479.root"
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38147p4717.root"
        else:
            sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3872.root"
    elif req_elID == 0:
        ##normal
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3850.root"
        sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3858.root"
        #cut out .04 < M_ee < .09
        #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3849.root" 
    else:
        print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
        exit()
elif year == 2023:
    #sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2023/ultraskimmed/bparking_2023sigMCtest432.root"
    #no pileup correction
    sigfile = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2023/ultraskimmed/bparking_2023sigMCtest433.root"
print("Opening sigfile %s"%sigfile) 
f = ROOT.TFile.Open(sigfile)
h = f.Get(distname)
#if incEE:
#    h.SetName(distname+"postEE")
#    fEE = ROOT.TFile.Open(sigEEfile)
#    hEE = fEE.Get(distname)
#    #hEE.Add(h, 28.9/38.0)
#    hEE.Add(h)
#    h = hEE
###temporary to get back to 2022 lumi only
##if req_elID == 1:
##    h.Scale(38.48/(38.48+28.89)) 

#try a few different rebin factors to check its effect on N2mu2e
rebin = 5 #6 #4 #
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

fitmin = .505 #.50  #.520 #
fitmax = .585 #.60  #.575 #

if mumu:
    #smod = "TripleGauss"
    smod = "DoubleGauss"
else:
    smod = "DoubleGauss"
    #smod = "Voigtian"
    #smod = "CB_Gauss"

#rrv = ROOT.RooRealVar("m_{2#mu2e}", "", .48, .7) #xmin, xmax) 
#rrv = ROOT.RooRealVar("m_{2#mu2e}", "", .48, .62) #xmin, xmax) 
rrv = ROOT.RooRealVar("m_{2#mu2e}", "", fitmin, fitmax) #xmin, xmax) 
#rrv.setRange('peak', 0.51, 0.63)
rrv.setRange('peak', fitmin, fitmax)
rrv.setRange("full", xmin, xmax)

sys.path.append("tm_analysis/analysis/python/")
import utils.fit_function_library as library

if mumu:
    sigModel = library.get_fit_function(smod, rrv)
    if smod == "DoubleGauss":
        nparam = 5
        if modmu and modnum > -1 and modnum < 6:
            sigModel.set_params( mg=library.Param(.533, .518, .548), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
        elif modmu and modnum > 5:
            sigModel.set_params( mg=library.Param(.553, .543, .588), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
        elif (not modmu) and modnum > -1 and modnum < 6:
            print("blblblblblblbblblblblblblblblblbl")
            sigModel.set_params( mg=library.Param(.545, .528, .553), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
        else:
            sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
elif req_elID == 2 or req_elID == 3:
    ##Breit-Wigner?
    #sigModel = library.get_fit_function('BreitWigner', rrv)
    #nparam = 3
    #sigModel.set_params( mb=library.Param(.548, .543, .553), wb=library.Param(.005, .0005, .05) )
    #DoubleGaussian??
    sigModel = library.get_fit_function(smod, rrv)
    if smod == "DoubleGauss":
        nparam = 5
        if modmu and modnum > -1 and modnum < 6 or (not modmu) and modnum > 0 and modnum < 6:
            sigModel.set_params( mg=library.Param(.538, .533, .553), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
        elif modmu and modnum > 5:
            sigModel.set_params( mg=library.Param(.553, .543, .558), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
        else:
            sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
    #CrystalBall + Gaussian??
    #sigModel = library.get_fit_function('CB_Gauss', rrv)
    elif smod == "CB_Gauss":
        nparam = 7
        sigModel.set_params( mcb=library.Param(.548, .543, .553), acb=library.Param(-13.7, -20.0, -1.0), ncb=library.Param(200, 10, 500), scb=library.Param(.010, .005, .025), sg=library.Param(0.01, 0.001, 0.1), CB_frac=library.Param(0.5, 0.01, 0.99) )
    #if do_pileup:
    #    if rebin == 5 or rebin == 4:
    #        sigModel.set_params( mcb=library.Param(.548, .543, .553), acb=library.Param(-13.7, -20.0, -1.0), ncb=library.Param(200, 10, 500), scb=library.Param(.010, .005, .015), sg=library.Param(0.01, 0.001, 0.1), CB_frac=library.Param(0.5, 0.01, 0.99) )
    #    elif rebin == 6:
    #        sigModel.set_params( mcb=library.Param(.548, .543, .553), acb=library.Param(-13.7, -20.0, -1.0), ncb=library.Param(200, 10, 500), scb=library.Param(.010, .008, .012), sg=library.Param(0.07, 0.005, 0.15), CB_frac=library.Param(0.5, 0.01, 0.99) )
    #else:
    #    sigModel.set_params( mcb=library.Param(.548, .543, .553), acb=library.Param(-13.7, -20.0, -1.0), ncb=library.Param(200, 10, 500), scb=library.Param(.010, .001, .012), sg=library.Param(0.01, 0.001, 0.1), CB_frac=library.Param(0.5, 0.01, 0.99) )
    ##TripleGaussian??
    #sigModel = library.get_fit_function('TripleGauss', rrv)
    #nparam = 7
    #sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.03, .001, .8), sg2=library.Param(.007, .001, .7), sg3=library.Param(.001, .00001, .5), sig1frac=library.Param(0.8, 0.1, 0.9), sig2frac=library.Param(0.2, 0.01, 0.95) )
    ##Voigtian?
    #sigModel = library.get_fit_function('Voigtian', rrv)
    elif smod == "Voigtian":
        nparam = 4
        sigModel.set_params( mv=library.Param(.548, .543, .553), wv=library.Param(.005, .0005, .05), sv=library.Param(.005, .0000001, .05) )
    elif smod == "TripleGauss":
        nparam = 7
        sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.02, .001, .1), sg2=library.Param(.008, .001, .05), sg3=library.Param(.005, .001, .1), sig1frac=library.Param(0.8, 0.05, 0.95), sig2frac=library.Param(0.05, 0.01, 0.99) )
        #sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.03, .001, .8), sg2=library.Param(.007, .001, .7), sg3=library.Param(.001, .00001, .5), sig1frac=library.Param(0.8, 0.1, 0.9), sig2frac=library.Param(0.2, 0.01, 0.95) )
    ##CrystalBall function?
    #sigModel = library.get_fit_function('CB', rrv)
    #nparam = 5
    #sigModel.set_params( mcb=library.Param(.548, .543, .553), acb=library.Param(-1.5, -100, -0.01), ncb=library.Param(200, 10, 5000), scb=library.Param(.010, .005, .05) )
elif req_elID == 1:
    #DoubleGaussian??
    sigModel = library.get_fit_function('DoubleGauss', rrv)
    nparam = 5
    sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.01, .001, .1), sg2=library.Param(.005, .001, .05), sig1frac=library.Param(0.8, 0.1, 0.9) )
elif req_elID == 0:
    #TripleGaussian??
    sigModel = library.get_fit_function('TripleGauss', rrv)
    nparam = 7
    sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.03, .001, .8), sg2=library.Param(.007, .001, .7), sg3=library.Param(.001, .00001, .5), sig1frac=library.Param(0.8, 0.1, 0.9), sig2frac=library.Param(0.2, 0.01, 0.95) )

##sigModel.set_params( mcb=library.Param(.548, .543, .553), acb=library.Param(-4.5, -6, -0.1), ncb=library.Param(40, .1, 400), scb=library.Param(.0195, .01, .05) )
##Voigtian?
#sigModel = library.get_fit_function('Voigtian', rrv)
#nparam = 4
#sigModel.set_params( mv=library.Param(.548, .543, .553), wv=library.Param(.005, .0005, .05), sv=library.Param(.005, .0000001, .05) )
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
##TripleGaussian??
#sigModel = library.get_fit_function('TripleGauss', rrv)
#nparam = 7
#sigModel.set_params( mg=library.Param(.548, .543, .553), sg1=library.Param(.03, .001, .5), sg2=library.Param(.007, .001, .5), sg3=library.Param(.001, .0001, .01), sig1frac=library.Param(0.8, 0.1, 0.9), sig2frac=library.Param(0.2, 0.05, 0.95) )


sigMC = ROOT.RooDataHist("sigMC", "sigMC", rrv, ROOT.RooFit.Import(h)) 
if mumu:
    nsig = ROOT.RooRealVar("nsig", "nsig", 50000, 1, 10000000)
else:
    nsig = ROOT.RooRealVar("nsig", "nsig", 50, 1, 10000)

model = ROOT.RooExtendPdf("model", "extended pdf", sigModel(), nsig)
fitres = model.fitTo(sigMC, ROOT.RooFit.Save())
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

#save the fit parameters to a text file (to be read in by fitInvariantMass_ROOFit.py
params = fitres.floatParsFinal()
#errors = fitres.GetParErrors()
outfname = "fit_results/%sMCParams_"%("sig" if not mumu else "mumu") + smod
if useLowPt:
    outfname += "_lowPt"
if new_wt:
    outfname += "_newWt%s"%(str(new_wt) if new_wt > 1 else "")
if toyTrig != 0:
    outfname += "_toyTrig%d"%toyTrig
if modnum > -1:
    outfname += "_mod%s"%("mu" if modmu else "")+str(modnum)
if req_muID:
    outfname += "_reqMuID"
if req_elID != 2:
    outfname += "_req%d"%req_elID
if not do_pileup:
    outfname += "_noPU"
if do_trigCor:
    outfname += "_trigCor"
if year != 2022:
    outfname += "_2023"
if rebin != 5:
    outfname += "_%dMeVbins"%rebin
if fitmin != .505:
    if fitmin > .505:
        outfname += "_smallRange"
    else:
        outfname += "_bigRange"
outfname += ".txt"
outf = open(outfname, "w") 
for i in range(len(params)):
    par = params.at(i).getVal()
    err = params.at(i).getError()
    nam = params.at(i).GetName()
    outf.write("%s\t%f\t%f\n"%(nam, par, err)) 
print("%s written."%outfname) 
outf.close()

#save the canvas
cname = "AN_Figures/MC_%s_"%("EtaTo2Mu2E" if not mumu else "EtaTo2Mu")
if useLowPt:
    cname += "lowPt_"
if new_wt:
    cname += "newWt%s_"%(str(new_wt) if new_wt > 1 else "")
if toyTrig != 0:
    cname += "toyTrig%d_"%toyTrig
if modnum > -1:
    cname += "mod%s"%("mu" if modmu else "")+str(modnum)+"_"
if req_muID:
    cname += "muID_"
if req_elID == 3:
    cname += "tightId_"
elif req_elID == 0:
    cname += "NoElID_"
if not do_pileup:
    cname += "NoPU_"
if do_trigCor:
    cname += "trigCor_"
if rebin != 5:  
    cname += "%dMeVbins_"%rebin
if fitmin > .505:
    cname += "smallRange_"
elif fitmin < .505:
    cname += "bigRange_"
cname += smod
cname += ".pdf"
c.SaveAs(cname)
