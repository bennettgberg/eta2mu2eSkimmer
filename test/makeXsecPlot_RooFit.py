import ROOT
from array import array
import sys
sys.path.append("tm_analysis/analysis/python/")
#import fitter
import utils.fitter as fitter

#if only testing, use smaller input file and do NOT overwrite output files
test = False #True
#set to true if want to remake all the plots for the individual pT bins
remake = True #False

#include the Run2 Scouting result on the same plot for comparison? or nah
inc_run2 = True #False

#include the corrections, or JUST show the uncorrected xsec ?
inc_corr = True

#make plot of xsec as a function of pT for eta->2mu!
if test:
    #f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/bparking_datatest335_1C.root", "read")
    #f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/bparking_datatest335_0.root", "read")
    #f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/bparking_datatest3613_1C.root", "read")
    f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/bparking_datatest3615_1C.root", "read")
else:
    #338: tight muID !!
    #f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest338_ALL.root", "read")
    #f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest376_ALL.root", "read")
    f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3818_ALL.root", "read")

h2d = f.Get("hMvsPtmumu")

#h2d.Draw("colz")
#input("Continue?")

ptmin = h2d.GetXaxis().GetXmin()
ptmax = h2d.GetXaxis().GetXmax()
nptbins = h2d.GetXaxis().GetNbins()

mmin = h2d.GetYaxis().GetXmin()
mmax = h2d.GetYaxis().GetXmax()
nMbins = h2d.GetYaxis().GetNbins()

etaMlo = .52
etaMhi = .58

lumi = 38.48 #fb^-1
bratio = 5.8e-6

#lumi in Run2 Scouting (needed for comparison)
oldlumi = 101

#first get Run2 result to be shown as reference
fRun2 = ROOT.TFile.Open("xsecs.root")
hXsecRun2 = fRun2.Get("corr_xsec")
hXsecRun2Unc = fRun2.Get("uncorr_xsec")

#make pT bins for new histogram (at low pT, width is 1 GeV; then gradually widen as pT goes up)
newbins = [i*1.0 for i in range(5, 31)]
for i in range(32, 41, 2):
    newbins.append(1.0*i)
for i in range(45, 56, 5):
    newbins.append(1.0*i)
newbins.append(70.0)
newbins.append(100.0)
newbins.append(101.0)
print("new pT bins: " + str(newbins)) 
#last entry in newbins is the upper limit of the last bin, not a new bin
newnptbins = len(newbins)-1
### new bins ^^

#construct the 1d hist of xsec vs pT -- UNCORRECTED
#hXsec = ROOT.TH1F("hXsecUncor", "hXsecUncor", nptbins, ptmin, ptmax)
hXsec = ROOT.TH1F("hXsecUncor", "hXsecUncor", newnptbins, array('d', newbins))
hXsec.Sumw2()
print("hXsec:") 
hXsec.Print()
#Xsec corrected for acceptance
#hXsecCor = ROOT.TH1F("hXsecCor", "hXsecCor", nptbins, ptmin, ptmax)
hXsecCor = ROOT.TH1F("hXsecCor", "hXsecCor", newnptbins, array('d', newbins))
hXsecCor.Sumw2()

#get efficiency histograms for correction purposes
#f2mu = ROOT.TFile.Open("bparking_mumuMCtest338.root")
#f2mu = ROOT.TFile.Open("bparking_mumuMCtest377.root")
#f2mu = ROOT.TFile.Open("bparking_mumuMCtest385.root")
f2mu = ROOT.TFile.Open("bparking_mumuMCtest388.root")
ptGen = f2mu.Get("hpTGenAll")
ptGen.SetName("hpTGenAllmumu")
ptGen.Rebin(5)
#ptAcc = f2mu.Get("hpTGenAccmumu")
ptAcc = f2mu.Get("hpTGenAccdRmumu")
#bins should be 1 GeV wide!!!
ptAcc.Rebin(5)
ptAcc.Sumw2()
ptAcc.Divide(ptGen)

#get efficiency histograms for correction purposes
f2mu2e = ROOT.TFile.Open("bparking_sigMCtest335.root")
ptGen2mu2e = f2mu2e.Get("hpTGenAll")
ptGen2mu2e.Rebin(5)
ptAcc2mu2e = f2mu2e.Get("hpTGenAccmmelel")
#ptAcc2mu2e = f2mu2e.Get("hpTGenAccmmlplp")
#bins should be 1 GeV wide!!!
ptAcc2mu2e.Rebin(5)
ptAcc2mu2e.Sumw2()
ptAcc2mu2e.Divide(ptGen2mu2e)

#counter for which newbin we're on
nb = 0
nextbin = newbins[nb]

if not test:
    fout = open("xsec_compare.txt", "w")
    fout.write("pT bin\tBParking Fitted Yield\tBParking Acc\tBParking XS\t\t\tRun2 Scouting Yield\tRun2 Scouting Acc\tRun2 Scouting XS\n")

#total number of eta mesons expected
totEta = 0.0
#sum of squares of uncertainty on total number of eta mesons expected (will take sqrt at end)
errEta = 0.0

import utils.fit_function_library as library

lowbin = -1
lowpt = -1
hipt = nextbin-1
hibin = int(hipt)-1
#need efficiency to correct!!
for j in range(nptbins+1): #-1):
    i = j+1
    #print("******i=%d, j=%d*********"%(i, j)) 
    #first few bins have acceptance of 0, so skip them altogether
    #if i < 7: # or i > 16: 
    #    #hXsec.SetBinContent(i, 0)
    #    #hXsec.SetBinError(i, 0)
    #    continue
    #skip any old pT bins that aren't the start of the new bin (these were used already).
    #if i != nextbin:
    if h2d.GetXaxis().GetBinLowEdge(i) != nextbin:
        print("i=%d; %f < pT < %f; continuing."%(i, h2d.GetXaxis().GetBinLowEdge(i), h2d.GetXaxis().GetBinLowEdge(i)+h2d.GetXaxis().GetBinWidth(i)))
        continue 
    #print("i=%d; %f < pT < %f; reached nextbin!"%(i, h2d.GetXaxis().GetBinLowEdge(i), h2d.GetXaxis().GetBinLowEdge(i)+h2d.GetXaxis().GetBinWidth(i)))
    #do fit to invariant mass spectrum within this pT bin!!
    lowpt = hipt
    lowbin = hibin+1
    hibin = i-1
    hipt = nextbin
    nextbin = newbins[nb+1]
    print("*********lowbin: %d, nextbin: %d, lowpt: %f, hipt: %f******************"%(lowbin, nextbin, lowpt, hipt)) 
    if remake and lowpt > 69: #28: #40: #3: #i > 30 and i < 32: 
        ctest = ROOT.TCanvas()
        ctest.cd()
    #projection = h2d.ProjectionY("hx"+str(i),i+1, i+1)
    #print("nb=%d: old pt bins %d thru %d"%(nb, lowbin, hibin)) 
    projection = h2d.ProjectionY("hx"+str(nb), lowbin, hibin)
    #factor by which to rebin
    if lowpt > 69:
        print("rebin factor is 80 :DDDDD")
        rebin = 4 #80
    #elif lowpt < 7 or lowpt > 37: #i > 31:
    elif lowpt < 7 or lowpt > 33: #i > 31:
        rebin = 2 #40
    else:
        rebin = 1 #20
    if test:
        rebin = 25 #*= 2
    projection.Rebin(rebin)
    #print("projection: " + str(i))
    projection.Print()
    #now do a fit to pol2 + gaussian
    #fit_func = ROOT.TF1("fit_func"+str(i), "[0] + [1]*x + [2]*x*x + [3]*exp(-0.5*((x-[4])/[5])**2)", 0.45, 0.65)
    #fit_func = ROOT.TF1("fit_func"+str(nb), "[0] + [1]*x + [2]*x*x + [3]*exp(-0.5*((x-[4])/[5])**2)", 0.45, 0.65)
    #fit_str = "[0] + [1]*x + [2]*x*x + [6]*x*x*x + [3]*exp(-0.5*((x-[4])/[5])**2)"
    #fit_func = ROOT.TF1("fit_func"+str(nb), fit_str, 0.45, 0.65)
    #print("Now i = %d. Low edge: %f, high edge: %f"%(i, lowpt, hipt)) 
    
    mlo = .45
    mhi = .65
    mbinwidth = (projection.GetXaxis().GetXmax() - projection.GetXaxis().GetXmin()) / projection.GetNbinsX()
    nmbins = int((mhi-mlo)/mbinwidth) 
    rrv = ROOT.RooRealVar("m_{2#mu}(p_{T}bin"+str(j)+")","",mlo,mhi)
    data = ROOT.RooDataHist("data"+str(j), "data"+str(j), rrv, ROOT.RooFit.Import(projection))
    #total number of fit parameters
    sigModel = 'Voigtian' #'CB' #'DoubleGauss' #'SingleGauss'
    if sigModel == 'SingleGauss':
        nparams = 6
    elif sigModel == 'DoubleGauss':
        nparams = 8
    elif sigModel == 'CB':
        nparams = 8
    elif sigModel == 'Voigtian':
        nparams = 7
    myfitter = fitter.fitter_2mu(mass=rrv, bkg_model='Cheb3', sig_model=sigModel) 
    if sigModel == 'SingleGauss':
        if lowpt > 41:
            myfitter.set_sig_params( mg=library.Param(.547, .544, .548), sg=library.Param(.01, .0001, .02) )
        elif lowpt > 27:
            myfitter.set_sig_params( mg=library.Param(.546, .547, .550), sg=library.Param(.01, .0001, .02) )
        elif lowpt > 22:
            myfitter.set_sig_params( mg=library.Param(.548, .545, .551), sg=library.Param(.005, .0001, .01) )
    elif sigModel == 'DoubleGauss':
        pass
    elif sigModel == 'CB':
        pass
    elif sigModel == 'Voigtian':
        pass
    myfitter.set_bkg_params( a1=library.Param(-.63, -1, 1), a2=library.Param(.92, -1., 1.), a3=library.Param(.001, -1., 1.) )
    fit_result = myfitter.model.fitTo(data, ROOT.RooFit.Save()) 
    #if test:
    #    fit_func.SetParameters(14000, .1, -.1, 1000, 0.548, 0.01, .1)
    ##elif nb == 1:
    ##    fit_func.SetParameters(400, .1, -.1, 50, 0.548, 0.01)
    #elif lowpt < 7:
    #    fit_func.SetParameters(20000, .1, -.1, 5000, 0.548, 0.01)
    #elif lowpt < 14:
    #    fit_func.SetParameters(100000, .1, -.1, 5000, 0.548, 0.01)
    #elif lowpt < 18:
    #    fit_func.SetParameters(25000, .1, -.1, 8000, 0.548, 0.01)
    #elif lowpt < 20:
    #    fit_func.SetParameters(10000, .1, -.1, 3000, 0.548, 0.01)
    #elif lowpt < 24:
    #    fit_func.SetParameters(6000, .1, -.1, 2000, 0.548, 0.01)
    #elif lowpt < 28:
    #    fit_func.SetParameters(3000, .1, -.1, 1000, 0.548, 0.01)
    #elif lowpt < 29:
    #    fit_func.SetParameters(2000, .1, -.1, 400, 0.548, 0.01)
    #elif lowpt < 31:
    #    fit_func.SetParameters(5000, .1, -.1, 600, 0.548, 0.01)
    #elif lowpt < 37:
    #    fit_func.SetParameters(4400, .1, -.1, 500, 0.548, 0.01)
    #elif lowpt < 38:
    #    fit_func.SetParameters(2000, .1, -.1, 200, 0.548, 0.01)
    #elif lowpt < 44:
    #    fit_func.SetParameters(3600, .1, -.1, 400, 0.548, 0.01)
    #else:
    #    fit_func.SetParameters(2000, .1, -.1, 250, 0.548, 0.01)
    #fit_func.SetParLimits(3, 0, 1000000)
    #fit_func.SetParLimits(4, .54, .56)
    #fit_func.SetParLimits(5, 0, .1)
    #fit_result = projection.Fit(fit_func, "RS")
    #now get the fit results, integrate to get the number of eta mesons
    #params = []
    #errors = []
    #paramList = fit_result.floatParsFinal()
    #for idx in range(nparams):
    #    par = paramList.at(idx)
    #    params.append( par.getVal() )
    #    errors.append( par.getError() )
    #print("params: %s"%str(params)) 
    #print("errors: %s"%str(errors)) 
    #s2pi = (2*3.14159265359)**0.5
    #nEta = params[3] * params[5] * s2pi / projection.GetBinWidth(1)
    ptbinwidth = hipt - lowpt
    #unct = ( (errors[3] / params[3])**2 + (errors[5] / params[5] )**2 )**0.5 * nEta
    #print("nEta = %f +/- %f ; ptbinwidth = %f"%(nEta, unct, ptbinwidth)) 
    #if test and i == 11 :
    ndata = projection.Integral(projection.FindBin(0.52), projection.FindBin(0.58))
    #ndata = h.Integral(h.FindBin(0.90), h.FindBin(1.0))
    argset = ROOT.RooArgSet(rrv)
    sig_int = myfitter.sig.func.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
    bkg_int = myfitter.bkg.func.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
    tot_int = myfitter.model.createIntegral(argset, ROOT.RooFit.NormSet(argset), ROOT.RooFit.Range("peak"))
    print("sigvals: %f, %f; bkgvals: %f, %f"%(sig_int.getVal(), myfitter.nsig.getVal(), bkg_int.getVal(), myfitter.nbkg.getVal())) 
    nsig = sig_int.getVal() * myfitter.nsig.getVal()
    nbkg = bkg_int.getVal() * myfitter.nbkg.getVal()
    ntot = tot_int.getVal() * (myfitter.nsig.getVal() + myfitter.nbkg.getVal())

    nEta = nsig
    unct = sig_int.getVal() * myfitter.nsig.getError()
    print("nEta = %f +/- %f ; ptbinwidth = %f"%(nEta, unct, ptbinwidth)) 

    if remake and lowpt > 69: #28: # 3: #30 and i < 32:
        frame = rrv.frame(ROOT.RooFit.Name(f"pullframe"), ROOT.RooFit.Title(" "))
        frame.SetTitle("Two-Muon Invariant Mass spectrum, %.1f < p_{T} < %.1f GeV"%(lowpt, hipt))
        frame.GetXaxis().SetTitle("m_{2#mu} [GeV]")
        binsize = mbinwidth
        frame.GetYaxis().SetTitle("Events / (%.3f GeV)"%(binsize))
        data.plotOn(frame, ROOT.RooFit.Name("Data"), ROOT.RooFit.DrawOption("PEZ"))
        myfitter.model.plotOn(frame, ROOT.RooFit.Name("Bkg"), ROOT.RooFit.Components('bkg'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(ROOT.kGreen-1))
        #myfitter.model.plotOn(frame, ROOT.RooFit.Name("Sig"), ROOT.RooFit.Components('sig'), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineStyle(4), ROOT.RooFit.LineColor(ROOT.kRed+1))
        myfitter.model.plotOn(frame, ROOT.RooFit.Name("Tot"), ROOT.RooFit.LineWidth(5), ROOT.RooFit.LineColor(ROOT.kBlue))
        datamin = 9999999 # data.GetMinimum()
        datahist = data.createHistogram("m_{2#mu}(p_{T}bin"+str(j)+")")
        #print("%d bins in the data hist!"%nmbins) 
        for ijk in range(1, nmbins+1):
            thisbin = datahist.GetBinContent(ijk)
            if thisbin < datamin:
                datamin = thisbin
        rchi2 = frame.chiSquare()
        ndf = nmbins - nparams
        chi2 = rchi2*ndf
        #datamax = 0.0
        #data.getRange(data, datamin, datamax)
        frame.SetMinimum( datamin / 1.2 ) 
        #print("trying to plot!!!!!!!!!!!!")
        #projection.SetMarkerStyle(20)
        #projection.Draw("PE")
        ##projection.GetXaxis().SetRangeUser(.45, .65) 
        #projection.GetXaxis().SetRangeUser(.5, .6) 
        #projection.GetXaxis().SetTitle("Invariant mass (GeV)")
        #projection.GetYaxis().SetTitle("Events / .0025 GeV")
        #projection.SetTitle("Invariant mass spectrum for %.1f < p_{T} < %.1f"%(lowpt, hipt)) 
        #leg = ROOT.TLegend()
        #header = fit_str + "\n"
        #for p in range(6):
        #    header += "Param %d: %.3f +/- %.3f"%(p, params[p], errors[p])
        #leg.SetHeader(header) 
        #leg.Draw("same") 
        ##print("nEta = %f"%nEta)
        #ctest.Modified()
        #ctest.Update()
        frame.Draw("AC")
        pav = ROOT.TPaveText(.49,datamin/1.2,.61,datamin/1.02, "NB")
        #pav = ROOT.TPaveText()
        pav.AddText("Fit Parameters:")
        #print("All sig pars: ") # + str(myfitter.sig.get_arg_list())) 
        for par in myfitter.sig.get_arg_list():
            pav.AddText("%s: %.5f +/- %.5f"%(par.GetName(), par.getValV(), par.getError())) 
            #print("%s: %.5f +/- %.5f"%(par.GetName(), par.getValV(), par.getError())) 
        print("All bkg pars: ") # + str(myfitter.sig.get_arg_list())) 
        for par in myfitter.bkg.get_arg_list():
            pav.AddText("%s: %.5f +/- %.5f"%(par.GetName(), par.getValV(), par.getError())) 
            #print("%s: %.5f +/- %.5f"%(par.GetName(), par.getValV(), par.getError())) 
        pav.AddText("#chi^{2}/ndf = %.3f / %d = %.3f"%(chi2, ndf, rchi2)) 
        pav.Draw("same")
        leg = ROOT.TLegend(0.55, 0.65, 0.95, 0.9)
        #leg.SetHeader("N: m_{2#mu2e} #in [0.51, 0.63] GeV")
        leg.SetHeader(".51 < m_{2#mu2e} < 0.63 GeV")
        #leg.SetHeader("N: m_{2#mu2e} #in [0.9, 1.0] GeV")
        leg.SetLineWidth(0)
        #leg.AddEntry("Sig", f"Signal ({sigModel:s}): N = {nsig:.1f}", "l")
        leg.AddEntry("Bkg", f"Background: N = {nbkg:.1f}", "l")
        leg.AddEntry("Tot", f"Sum: N = {ntot:.1f}", "l")
        leg.AddEntry("Data", f"Data, N = {ndata:n}", "lep")
        leg.Draw("same")
        ctest.Modified()
        ctest.Update()
        if not test:
            ctest.SaveAs("fit_plots_%s/bin%d.pdf"%(sigModel, lowbin-1)) 
        if lowpt > 5:
            input("wait...")
    #if nEta < 0:
    #    print("nEta < 0 ????????????")
    #    projection.Draw("hist")
    #    input("continue?")
    if unct > nEta:
        #hXsec.SetBinContent(nb, 0)
        #hXsec.SetBinError(nb, 0)
        ##projection.Draw("hist")
        ##input("continue?")
        #nb += 1
        #continue
        pass
    #print("old pt bins %d to %d; setting bin content for bin %d"%(i, nextbin, nb))
    ##nEta = 0
    ##for j in range(nMbins):
    ##    center = h2d.GetYaxis().GetBinCenter(j)
    ##    if center >= etaMlo and center <= etaMhi:
    ##        nEta += h2d.GetBinContent(i, j)
    scale = 1.0 / (lumi * bratio * ptbinwidth)
    xsec = nEta*scale
    #hXsec.SetBinContent(i, xsec)
    hXsec.SetBinContent(nb, xsec)
    xerr = unct*scale
    #hXsec.SetBinError(i, xerr)
    hXsec.SetBinError(nb, xerr)
    print("hXsec: setting bin %d content to %f, error to %f"%(nb, xsec, xerr)) 
    hXsec.Print()

    #now do the corrected one
    acc = ptAcc.GetBinContent(j)
    accErr = ptAcc.GetBinError(j)

    print("Acceptance histogram j=%d, lowpt: %f, binwidth: %f"%(j, ptAcc.GetBinLowEdge(j), ptAcc.GetBinWidth(j))) 
    
    if acc > 0.:
        Ncor = nEta / acc
        Nerr = Ncor * ( (unct / nEta)**2 + (accErr / acc)**2 )**0.5
    else:
        Ncor = 0
        Nerr = 0
    xcor = Ncor*scale

    #need to divide by mumu acceptance and multiply by 2mu2e acceptance
    acc2mu2e = ptAcc2mu2e.GetBinContent(j)
    print("pt bin %d; acc: %f +/- %f; acc2mu2e: %f"%(j, acc, accErr, acc2mu2e)) 
    if acc > 0:
        totEta += (nEta / acc)*acc2mu2e
        #errEta += ((unct / acc)*acc2mu2e)**2

    #print fitted yield, acceptance, and cross section for both 2022 BParking AND Run2 Scouting
    oldxs = hXsecRun2.GetBinContent(nb) #+6)
    oldxsUnc = hXsecRun2Unc.GetBinContent(nb) #+6)
    print("oldxs j=%d, nb=%d, lowpt = %f"%(j, nb, hXsecRun2.GetBinLowEdge(nb))) 
    if oldxs == 0:
        oldacc = 0
    else:
        oldacc = oldxsUnc / oldxs
    oldyield = oldxsUnc * oldlumi * bratio
    if not test:
        fout.write("%d-%d\t%.2f\t%.4f\t%.2f\t\t\t%.2f\t%.4f\t%.2f\n"%(i, int(nextbin), nEta, acc, xcor, oldyield, oldacc, oldxs))

    #hXsecCor.SetBinContent(i, xcor)
    hXsecCor.SetBinContent(nb, xcor)
    xerC = Nerr*scale
    #hXsecCor.SetBinError(i, xerC)
    hXsecCor.SetBinError(nb, xerC)
    nb += 1

#errEta = errEta**0.5
if not test:
    fout.close()
c1 = ROOT.TCanvas()
c1.cd()
print("hXsec:") 
hXsec.Print()
hXsecCor.GetXaxis().SetTitle("#eta meson p_{T} (GeV)")
hXsecCor.GetYaxis().SetTitle("d#sigma/dp_{T} [fb GeV^{-1}") 
hXsecCor.SetTitle("UNCORRECTED AND CORRECTED Xsections vs. p_{T} in 2022 BParking mu-mu")

hXsec.SetLineColor(ROOT.kMagenta)
hXsec.SetLineWidth(2)
c1.SetLogx()
c1.SetLogy()
ROOT.gStyle.SetOptStat(0)
hXsecCor.GetYaxis().SetRangeUser(10000, 100000000000)
if inc_corr:
    #first draw the corrected one (it's higher)
    hXsecCor.SetLineColor(ROOT.kGreen)
    hXsecCor.SetLineWidth(2)
    hXsecCor.Draw()

    #now do a fit to the corrected xsec
    fX = ROOT.TF1("fX", "[0] / x**[1]", 6, 100)
    #fX = ROOT.TF1("fX", "[0] * exp( - [1] * x ) + [2]", 6, 100)
    fX.SetParameter(0, 5.0*10**10)
    fX.SetParameter(1, 0.17)
    x_res = hXsecCor.Fit(fX, "RS", "", 10, 100)
    xparams = x_res.Parameters()
    xerrors = x_res.GetErrors()
    print("****fit result: params %s, errors %s****"%(str(xparams), str(xerrors))) 
    pav = ROOT.TPaveText(25, 500000000, 95, 80000000000)
    pav.AddText("Best fit: y = %.4f / x^{%.4f}"%(xparams[0], xparams[1])) 
    pav.Draw("same") 

    print("hXsec:") 
    hXsec.Print()
    hXsec.Draw("same")
    #is this necessary??
    fX.Draw("same")
else:
    hXsec.GetXaxis().SetTitle("#eta meson p_{T} (GeV)")
    hXsec.GetYaxis().SetTitle("d#sigma/dp_{T} [fb GeV^{-1}") 
    hXsec.SetTitle("UNCORRECTED Xsections vs. p_{T} in 2022 BParking mu-mu")
    hXsec.Draw()

if inc_run2:
    hXsecRun2.SetLineColor(ROOT.kBlue)
    hXsecRun2.SetLineWidth(2)
    hXsecRun2.Draw("same") 
    hXsecRun2Unc.SetLineColor(ROOT.kCyan)
    hXsecRun2Unc.SetMarkerColor(ROOT.kCyan)
    hXsecRun2Unc.SetLineWidth(2)
    hXsecRun2Unc.Draw("same") 

if inc_corr:
    leg = ROOT.TLegend()
    leg.AddEntry(hXsecCor, "Corrected by CMS efficiency")
    leg.AddEntry(hXsec, "Uncorrected")
    if inc_run2:
        leg.AddEntry(hXsecRun2, "Run2 corrected")
        leg.AddEntry(hXsecRun2Unc, "Run2 Uncorrected")
    leg.Draw("same")

c1.Modified()
c1.Update()
input("h")

if ( not test ) and inc_corr :
    outfile = ROOT.TFile.Open("xsec2022.root", "recreate")
    hXsec.Write()
    hXsecCor.Write()
    fX.Write()
    outfile.Close()

print("**TOTAL predicted eta->2mu2e decays seen for 2022 BParking: %f**"%(totEta/3)) #, errEta)) 
