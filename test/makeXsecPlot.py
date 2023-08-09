import ROOT
from array import array

#if only testing, use smaller input file and do NOT overwrite output files
test = False

#make plot of xsec as a function of pT for eta->2mu!
if test:
    #f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/bparking_datatest335_1C.root", "read")
    f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/bparking_datatest335_0.root", "read")
else:
    f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest338_ALL.root", "read")

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
newbins = [i*1.0 for i in range(6, 31)]
for i in range(32, 41, 2):
    newbins.append(1.0*i)
for i in range(45, 56, 5):
    newbins.append(1.0*i)
newbins.append(70.0)
newbins.append(100.0)
print("new pT bins: " + str(newbins)) 
#last entry in newbins is the upper limit of the last bin, not a new bin
newnptbins = len(newbins)-1
### new bins ^^

#construct the 1d hist of xsec vs pT -- UNCORRECTED
#hXsec = ROOT.TH1F("hXsecUncor", "hXsecUncor", nptbins, ptmin, ptmax)
hXsec = ROOT.TH1F("hXsecUncor", "hXsecUncor", newnptbins, array('d', newbins))
hXsec.Sumw2()
#Xsec corrected for acceptance
#hXsecCor = ROOT.TH1F("hXsecCor", "hXsecCor", nptbins, ptmin, ptmax)
hXsecCor = ROOT.TH1F("hXsecCor", "hXsecCor", newnptbins, array('d', newbins))
hXsecCor.Sumw2()

#get efficiency histograms for correction purposes
f2mu = ROOT.TFile.Open("bparking_mumuMCtest338.root")
ptGen = f2mu.Get("hpTGenAll")
ptGen.SetName("hpTGenAllmumu")
ptGen.Rebin(5)
ptAcc = f2mu.Get("hpTGenAccmumu")
#bins should be 1 GeV wide!!!
ptAcc.Rebin(5)
ptAcc.Sumw2()
ptAcc.Divide(ptGen)

#get efficiency histograms for correction purposes
f2mu2e = ROOT.TFile.Open("bparking_sigMCtest335.root")
ptGen2mu2e = f2mu2e.Get("hpTGenAll")
ptGen2mu2e.Rebin(5)
#ptAcc2mu2e = f2mu2e.Get("hpTGenAccmmelel")
ptAcc2mu2e = f2mu2e.Get("hpTGenAccmmlplp")
#bins should be 1 GeV wide!!!
ptAcc2mu2e.Rebin(5)
ptAcc2mu2e.Sumw2()
ptAcc2mu2e.Divide(ptGen2mu2e)

#counter for which newbin we're on
nb = 1
nextbin = newbins[nb]

if not test:
    fout = open("xsec_compare.txt", "w")
    fout.write("pT bin\tBParking Fitted Yield\tBParking Acc\tBParking XS\t\t\tRun2 Scouting Yield\tRun2 Scouting Acc\tRun2 Scouting XS\n")

#total number of eta mesons expected
totEta = 0.0
#sum of squares of uncertainty on total number of eta mesons expected (will take sqrt at end)
errEta = 0.0

#need efficiency to correct!!
for j in range(nptbins-1):
    i = j+1
    #first few bins have acceptance of 0, so skip them altogether
    #if i < 7: # or i > 16: 
    #    #hXsec.SetBinContent(i, 0)
    #    #hXsec.SetBinError(i, 0)
    #    continue
    #skip any old pT bins that aren't the start of the new bin (these were used already).
    if i != nextbin:
        print("i=%d; %f < pT < %f; continuing."%(i, h2d.GetXaxis().GetBinLowEdge(i), h2d.GetXaxis().GetBinLowEdge(i)+h2d.GetXaxis().GetBinWidth(i)))
        continue
    #do fit to invariant mass spectrum within this pT bin!!
    nextbin = newbins[nb+1]
    #projection = h2d.ProjectionY("hx"+str(i),i+1, i+1)
    print("nb=%d: old pt bins %d thru %d"%(nb, i, int(nextbin)-1)) 
    projection = h2d.ProjectionY("hx"+str(nb), i, int(nextbin)-1)
    #factor by which to rebin
    if nb == 1 or i > 31:
        rebin = 40
    else:
        rebin = 20
    if test:
        rebin *= 2
    projection.Rebin(rebin)
    #print("projection: " + str(i))
    projection.Print()
    #now do a fit to pol2 + gaussian
    #fit_func = ROOT.TF1("fit_func"+str(i), "[0] + [1]*x + [2]*x*x + [3]*exp(-0.5*((x-[4])/[5])**2)", 0.45, 0.65)
    fit_func = ROOT.TF1("fit_func"+str(nb), "[0] + [1]*x + [2]*x*x + [3]*exp(-0.5*((x-[4])/[5])**2)", 0.45, 0.65)
    if test:
        fit_func.SetParameters(14000, .1, -.1, 1000, 0.548, 0.01)
    elif nb == 1:
        fit_func.SetParameters(400, .1, -.1, 50, 0.548, 0.01)
        #fit_func.SetRange(.5, .6)
    elif i < 16:
        fit_func.SetParameters(20000, .1, -.1, 5000, 0.548, 0.01)
    elif i < 32:
        fit_func.SetParameters(3000, .1, -.1, 1000, 0.548, 0.01)
    elif i < 38:
        fit_func.SetParameters(2000, .1, -.1, 500, 0.548, 0.01)
    else:
        fit_func.SetParameters(2000, .1, -.1, 100, 0.548, 0.01)
    fit_func.SetParLimits(3, 0, 1000000)
    fit_func.SetParLimits(5, 0, .1)
    fit_result = projection.Fit(fit_func, "RS")
    #now get the fit results, integrate to get the number of eta mesons
    params = fit_result.Parameters()
    errors = fit_result.GetErrors()
    print("params: %s"%str(params)) 
    print("errors: %s"%str(errors)) 
    s2pi = (2*3.14159265359)**0.5
    nEta = params[3] * params[5] * s2pi / projection.GetBinWidth(1)
    ptbinwidth = nextbin - 1.0*i
    if test and i == 11:
        projection.Draw()
        print("nEta = %f"%nEta)
        input("wait...")
        input("wait...")
    unct = ( (errors[3] / params[3])**2 + (errors[5] / params[5] )**2 )**0.5 * nEta
    print("nEta = %f +/- %f"%(nEta, unct)) 
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
    print("old pt bins %d to %d; setting bin content for bin %d"%(i, nextbin, nb))
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

    #now do the corrected one
    acc = ptAcc.GetBinContent(i)
    accErr = ptAcc.GetBinError(i)

    if acc > 0.:
        Ncor = nEta / acc
        Nerr = Ncor * ( (unct / nEta)**2 + (accErr / acc)**2 )**0.5
    else:
        Ncor = 0
        Nerr = 0
    xcor = Ncor*scale

    #need to divide by mumu acceptance and multiply by 2mu2e acceptance
    acc2mu2e = ptAcc2mu2e.GetBinContent(i)
    print("pt bin %d; acc: %f; acc2mu2e: %f"%(i, acc, acc2mu2e)) 
    if acc > 0:
        totEta += (nEta / acc)*acc2mu2e
        #errEta += ((unct / acc)*acc2mu2e)**2

    #print fitted yield, acceptance, and cross section for both 2022 BParking AND Run2 Scouting
    oldxs = hXsecRun2.GetBinContent(nb+6)
    oldxsUnc = hXsecRun2Unc.GetBinContent(nb+6)
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
hXsec.Print()
hXsecCor.GetXaxis().SetTitle("#eta meson pT (GeV)")
hXsecCor.GetYaxis().SetTitle("d#sigma/dp_{T} [fb GeV^{-1}") 
hXsecCor.SetTitle("UNCORRECTED AND CORRECTED Xsections vs. pT in 2022 BParking mu-mu")

#first draw the corrected one (it's higher)
hXsecCor.SetLineColor(ROOT.kGreen)
hXsecCor.SetLineWidth(2)
hXsecCor.Draw()

#now do a fit to the corrected xsec
fX = ROOT.TF1("fX", "[0] / x**[1]", 6, 100)
x_res = hXsecCor.Fit(fX, "RS")
xparams = x_res.Parameters()
xerrors = x_res.GetErrors()
print("****fit result: params %s, errors %s****"%(str(xparams), str(xerrors))) 

hXsec.SetLineColor(ROOT.kMagenta)
hXsec.SetLineWidth(2)
hXsec.Draw("same")
c1.SetLogx()
c1.SetLogy()
hXsecCor.GetYaxis().SetRangeUser(10000, 100000000000)

hXsecRun2.SetLineColor(ROOT.kBlue)
hXsecRun2.SetLineWidth(2)
hXsecRun2.Draw("same") 
hXsecRun2Unc.SetLineColor(ROOT.kCyan)
hXsecRun2Unc.SetMarkerColor(ROOT.kCyan)
hXsecRun2Unc.SetLineWidth(2)
hXsecRun2Unc.Draw("same") 

leg = ROOT.TLegend()
leg.AddEntry(hXsecCor, "Corrected by CMS efficiency")
leg.AddEntry(hXsec, "Uncorrected")
leg.AddEntry(hXsecRun2, "Run2 corrected")
leg.AddEntry(hXsecRun2Unc, "Run2 Uncorrected")
leg.Draw("same")

input("h")

if not test:
    outfile = ROOT.TFile.Open("xsec2022.root", "recreate")
    hXsec.Write()
    hXsecCor.Write()
    fX.Write()
    outfile.Close()

print("**TOTAL predicted eta->2mu2e decays seen for 2022 BParking: %f**"%(totEta/3)) #, errEta)) 
