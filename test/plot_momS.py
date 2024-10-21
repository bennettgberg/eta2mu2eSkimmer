import ROOT
import array

etamass = .547862

#plot either the eta->2mu2e yield or the fitted eta meson mass as a function of the momentum scale factor s
#get the values from the saved txt files

#plot mass peak instead of relative yield?
plot_m = True

#do for mumu instead of sig?
mumu = True
#plot both mumu AND sig?
both = True

#require muon ID, or nah?
req_muID = True

#including trigger eff corrections, or nah?
do_trigCor = True

#do correction for muons instead of electrons?
do_mu = True

testnummu = "38171p4737"
testnumel = "38172p4738"
if do_mu:
    #testnum = "38137p478"
    testnum = testnummu 
else:
    #testnum = "38145p4715"
    testnum = testnumel

#how many points of s were tried?
npoints = 12 + 1
#what is the maximum value of s tried?
if do_mu:
    smax = 0.003
else:
    smax = 0.012 #0.03
step = 2*smax / (npoints-1)

if not plot_m:
    fname1 = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest%s.root"%("mumu" if mumu else "sig", testnum)
    print("Opening %s"%(fname1))
    f1 = ROOT.TFile.Open(fname1)
    if both:
        testnum = testnumel if do_mu else testnummu
        fname2 = "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_%sMCtest%s.root"%("sig" if mumu else "mumu", testnum)
        f2 = ROOT.TFile.Open(fname2)
#xvals
sp1s = array.array('f', [0.0 for i in range(npoints)])
yvals = array.array('f', [0.0 for i in range(npoints)])
yerrs  = array.array('f', [0.0 for i in range(npoints)])
if both:
    yvals2 = array.array('f', [0.0 for i in range(npoints)])
    yerrs2 = array.array('f', [0.0 for i in range(npoints)])
    
ymin = 9999
ymax = 0
for np in range(npoints):
    if np == 6:
        #6 is nominal
        n = -1
    elif np > 5:
        n = np-1
    else:
        n = np
    
    sp1s[np] = 1.0-smax + step*np
    #open the file
    #yvals[np] = 1.0+step*(np-6)
    if plot_m:
        fname = "fit_results/%sMCParams_DoubleGauss_newWt"%("mumu" if mumu else "sig")
        if both:
            fname2 = "fit_results/%sMCParams_DoubleGauss_newWt"%("mumu" if not mumu else "sig")
        if np != 6:
            fname += "_mod%s%d"%("mu" if do_mu else "", n) 
            if both:
                fname2 += "_mod%s%d"%("mu" if do_mu else "", n) 
        fname += "%s%s.txt"%("_reqMuID" if req_muID else "", "_trigCor" if do_trigCor else "")
        if both:
            fname2 += "%s%s.txt"%("_reqMuID" if req_muID else "", "_trigCor" if do_trigCor else "")
        f = open(fname, "r")
        if both:
            f2 = open(fname2, "r")
        for line in f:
            words = line.split()
            if plot_m and "mg" in words[0]:
                yvals[np] = float(words[1]) / etamass
                yerrs[np] = float(words[2]) / etamass
                if yvals[np] > ymax:
                    ymax = yvals[np]
                if yvals[np] < ymin:
                    ymin = yvals[np]
                break
            elif (not plot_m) and "nsig" in words[0]:
                yvals[np] = float(words[1])
                yerrs[np] = float(words[2])
                break
        if both:
            for line in f2:
                words = line.split()
                if plot_m and "mg" in words[0]:
                    yvals2[np] = float(words[1]) / etamass
                    yerrs2[np] = float(words[2]) / etamass
                    break
                elif (not plot_m) and "nsig" in words[0]:
                    yvals2[np] = float(words[1])
                    yerrs2[np] = float(words[2])
                    break
    else:
        #just get the total sum of weights from the root file
        if np == 6:
            if mumu:
                distname = "hMmumu"
            else:
                distname = "hMmmelel"
        else:
            distname = "hMMod%s"%("mu" if do_mu else "") + str(n)
        h = f1.Get(distname)
        yvals[np] = h.GetSumOfWeights()
        yerrs[np] = yvals[np]**0.5
        #print("yvals[%d] = %f"%(np, yvals[np])) 
        if both:
            #do the other way now!
            if np == 6:
                if mumu:
                    distname2 = "hMmmelel"
                else:
                    distname2 = "hMmumu"
            else:
                distname2 = "hMModmu" + str(n)
            h2 = f2.Get(distname2)
            yvals2[np] = h2.GetSumOfWeights()
            yerrs2[np] = h2.GetSumOfWeights()**0.5
            print("yvals2[%d] = %f"%(np, yvals2[np])) 
nomval = yvals[6]
print("0 nomval: %f"%nomval)
if both:
    nomval2 = yvals2[6]
if not plot_m:
    for np in range(npoints):
        yvals[np] /= nomval
        yerrs[np] /= nomval
        print("new yvals[%d] = %f"%(np, yvals[np]))
        if both:
            yvals2[np] /= nomval2
            yerrs2[np] /= nomval2
            print("new yvals2[%d] = %f"%(np, yvals2[np]))
print("1 nomval: %f"%nomval)
if plot_m:
#if True:
    tg = ROOT.TGraphErrors(npoints, sp1s, yvals, 0, yerrs)
    tg.SetLineColor(ROOT.kBlue)
else:
    tg = ROOT.TGraph(npoints, sp1s, yvals)
tg.GetXaxis().SetTitle("(1+s_{%s})"%("#mu" if do_mu else ""))
if plot_m:
    if both:
        ytitle = "M_{#eta}^{fitted}/M_{#eta}^{PDG}"
    else:
        ytitle = "M_{2#mu%s}/M_{#eta}"%("" if mumu else "2e")
else:
    if both:
        ytitle = "Relative #eta Yield"
    elif mumu:
        ytitle = "Relative #eta #rightarrow 2#mu Yield"
    else:
        ytitle = "Relative #eta #rightarrow 2#mu2e Yield"
tg.GetYaxis().SetTitle(ytitle)
tg.SetTitle("")
tg.SetName("tg1")
tg.SetMarkerStyle(8)
tg.SetMarkerColor(ROOT.kBlue)
tg.Draw("APE")
xmin = 1.0-smax - step
xmax = 1.0+smax + step
if plot_m and (both or not mumu):
    #draw lines for the data
    fdata = open("fit_results/twoMu2EDataFit_newWt_modmu3_muID_trigCor_DoubleGauss_Threshold2mu2eBreitWignerBkg.txt", "r")
    for line in fdata:
        words = line.split()
        if "mg" in words[0]:
            mdata = float(words[1])
            edata = float(words[2])
            break
    lcent = ROOT.TLine(xmin, mdata/etamass, xmax, mdata/etamass)
    lcent.SetLineColor(ROOT.kBlack)
    lcent.SetLineWidth(ROOT.kBlack)
    lcent.Draw("C same")
    lup = ROOT.TLine(xmin, (mdata+edata)/etamass, xmax, (mdata+edata)/etamass)
    lup.SetLineColor(ROOT.kBlack)
    lup.SetLineStyle(9) #dashed
    lup.Draw("same")
    ldn = ROOT.TLine(xmin, (mdata-edata)/etamass, xmax, (mdata-edata)/etamass)
    ldn.SetLineColor(ROOT.kBlack)
    ldn.SetLineStyle(9) #dashed
    ldn.Draw("same")
if plot_m and not both:
    #draw the line for the chosen correction.
    if do_mu:
        cor = 1 - 0.0014 #0.0015
    else:
        cor = 1 + 0.001
    lcor = ROOT.TLine(cor, ymin-.002, cor, ymax+.002)
    lcor.SetLineColor(ROOT.kGreen)
    lcor.Draw("C same")
if plot_m and (both or mumu):
    #now draw the mumu data line
    lmumu = ROOT.TLine(xmin, 5.4739e-01/etamass, xmax, 5.4739e-01/etamass)
    lmumu.SetLineWidth(2)
    lmumu.SetLineColor(ROOT.kViolet)
    lmumu.Draw("same")
if both:
    if plot_m:
        tg2 = ROOT.TGraphErrors(npoints, sp1s, yvals2, 0, yerrs2) 
        tg2.SetLineColor(ROOT.kRed)
    else:
        tg2 = ROOT.TGraph(npoints, sp1s, yvals2) 
    tg2.SetName("tg2")
    tg2.SetMarkerStyle(8)
    tg2.SetMarkerColor(ROOT.kRed)
    tg2.Draw("PE same")
    leg = ROOT.TLegend()
    leg.AddEntry(tg, "#eta #rightarrow #mu#mu%s MC"%("" if mumu else "ee"))
    leg.AddEntry(tg2, "#eta #rightarrow #mu#mu%s MC"%("ee" if mumu else ""))
    if plot_m:
        leg.AddEntry(lcent, "#mu#muee Data", "l") 
        leg.AddEntry(lmumu, "#mu#mu Data", "l") 
        tg.GetYaxis().SetRangeUser(.995, 1.004)
    leg.Draw("same")

#if plot_m and do_mu and not both:
#    #fit to a straight line
#    tf = ROOT.TF1("linfit", "[0]+[1]*x", 1.0-smax, 1.0+smax)
#    #tf = ROOT.TF1("linfit", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 1.0-smax, 1.0+smax)
#    tg.Fit(tf)
#    tf.Draw("same")
#    tline = ROOT.TLine(1.0, 1.0-smax, 1.0, 1.0+smax)
#    tline.SetLineColor(ROOT.kGreen)
#    tline.Draw("same")

input("h...") 
