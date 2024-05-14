import ROOT
import array

etamass = .547862

#plot either the eta->2mu2e yield or the fitted eta meson mass as a function of the momentum scale factor s
#get the values from the saved txt files

#plot mass peak instead of relative yield?
plot_m = False #True

#require muon ID, or nah?
req_muID = True

#including trigger eff corrections, or nah?
do_trigCor = True

#how many points of s were tried?
npoints = 12 + 1
#what is the maximum value of s tried?
smax = 0.03
step = 2*smax / (npoints-1)

#xvals
sp1s = array.array('f', [0.0 for i in range(npoints)])
yvals = array.array('f', [0.0 for i in range(npoints)])
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
    fname = "fit_results/sigMCParams_DoubleGauss_newWt"
    if np != 6:
        fname += "_mod%d"%(n) 
    fname += "%s%s.txt"%("_reqMuID" if req_muID else "", "_trigCor" if do_trigCor else "")
    f = open(fname, "r")
    for line in f:
        words = line.split()
        if plot_m and "mg" in words[0]:
            yvals[np] = float(words[1]) / etamass
            break
        elif (not plot_m) and "nsig" in words[0]:
            yvals[np] = float(words[1])
            break
nomval = yvals[6]
if not plot_m:
    for np in range(npoints):
        yvals[np] /= nomval
tg = ROOT.TGraph(npoints, sp1s, yvals)
tg.GetXaxis().SetTitle("(1+s)")
if plot_m:
    ytitle = "M_{2#mu2e}/M_{#eta}"
else:
    ytitle = "Relative #eta #rightarrow 2#mu2e Yield"
tg.GetYaxis().SetTitle(ytitle)
tg.SetTitle("")
tg.SetMarkerStyle(8)
tg.SetMarkerColor(ROOT.kBlue)
tg.Draw("APE")

if plot_m and False:
    #fit to a straight line
    tf = ROOT.TF1("linfit", "[0]+[1]*x", 1.0-smax, 1.0+smax)
    #tf = ROOT.TF1("linfit", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 1.0-smax, 1.0+smax)
    tg.Fit(tf)
    tf.Draw("same")

input("h...") 
