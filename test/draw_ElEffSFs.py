import ROOT
import numpy as np
import matplotlib.pyplot as plt

do_id = False
do_unc = False #True

#f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38174p4740.root")
f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest38222p4783.root")

efftype = "ID" if do_id else "Rec"
h = f.Get("hEl%ssfmmelel"%(efftype))
hUp = f.Get("hEl%ssfUpmmelel"%(efftype))
hDn = f.Get("hEl%ssfDnmmelel"%(efftype))

h.Scale(0.5)
hUp.Scale(0.5)
hDn.Scale(0.5)

nbins = h.GetNbinsX()*100
bin_edges = np.linspace(-2.5, 2.5, nbins + 1)
hvals = np.zeros(nbins)
hUpvals = np.zeros(nbins)
hDnvals = np.zeros(nbins)

for i in range(nbins):
    hvals[i] = h.GetBinContent(int(i/100)) 
    hUpvals[i] = hUp.GetBinContent(int(i/100))
    hDnvals[i] = hDn.GetBinContent(int(i/100)) 
    #some might accidentally only have a value in one file instead of both, oops
    if hvals[i] < .59:
        hvals[i] *= 2
        hUpvals[i] *= 2
        hDnvals[i] *= 2
    if hvals[i] == 0 and i > 0:
        hvals[i] = hvals[i-1]
        hUpvals[i] = hUpvals[i-1]
        hDnvals[i] = hDnvals[i-1]

fig, ax = plt.subplots(figsize=(10, 6))
if do_unc:
    ax.fill_between(bin_edges[:-1], hDnvals, hUpvals, color='blue', alpha=0.4, label='Uncertainties')
ax.step(bin_edges[:-1], hvals, where='mid', color='black', linewidth=1.5, linestyle='-')
if do_id:
    title = 'Electron ID Scale factors'
else:
    title = 'Electron Reco Scale factors'
if do_unc:
    title += ' with Uncertainties'
ax.set_title(title)
ax.set_xlabel('Pseudorapidity', fontsize=16)
ax.set_ylabel('Scale factor', fontsize=16)
#if not (do_unc and do_id):
#    ax.set_ylim(0.9, 1.1)
#else:
ax.set_ylim(0.8, 1.1)
plt.show()

