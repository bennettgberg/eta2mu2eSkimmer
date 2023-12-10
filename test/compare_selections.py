import ROOT

#Loosest cuts
#fNoCuts = ROOT.TFile.Open("bparking_bkgMCtest360.root")
#fNoCuts = ROOT.TFile.Open("bparking_bkgMCtest3811.root")
fNoCuts = ROOT.TFile.Open("bparking_bkgMCtest3818.root")
hNoCuts = fNoCuts.Get("hMmmelel")
hNoCuts.SetName("hNoCuts")
hNoCuts.SetTitle("Resonant Background Distribution Comparison")

#Strictest cuts
#fIdCuts = ROOT.TFile.Open("bparking_bkgMCtest389.root")
fIdCuts = ROOT.TFile.Open("bparking_bkgMCtest3819.root")
hIdCuts = fIdCuts.Get("hMmmelel")
hIdCuts.SetName("hIdCuts")
hIdCuts.SetTitle("hIdCuts")

#In-between cuts
#fBasicCuts = ROOT.TFile.Open("bparking_bkgMCtest3812.root")
fBasicCuts = ROOT.TFile.Open("bparking_bkgMCtest3820.root")
hBasicCuts = fBasicCuts.Get("hMmmelel")
hBasicCuts.SetName("hBasicCuts")

hNoCuts.SetLineWidth(2)
hIdCuts.SetLineWidth(2)
hBasicCuts.SetLineWidth(2)

hNoCuts.SetLineColor(ROOT.kBlack)
hIdCuts.SetLineColor(ROOT.kGreen)
hBasicCuts.SetLineColor(ROOT.kRed)

rebin = 5
hNoCuts.Rebin(rebin)
hIdCuts.Rebin(rebin)
hBasicCuts.Rebin(rebin)

hNoCuts.GetXaxis().SetTitle("m_{#mu#muee} (GeV)")
hNoCuts.GetYaxis().SetTitle("Events / .005 GeV")

c = ROOT.TCanvas()
c.cd()
hNoCuts.Draw()
hBasicCuts.Draw("same")
hIdCuts.Draw("same")

leg = ROOT.TLegend()
#leg.AddEntry(hNoCuts, "No selections")
leg.AddEntry(hNoCuts, "vProb > .1 ; nMissingHits <= 3")
leg.AddEntry(hBasicCuts, "vProb > .5 ; nMissingHits <= 3")
leg.AddEntry(hIdCuts, "vProb > .5 ; nMissingHits == 0")
leg.Draw("same")

c.Modified()
c.Update()
input("h")
