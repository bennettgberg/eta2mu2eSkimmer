import ROOT

#plot invariant mass SQUARED or just invariant mass?
square = False

m2min = 0.0
m2max = .6
m2bins = 120

if square:
    m2max = .4
    m2bins = 320 #80

#regions: peak, LSide (lower sideband), RSide (upper sideband) in invariant mass spectrum
region = "peak" #"RSide"
if region == "peak": region = ""

#f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest36_ALL.root")
#f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest366_ALL.root")
#f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest369_ALL.root")
#f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3610_ALL.root")
f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3611_ALL.root")
#f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3613_ALL.root")

hFull = f.Get("hMmmelel")

sq = ""
if square: sq = "2"
#make new hist to draw hFull only from .52 to .58 GeV
hmmelel2 = ROOT.TH1F("hM"+sq, "hM"+sq, m2bins, m2min, m2max) 
for i in range(hFull.GetNbinsX()+1):
    center = hFull.GetBinCenter(i)
    #if (region == "" and center >= .52 and center <= .58) or (region == "LSide" and center < .52 and i != 0) or (region == "RSide" and center > .58):
    if (region == "" and center >= .52 and center <= .58) or (region == "LSide" and center < .45 and i != 0) or (region == "RSide" and center > .65 and center < .75):
        if square: hmmelel2.Fill(center*center, hFull.GetBinContent(i)) 
        else: hmmelel2.Fill(center, hFull.GetBinContent(i)) 

hNoEl = f.Get("hMNoEl"+region)
hNoEl2 = ROOT.TH1F("hMNoEl2","m^{%s} with no electrons"%sq, m2bins, m2min, m2max)
for i in range(hNoEl.GetNbinsX()+1):
    center = hNoEl.GetBinCenter(i)
    m2 = center*center
    if not square: m2 = center
    hNoEl2.Fill(m2, hNoEl.GetBinContent(i)) 
hNoEl2.SetLineWidth(2)
hNoEl2.SetLineColor(ROOT.kRed)

hNoMu = f.Get("hMNoMu"+region)
hNoMu2 = ROOT.TH1F("hMNoMu2","m^{%s} with no muons"%sq, m2bins, m2min, m2max)
for i in range(hNoMu.GetNbinsX()+1):
    center = hNoMu.GetBinCenter(i)
    m2 = center*center
    if not square: m2 = center
    hNoMu2.Fill(m2, hNoMu.GetBinContent(i)) 
hNoMu2.SetLineWidth(2)
hNoMu2.SetLineColor(ROOT.kGreen)

hmmelel2.SetLineColor(ROOT.kBlack)
hmmelel2.SetLineWidth(2)

c = ROOT.TCanvas()

hmmelel2.Draw("hist")
hNoEl2.Draw("hist same")
hNoMu2.Draw("hist same")

hmmelel2.GetXaxis().SetTitle("%sInvariant mass (GeV^{%s})"%("Squared " if square else "", sq))
hmmelel2.GetYaxis().SetTitle("Events / %f GeV^{%s}"%(.005 if not square else .00125, sq))
if square:
    hmmelel2.GetYaxis().SetRangeUser(0, 1300)

leg = ROOT.TLegend()
leg.AddEntry(hmmelel2, "Full #mu#muee")
leg.AddEntry(hNoEl2, "#mu#mu only")
leg.AddEntry(hNoMu2, "ee only")
leg.Draw("same")

c.Modified()
c.Update()

input("holup") 
