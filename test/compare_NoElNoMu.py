import ROOT

#plot invariant mass SQUARED or just invariant mass?
square = False

m2min = 0.0
m2max = 1.0 #.6
m2bins = 200 #120

if square:
    m2max = .4
    m2bins = 320 #80

#regions: peak, LSide (lower sideband), RSide (upper sideband) in invariant mass spectrum
region = "peak" #"Side"
if region == "peak": region = ""

subtract_side = False

##WP90 req'd on both electrons
#f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3837_ALL.root")
###singleVert only!
##f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3840_ALL.root")
##WP90 req'd on one electron only
#f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3838_ALL.root")
#No elID required
f = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest3850_ALL.root")

hFull = f.Get("hMmmelel")

sq = ""
if square: sq = "2"

#make new hist to draw hFull only from .52 to .58 GeV
hmmelel2 = ROOT.TH1F("hM"+sq, "hM"+sq, m2bins, m2min, m2max) 
if subtract_side:
    hside = ROOT.TH1F("hMside", "hMside", m2bins, m2min, m2max) 

for i in range(hFull.GetNbinsX()+1):
    center = hFull.GetBinCenter(i)
    #if (region == "" and center >= .52 and center <= .58) or (region == "LSide" and center < .52 and i != 0) or (region == "RSide" and center > .58):
    if (region == "" and center >= .52 and center <= .58) or (region == "LSide" and center < .45 and i != 0) or (region == "RSide" and center > .65 and center < .75) or (region == "Side" and ((center > .47 and center < .51) or (center > .60 and center < .65))):
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

rebin = 1 #2
if rebin > 1:
    hNoEl2.Rebin(rebin)
    hNoMu2.Rebin(rebin)
    hmmelel2.Rebin(rebin)

if subtract_side:
    hNoElSide = f.Get("hMNoElSide")
    hNoMuSide = f.Get("hMNoMuSide")
    if rebin > 1:
        hNoElSide.Rebin(rebin)
        hNoMuSide.Rebin(rebin)
    hNoEl2.Add(hNoElSide, -1.0)
    hNoMu2.Add(hNoMuSide, -1.0)

hNoMu2.SetLineWidth(2)
hNoMu2.SetLineColor(ROOT.kGreen+1)

hmmelel2.SetLineColor(ROOT.kBlack)
hmmelel2.SetLineWidth(2)

c = ROOT.TCanvas()

hmmelel2.SetStats(ROOT.kFALSE)
hmmelel2.SetTitle("")
hmmelel2.Draw("hist")
hNoEl2.Draw("hist same")
hNoMu2.Draw("hist same")

hmmelel2.GetXaxis().SetTitle("%sInvariant mass (GeV^{%s})"%("Squared " if square else "", sq))
hmmelel2.GetYaxis().SetTitle("Events / %s GeV^{%s}"%(("0.005" if rebin == 1 else "0.01") if not square else ".00125", sq))
if square:
    hmmelel2.GetYaxis().SetRangeUser(0, 1300)

leg = ROOT.TLegend()
leg.AddEntry(hmmelel2, "Full #mu#muee")
leg.AddEntry(hNoEl2, "#mu#mu only"+(" - sideband" if subtract_side else ""))
leg.AddEntry(hNoMu2, "ee only"+(" - sideband" if subtract_side else ""))
leg.Draw("same")

c.Modified()
c.Update()

input("holup") 
