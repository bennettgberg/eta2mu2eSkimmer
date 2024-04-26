import ROOT

#no elID
fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3850.root")
##WP90 req'd on both
#fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3837.root")

#distname = "hMNoMu" #"mmelel"
distname = "hMNoEl" #"mmelel"
#hsigAll = fsig.Get("hpTGenAll")
#hsigAll.SetName("hsigAll")
#hsigAcc = fsig.Get("hpTGenAcc%s"%distname)
#hsigAcc.SetName("hsigAcc")
#
#hbkgAll = fbkg.Get("hpTGenAll")
#hbkgAcc = fbkg.Get("hpTGenAcc%s"%distname) 

hsig = fsig.Get(distname)
hsig.SetName("hsig")

#no elID
fbkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3850.root")
##WP90 req'd on both
#fbkg = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_bkgMCtest3837.root")
hbkg = fbkg.Get(distname)
hbkg.SetName("hbkg")

hsig.SetLineColor(ROOT.kBlue)
hbkg.SetLineColor(ROOT.kRed)
hsig.SetLineWidth(2)
hbkg.SetLineWidth(2)
hsig.SetStats(ROOT.kFALSE)
#hsig.SetTitle("Dielectron invariant mass of events in the signal peak: 0.52 < M_{#mu#muee} < 0.58") 
hsig.SetTitle("Dimuon invariant mass of events in the signal peak: 0.52 < M_{#mu#muee} < 0.58") 

binwidth = (hsig.GetXaxis().GetXmax() - hsig.GetXaxis().GetXmin()) / hsig.GetNbinsX()
hsig.GetXaxis().SetTitle("M_{ee} (GeV)")
hsig.GetYaxis().SetTitle("Events / %.3f GeV"%binwidth) 
hsig.Draw("hist")
hbkg.Draw("hist same")

leg = ROOT.TLegend()
leg.AddEntry(hsig, "Signal #eta#rightarrow#mu#muee")
leg.AddEntry(hbkg, "Resonant background #eta#rightarrow#mu#mu#gamma")
leg.Draw("same")

input("kjkjk")
