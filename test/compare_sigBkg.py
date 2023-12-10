import ROOT

fsig = ROOT.TFile.Open("bparking_sigMCtest387.root")
fbkg = ROOT.TFile.Open("bparking_bkgMCtest389.root")

distname = "mmelel"
hsigAll = fsig.Get("hpTGenAll")
hsigAll.SetName("hsigAll")
hsigAcc = fsig.Get("hpTGenAcc%s"%distname) 
hsigAcc.SetName("hsigAcc") 

hbkgAll = fbkg.Get("hpTGenAll")
hbkgAcc = fbkg.Get("hpTGenAcc%s"%distname) 
