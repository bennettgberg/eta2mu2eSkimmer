import ROOT

mumu = False #True

#which version: central pseudorapidity only (1), outer only (2), or combined (0)?
region = 0

##factors by which to rebin in x and y
#they have different nbins actually so don't use this
#rebinx = 4
#rebiny = 2

if mumu:
    fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_mumuMCtest3886.root")
else:
    fsig = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_sigMCtest3884.root")
hsigAll = fsig.Get("hsubVdRGenAll") 
hsigAll.SetName("hsigAll")
hsigTrig = fsig.Get("hsubVdRGenTrig")
hsigTrig.SetName("hsigTrig")

#denominator, numerator come from different files for data
#fdataall = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/trigEff/trigEff_ALL2022.root") 
fdataall = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/trigEff_2GeVmaxFR/trigEffD_All2022.root") 
#fdataone = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/trigEff/trigEff0_ALL2022.root") 
fdataone = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/trigEff_2GeVmaxFR/trigEff0D_0_2022.root") 
#hdataAll = fdataone.Get("hsubVdRAll") 
if region == 0 or region == 1:
    hdataAll = fdataone.Get("hsubVdRAllCenter") 
    hdataTrig = fdataall.Get("hsubVdRTrigCenter") 
    if region == 1:
        hdataAll.SetTitle("|#eta| < 1.4")
elif region == 2:
    hdataAll = fdataone.Get("hsubVdRAllOuter") 
    hdataTrig = fdataall.Get("hsubVdRTrigOuter") 
    hdataAll.SetTitle("|#eta| > 1.4")
if region == 0:
    hdataAllOut = fdataone.Get("hsubVdRAllOuter") 
    hdataAll.Add(hdataAllOut)
    hdataTrigOut = fdataall.Get("hsubVdRTrigOuter") 
    hdataTrig.Add(hdataTrigOut)
    hdataAll.SetTitle("Inclusive in pseudorapidity")

#if rebinx > 1:
if 0 == 0:
    hsigAll.RebinX(10)
    hsigTrig.RebinX(10)
    hsigAll.RebinY(20)
    hsigTrig.RebinY(20)
    hdataAll.RebinX(2)
    hdataTrig.RebinX(2)
    hdataAll.RebinY(2)
    hdataTrig.RebinY(2)

hdataAll.Sumw2()
hsigAll.Sumw2()
hsigTrig.Sumw2()
hdataTrig.Sumw2()

#Divide to get the signal and data efficiencies
hsigTrig.Divide(hsigAll)
hsigTrig.Draw("colz text E")

input("Display data eff?")

hdataTrig.Divide(hdataAll)
hdataTrig.Draw("colz text E")

input("Display ratio (correction)?")

#Now divide the two efficiencies to get the corrections
hdataTrig.Divide(hsigTrig)

#Draw the 2d histogram of corrections
c = ROOT.TCanvas()
hdataTrig.Draw("colz text E")
hdataTrig.SetStats(ROOT.kFALSE)
c.SetLogz()
input("Write to output file? Ctrl-C to cancel.")

#Write the result to an output file
outname = "trigger_corrections_%s2022.root"%("mumuMC" if mumu else "sigMC")
fout = ROOT.TFile.Open(outname, "recreate")
hdataTrig.SetName("%sMCcorrection"%("mumu" if mumu else "sig"))
hdataTrig.Write()
print("Wrote to output file %s."%(outname)) 
fout.Close()
