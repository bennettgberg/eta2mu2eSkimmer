import ROOT

#todo: finish the last few lines of this program.

let = "C"
setn = 0
maxn = 1278
allf = ["root://cmseos.fnal.gov//store/user/lpcdisptau/eta2mu2e/BParking_2022/ParkingDoubleMuonLowMass%d%s/all_ParkingDoubleMuonLowMass%d%s_%s.root"%(setn, let, setn, let, str(i+1).zfill(3)) for i in range(maxn)]
#allf = ["root://cmseos.fnal.gov//store/user/lpcdisptau/eta2mu2e/BParking_2022_noL1/ParkingDoubleMuonLowMass%d%s/all_ParkingDoubleMuonLowMass%d%s_%s.root"%(setn, let, setn, let, str(i+1).zfill(3)) for i in range(0, maxn, 10)]

outf = open("%d%s_nEvents.txt"%(setn, let), "w")
totn = 0
#total events in the group of 10 files (to compare)
tot10 = 0
for i,f in enumerate(allf):
    if i % 10 == 0:
        outf.write("%d: %d\n"%(i+1, tot10))
        tot10 = 0
    print(f)
    ff = ROOT.TFile.Open(f)
    t = ff.Get("ntuples/recoT")
    n = t.GetEntries()
    if n == 0:
        print("Error: 0 events")
        exit()
    totn += n


print("%s %d totn: %d"%(let, setn, totn)) 
outf.close()
