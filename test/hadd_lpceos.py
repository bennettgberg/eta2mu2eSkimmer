import ROOT
import os

testnum = 4734
#cmd = "hadd -f root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/sigMC_test44_0.root "
#cmd = "hadd -f root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/dataDoubleElectron_test44_0C.root "
#cmd = "hadd -f root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/dataDoubleElectron_test44_0C_1.root "
#cmd = "hadd -f root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/dataDoubleElectron_test44_0C_2.root "
#for i in range(1, 391):
#for i in range(1, 1279):
#for i in range(301, 901):
#for i in range(901, 1279):
    #cmd += "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/CRAB_UserFiles/test44MC_EtaTo2Mu2E_0/240229_003322/0000/output_%d.root "%i
    #cmd += "root://cmseos.fnal.gov//store/user/bgreenbe/BParking_2022/ParkingDoubleMuonLowMass0C/all_ParkingDoubleMuonLowMass0C_%s.root "%(str(i).zfill(3)) 
lets = ["C", "Dv1", "Dv2", "E", "F", "G"]
for i in range(8):
    for j in lets :
        #TODO: delet this line!!!!!
        #if i < 5 and j != "G": continue
        #if i == 0 or (i == 1 and j in ['C', 'Dv1', 'Dv2']): continue
        #if not ( i == 5 and j == 'G'): continue
        cmd = "hadd -f root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest%d_%d%s.root "%(testnum, i, j)
        for k in range(32):
            if os.path.exists("filelists/flist_%d%s_%d.txt"%(i, j, k)):
                cmd += "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_test%d_%s%d_%d.root "%(testnum, j, i, k)

        print("cmd: ")
        print(cmd)
        os.system(cmd)

#now hadd errything
cmd = "hadd -f root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest%d_ALL.root "%testnum
for i in range(8):
    for j in lets:
        cmd += "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest%d_%d%s.root "%(testnum, i, j)
print("cmd: ")
print(cmd)
os.system(cmd)

