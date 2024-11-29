import ROOT
import os

testnum = 4771 
lets = ["C", "Dv1", "Dv2", "E", "F", "G"]

for j in lets :
    for i in range(8):
        #TODO: delet this line!!!!!
        if not (j == "F" and i == 6): continue
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

