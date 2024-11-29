
import os
import ROOT

testn = 480

print("Nonexistent files:")
#for let in ['C', 'Dv1', 'Dv2', 'E', 'F', 'G']:
for let in ['C']:
#for let in ['E', 'F', 'G']:
    for setn in range(1): #8):
        for num in range(32):
            #if os.path.exists("/afs/cern.ch/user/b/bgreenbe/private/filelists/flist_%s%d_%d.txt"%(let, setn, num)):
            if os.path.exists("/afs/cern.ch/user/b/bgreenbe/private/filelists/flist_%d%s_%d.txt"%(setn, let, num)):
                #os.system("eosls /store/user/bgreenbe/BParking2022/bparking_datatest%d_%s%d_%d.root >>temp.txt"%(testn, let, setn, num)) 
                #print("eos root://cmseos.fnal.gov/ ls /store/user/bgreenbe/BParking2022/ultraskimmed/bparking_test%d_%s%d_%d.root >>temp.txt"%(testn, let, setn, num)) 
                os.system("eos root://cmseos.fnal.gov/ ls /store/user/bgreenbe/BParking2022/ultraskimmed/bparking_test%d_%s%d_%d.root >>temp.txt"%(testn, let, setn, num)) 
                #fnew = ROOT.TFile.Open("root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_test%d_%s%d_%d.root"%(38130, let, setn, num)) 
                #hNnew = fnew.Get("hNevt")
                #nNew = hNnew.GetEntries()
                #if nNew == 0:
                #    print("%s %d %d : 0 events!!"%(let, setn, num))
