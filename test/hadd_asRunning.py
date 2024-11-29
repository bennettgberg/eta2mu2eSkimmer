import ROOT
import os

testnum = 4787
backup = False

lets = ['C', 'Dv1', 'Dv2', 'E', 'F', 'G']
hadded = {} 
running = {}
for let in lets:
    hadded[let] = [False for i in range(8)] 
    running[let] = [False for i in range(8)] 

os.system("condor_q -nobatch >current_condor.txt")

curf = open("current_condor.txt", "r")

for line in curf:
    words = line.strip().split()
    if len(words) < 12: continue
    if words[1] != "bgreenbe": continue
    arg = int(words[-1])
    num = int(words[-2])
    let = words[-3]
    script = words[-5]
    if "sub_jobs" not in script:
        continue
    if backup and "backup" not in script:
        continue
    elif not backup and "backup" in script:
        continue
    running[let][num] = True

#now find out what has already been hadded
os.system("eos root://cmseos.fnal.gov/ ls /store/user/bgreenbe/BParking2022/ultraskimmed/ | grep datatest%d >hadded.txt"%(testnum)) 
curh = open("hadded.txt", "r")
for line in curh:
    if "bparking_datatest" not in line:
        continue
    #line should be: bparking_datatest[num]_[let][num].root
    words = line.strip().split("_")
    data = words[-1].split(".")[0]
    num = int(data[0]) 
    let = data[1:]
    hadded[let][num] = True

#now do the hadding for whatever is needed
for let in lets:
    for num in range(8):
        if not hadded[let][num] and not running[let][num]:
            #if not yet hadded and not running anymore, do the hadding!
            cmd = "hadd -f root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_datatest%d_%d%s.root "%(testnum, num, let)
            nfiles = 32
            if let == 'F':
                nfiles = 64
            for k in range(32):
                if os.path.exists("filelists/flist_%d%s_%d.txt"%(num, let, k)):
                    cmd += "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/ultraskimmed/bparking_test%d_%s%d_%d.root "%(testnum, let, num, k)

            print("cmd: ")
            print(cmd)
            os.system(cmd)
