
import os

print("Nonexistent files:")
for let in ['C', 'D', 'E', 'F', 'G']:
    for setn in range(8):
        for num in range(32):
            if os.path.exists("filelists/flist_%s%d_%d.txt"%(let, setn, num)):
                os.system("eosls /store/user/bgreenbe/BParking2022/bparking_datatest3860_%s%d_%d.root >>temp.txt"%(let, setn, num)) 