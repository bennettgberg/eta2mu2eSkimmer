
import os

testn = 38126

print("Nonexistent files:")
for let in ['C', 'D', 'E', 'F', 'G']:
    for setn in range(8):
        for num in range(32):
            if os.path.exists("filelists/flist_%s%d_%d.txt"%(let, setn, num)):
                os.system("eosls /store/user/bgreenbe/BParking2022/bparking_datatest%d_%s%d_%d.root >>temp.txt"%(testn, let, setn, num)) 
