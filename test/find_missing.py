
import os

print("Nonexistent files:")
#for let in ['C', 'D', 'E', 'F', 'G']:
for let in ['C', 'Dv1', 'Dv2', 'E', 'F', 'G']:
    for num in range(8):
        fnum = 1
        #how many root files per job were used
        fperj = 1 
        #if let == 'C' and num == 0: fperj = 5
        while True:
            #if let == 'E' or let == 'F':
            if 2+2 == 5:
                #fname = "/store/user/bgreenbe/BParking_2022/ParkingDoubleMuonLowMass%d%s/all_ParkingDoubleMuonLowMass%d%s_%s.root"%(num, let, num, let, str(fnum).zfill(3))
                fname = "/store/user/bgreenbe/BParking_2022/ParkingDoubleElectronLowMass%d%s/all_ParkingDoubleElectronLowMass%d%s_%s.root"%(num, let, num, let, str(fnum).zfill(3))
            else:
                fname = "/store/user/lpcdisptau/eta2mu2e/BParking_2022/ParkingDoubleMuonLowMass%d%s/all_ParkingDoubleMuonLowMass%d%s_%s.root"%(num, let, num, let, str(fnum).zfill(3))
            #first make sure the jdl file exists
            if not os.path.exists("BParking_2022/ParkingDoubleMuonLowMass%d%s/ParkingDoubleMuonLowMass%d%s_%s.jdl"%(num, let, num, let, str(fnum).zfill(3))):
            #if not os.path.exists("BParking_2022/ParkingDoubleElectronLowMass%d%s/ParkingDoubleElectronLowMass%d%s_%s.jdl"%(num, let, num, let, str(fnum).zfill(3))):
                break
            #now if so, see if the root file exists
            try:
                os.system("xrdfs root://cmseos.fnal.gov ls %s >>allfiles.txt"%fname)
            except:
                print(fname)
                print("kjlk;kjk;lkjlkjhkj;lkj;lkjlkjahjdflkjhdjfhaj;hfkjdlhfaj;") 
            fnum += fperj
