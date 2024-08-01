import os,sys
import math

#run letters of 2022 to go over
lets = ['C', 'Dv1', 'Dv2', 'E', 'F', 'G']
#lets = ['F']
#which sets of LowMass BParking data to go over
nums = [i for i in range(8)]
##which testnumber ?
#testnum = 36

njobs = 32

for let in lets:
    #for num in nums:
    for num in [1]:
        nfiles = 0
        flist = []
        #dpath = "BParking_2022/ParkingDoubleMuonLowMass%d%s"%(num, let)
        dpath = "BParking_2022_old/ParkingDoubleMuonLowMass%d%s"%(num, let)
        if not os.path.exists(dpath):
            continue
        while True:
            fpath = "%s/ParkingDoubleMuonLowMass%d%s_%s.jdl"%(dpath, num, let, str(nfiles+1).zfill(3))
            if not os.path.exists(fpath):
                break

            rpath = "/store/group/lpcdisptau/eta2mu2e/BParking_2022/ParkingDoubleMuonLowMass%d%s/all_ParkingDoubleMuonLowMass%d%s_%s.root"%(num, let, num, let, str(nfiles+1).zfill(3))
            flist.append(rpath)
            nfiles += 1

        #min number of files per job
        nf = int(nfiles/njobs)
        #number of jobs with one extra file
        nextra = nfiles % njobs
        kfiles = 0
        for j in range(njobs):
            outn = "filelists_new/flist_%s%d_%d.txt"%(let, num, j)
            outf = open(outn, 'w')
            jfiles = nf
            if j < nextra:
                jfiles += 1 
            for k in range(jfiles):
                outf.write(flist[kfiles]+'\n')
                kfiles += 1
            outf.close()
        if kfiles != nfiles:
            print("aaaaa kfiles=%d at end but nfiles=%d aaaaaaaaa"%(kfiles, nfiles)) 
