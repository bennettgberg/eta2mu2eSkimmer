import ROOT
import sys

#reduce bkgMC files to include only events with at least two good muons and two good (regular) electrons.

arg = int(sys.argv[1]) 

#inf = open("bkgMCList%d.txt"%arg, "r")
inf = open("bkgMCFullLists/bkgMCList%d.txt"%arg, "r")
#inf = open("bkgMCCentralList%d.txt"%arg, "r")

for line in inf:
    fullpath = line.strip()
    words = fullpath.split("/")
    fname = words[-1]

    print("Starting %s"%fname) 

    partpath = "/".join(words[:-1])

    #outname = "root://cmseos.fnal.gov/"+partpath + "/reduced_" + fname
    outname = "root://cmseos.fnal.gov/"+partpath + "/reducedlpt_" + fname
    print("Opening outfile %s"%outname) 
    fnew = ROOT.TFile.Open(outname, "recreate") 

    f = ROOT.TFile.Open("root://cmseos.fnal.gov/" + fullpath)
    if not f:
        print("lkjdsalfhakjdshflkjdsahflkjhdslkufhdskjhflkjfuckalkjdflkjdhlkjfdhlkjhalkjhdfkjhd")
        exit()

    hpTGenAllEta = f.Get("ntuples/allGenPtEta") 
    hpTGenAllEtaClone = hpTGenAllEta.Clone()
    hpTGenAllEtaClone.SetDirectory(fnew)

    recoT = f.Get("ntuples/recoT")
    genT = f.Get("ntuples/genT")

    newRecoT = recoT.CloneTree(0)
    newRecoT.SetDirectory(fnew)
    newGenT = genT.CloneTree(0)
    newGenT.SetDirectory(fnew)

    nTot = recoT.GetEntries()
    for i,e in enumerate(recoT):
        if i%10000 == 0: print("Processing event %d / %d"%(i, nTot)) 
        #if ord(e.nGoodMuon) > 1 and ord(e.nGoodElectron) > 1:
        if ord(e.nGoodMuon) > 1 and ord(e.nGoodLowPtElectron) > 1:
            newRecoT.Fill()
            genT.GetEntry(i)
            newGenT.Fill()

    fnew.cd()
    newRecoT.Write()
    newGenT.Write()
    hpTGenAllEtaClone.Write()
    fnew.Close()

print("Done reducing.") 
