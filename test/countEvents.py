import ROOT

lets = ['C', 'D', 'E', 'F', 'G']
ntot = {}

for let in lets:
    ntot[let] = 0

for let in lets:
    allf = open("allfiles_%s.txt"%(let), "r")
    for fn in allf:
        fname = "root://cmseos.fnal.gov/" + fn.strip()
        print("Opening " + fname)
        f = ROOT.TFile.Open(fname)
        t = f.Get("ntuples/recoT")
        n = t.GetEntries()
        ntot[let] += n
    print("Total for %s: %d"%(let, ntot[let])) 

ttot = 0
for let in lets:
    print("Tot for %s: %d"%(let, ntot[let]))
    ttot += ntot[let]
print("Grand total: %d"%(ttot))
