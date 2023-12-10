
def printEvent(entry, isMC=False, gentry=None) :

    print("** Run={0:d} LS={1:d} Event={2:d}".format(entry.run,entry.lumi_sec,entry.evt))
    #if entry.nMuon > 0 :
    if ord(entry.nGoodMuon) > 0 :
        print("Muons\n # Q     Pt     Eta     Phi     Id")
        #for j in range(entry.nMuon) :
        for j in range(ord(entry.nGoodMuon)) :
            muSign = '+'
            if ord(entry.Muon_charge[j]) != 1 : muSign = '-'
            print("{0:2d} {1:2s}   {2:5.1f}   {3:6.2f}   {4:6.2f}   {5:d}".format(
                j,muSign,entry.Muon_pt[j],entry.Muon_eta[j],entry.Muon_phi[j],
                ord(entry.Muon_id[j])
                ))
    if ord(entry.nGoodElectron) > 0 :
        print("Electrons\n #  Q    Pt        Eta         Phi      Id   convVeto  nMissingHits" )
        #for j in range(entry.nElectron) :
        for j in range(ord(entry.nGoodElectron)) :
            elSign = '+'
            if ord(entry.Electron_charge[j]) != 1 : elSign = '-'
            print("{0:2d}{1:2s}   {2:5.1f}   {3:6.2f}   {4:6.2f}   {5:d}   {6:d}   {7:d}".format(
                j,elSign,entry.Electron_pt[j],entry.Electron_eta[j],entry.Electron_phi[j],
                ord(entry.Electron_id[j]), ord(entry.Electron_convVeto[j]), ord(entry.Electron_nMissingHits[j])
                ))

    if ord(entry.nGoodPhoton) > 0 :
        print("Photons\n # Pt   Eta   Phi ")
        for j in range(ord(entry.nGoodPhoton)) :
            print("{0:2d} {1:5.1f}{2:6.2f}{3:6.2f}".format(
                j,entry.Photon_pt[j],entry.Photon_eta[j],entry.Photon_phi[j]))

    if len(entry.Vertex_mmelel_dR) > 0:
        print("Vertices mu+mu-e+e-\n #   vx    vy    vz    sigmaVxy    chi2     ndof     dR    mu+    mu-    e+    e-")
        for j in range(len(entry.Vertex_mmelel_dR)):
            print("{0:d}  {1:.3f}   {2:.3f}   {3:.3f}   {4:.3f}   {5:.3f}   {11:.1f}   {6:.3f}   {7:d}   {8:d}   {9:d}   {10:d}".format(
                j,entry.Vertex_mmelel_vx[j], entry.Vertex_mmelel_vy[j], entry.Vertex_mmelel_vz[j], entry.Vertex_mmelel_sigmaVxy[j],
                entry.Vertex_mmelel_chi2[j], entry.Vertex_mmelel_dR[j], 
                ord(entry.Vertex_mmelel_muP[j]), ord(entry.Vertex_mmelel_muN[j]), ord(entry.Vertex_mmelel_eleP[j]), ord(entry.Vertex_mmelel_eleN[j]),
                entry.Vertex_mmelel_ndof[j]
                ))

    if len(entry.Vertex_mumu_dR) > 0:
        print("Vertices mu+mu-\n #      vx       vy       vz     sigmaVxy    chi2    ndof    dR     mu+     mu-")
        for j in range(len(entry.Vertex_mumu_dR)):
            print("{0:d}  {1:.3f}    {2:.3f}     {3:.3f}     {4:.3f}      {5:.3f}   {9:.1f}    {6:.3f}      {7:d}      {8:d}".format(
                j,entry.Vertex_mumu_vx[j], entry.Vertex_mumu_vy[j], entry.Vertex_mumu_vz[j], entry.Vertex_mumu_sigmaVxy[j],
                entry.Vertex_mumu_chi2[j], entry.Vertex_mumu_dR[j], 
                ord(entry.Vertex_mumu_muP[j]), ord(entry.Vertex_mumu_muN[j]), entry.Vertex_mumu_ndof[j]
                ))

    if len(entry.Vertex_elel_dR) > 0:
        print("Vertices e+e-\n #      vx       vy       vz     sigmaVxy    chi2    ndof    dR     el+     el-")
        for j in range(len(entry.Vertex_elel_dR)):
            print("{0:d}  {1:.3f}    {2:.3f}     {3:.3f}     {4:.3f}      {5:.3f}   {9:.1f}    {6:.3f}      {7:d}      {8:d}".format(
                j,entry.Vertex_elel_vx[j], entry.Vertex_elel_vy[j], entry.Vertex_elel_vz[j], entry.Vertex_elel_sigmaVxy[j],
                entry.Vertex_elel_chi2[j], entry.Vertex_elel_dR[j], 
                ord(entry.Vertex_elel_eleP[j]), ord(entry.Vertex_elel_eleN[j]), entry.Vertex_elel_ndof[j]
                ))
    if isMC:
        #print out the gen info too!
        if gentry.nGenPart > 0:
            print("GenParts\n #   pdgId    pT    eta     phi     pz    vxy    vz    mass")
            for j in range(gentry.nGenPart):
                print("{0:d}      {1:d}    {2:.3f}    {3:.3f}    {4:.3f}    {5:.3f}    {6:.3f}    {7:.3f}   {8:.3f}".format(
                    j, gentry.GenPart_pdgId[j], gentry.GenPart_pt[j], gentry.GenPart_eta[j], gentry.GenPart_phi[j], gentry.GenPart_pz[j],
                    gentry.GenPart_vxy[j], gentry.GenPart_vz[j], gentry.GenPart_mass[j]
                    ))
    return
