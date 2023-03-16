#ifndef NTUPLECONTAINER_HH
#define NTUPLECONTAINER_HH

#include <vector>
#include <map>
#include <string>
#include <sstream>
using std::vector;

#include <TTree.h>

class NtupleContainer {

public:
    //types of vertices
    // elec-elec, packedCand-packedCand, mu-mu-e-e (where e's are packedCandidates)
    NtupleContainer();
    vector<std::string> vtxTypes = { "elel", "pcpc", "mmee", "mmpp" };
    vector<std::string> photonIDs = {"mvaPhoID-RunIIFall17-v2-wp80", "mvaPhoID-RunIIFall17-v2-wp90" };
    vector<std::string> electronIDs = {"mvaEleID-Fall17-noIso-V2-wp80", "mvaEleID-Fall17-noIso-V2-wp90", "mvaEleID-Fall17-noIso-V2-wpLoose" };
    virtual ~NtupleContainer();
    void SetRecoTree(TTree *tree);
    void SetGenTree(TTree *tree);
    void CreateTreeBranches();
    void ClearTreeBranches();

    // Trigger and event-level branches
    //unsigned int fired_;
    unsigned long long eventNum_;
    unsigned long long runNum_;
    unsigned long long lumiSec_;

    //uint32_t METFiltersFailBits_;

    // Gen branches
    
    // Pileup and gen weight
    int genpuobs_;
    int genputrue_;
    //float genwgt_;
    int npv_;

    // Gen particles
    int nGen_;
    vector<int> genID_;
    // Only save hard-process gen particles--haha jk
    vector<bool> genHardProcess_; //--> would always be 1  ??? 0??
    vector<int> genCharge_;
    vector<float> genPt_;
    vector<float> genEta_;
    vector<float> genPhi_;
    vector<float> genPz_;
    vector<float> genEn_;
    vector<float> genVxy_;
    vector<float> genVz_;
    vector<float> genMass_;
    
    // Reco branches
 
   // // Reco GM branches
   // int recoNGM_;
   // int recoNGoodGM_;
   // vector<float> recoGMPt_;
   // vector<float> recoGMPtError_;
   // vector<float> recoGMEta_;
   // vector<float> recoGMEtaError_;
   // vector<float> recoGMPhi_;
   // vector<float> recoGMPhiError_;
   // vector<float> recoGMDxy_;
   // vector<float> recoGMDxyError_;
   // vector<float> recoGMDz_;
   // vector<float> recoGMDzError_;
   // vector<int> recoGMCharge_;
   // vector<float> recoGMTrkChi2_;
   // vector<float> recoGMTrkNumPlanes_;
   // vector<float> recoGMTrkNumHits_;
   // vector<float> recoGMTrkNumDTHits_;
   // vector<float> recoGMTrkNumCSCHits_;
   // vector<bool> recoGMIsPF_;
   // vector<float> recoGMPFIso_;
   // vector<float> recoGMTrkIso_;
   // vector<float> recoGMTrkNumPixelHit_;
   // vector<float> recoGMTrkNumTrkLayers_;
    
    // Reco electron branches
    int recoNElectron_;
    int recoNGoodElectron_;
    vector<float> recoElectronPt_;
    vector<float> recoElectronEta_;
    vector<float> recoElectronPhi_;
    vector<float> recoElectronVxy_;
    vector<float> recoElectronVz_;
    vector<int> recoElectronCharge_;
    //vector<float> recoElectronIDResult_;
    std::map<std::string, vector<float>> recoElectronIDResult_; 
    
    // Reco photon branches
    int recoNPhoton_;
    int recoNGoodPhoton_;
    vector<float> recoPhotonPt_;
    vector<float> recoPhotonEta_;
    vector<float> recoPhotonPhi_;
    //vector<float> recoPhotonIDResult_;
    std::map<std::string, vector<float>> recoPhotonIDResult_;

    // Reco muon branches
    int recoNMuon_;
    int recoNGoodMuon_;
    vector<float> recoMuonPt_;
    vector<float> recoMuonEta_;
    vector<float> recoMuonPhi_;
    vector<int> recoMuonIDResult_;
    vector<int> recoMuonCharge_;
    vector<float> recoMuonVxy_;
    vector<float> recoMuonVz_;

    // Reco PackedCandidate track branches
    int recoNTrk_;
    int recoNGoodTrk_;
    vector<float> recoTrkPt_;
    vector<float> recoTrkEta_;
    vector<float> recoTrkPhi_;
    vector<int> recoTrkCharge_;

    // GsfElectron track branches
    int gsfNTrk_;
    int gsfNGoodTrk_;
    vector<float> gsfTrkPt_;
    vector<float> gsfTrkEta_;
    vector<float> gsfTrkPhi_;
    vector<float> gsfTrkCharge_;

    // Vertex branches
    float pvx_;
    float pvy_;
    float pvz_;
    //float recoVtxVxy_;
    //float recoVtxVz_;
    //float recoVtxSigmaVxy_;
    //float recoVtxReducedChi2_;
    //float recoVtxDr_;
    vector<float> gmgm_recoVtxVxy_;
    vector<float> gmgm_recoVtxVz_;
    vector<float> gmgm_recoVtxSigmaVxy_;
    vector<float> gmgm_recoVtxReducedChi2_;
    vector<float> gmgm_recoVtxDr_;
    std::map<std::string, vector<float>> recoVtxVxy_;
    std::map<std::string, vector<float>> recoVtxVz_;
    std::map<std::string, vector<float>> recoVtxSigmaVxy_;
    std::map<std::string, vector<float>> recoVtxReducedChi2_;
    std::map<std::string, vector<float>> recoVtxDr_;
    std::map<std::string, vector<float>> recoVtxPt_;
    std::map<std::string, vector<float>> recoVtxM_;
    ////electron vertices
    //vector<float> elel_recoVtxVxy_;
    //vector<float> elel_recoVtxVz_;
    //vector<float> elel_recoVtxSigmaVxy_;
    //vector<float> elel_recoVtxReducedChi2_;
    //vector<float> elel_recoVtxDr_;
    ////packed candidate vertices
    //vector<float> pcpc_recoVtxVxy_;
    //vector<float> pcpc_recoVtxVz_;
    //vector<float> pcpc_recoVtxSigmaVxy_;
    //vector<float> pcpc_recoVtxReducedChi2_;
    //vector<float> pcpc_recoVtxDr_;

    //float recoMmumu_;

private:
    // Reco and gen TTrees
    TTree * recoT;
    TTree * genT;
    bool isData_;

};


#endif // NTUPLECONTAINER_HH
