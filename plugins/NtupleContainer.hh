#ifndef NTUPLECONTAINER_HH
#define NTUPLECONTAINER_HH

#include <vector>
#include <string>
#include <sstream>
using std::vector;

#include <TTree.h>

class NtupleContainer {

public:
    NtupleContainer();
    //list of photonID names to save
    vector<std::string> photonIDs = {"mvaPhoID-RunIIFall17-v2-wp80", "mvaPhoID-RunIIFall17-v2-wp90" };
    //list of electron ID names to save
    vector<std::string> electronIDs = {"mvaEleID-Fall17-noIso-V2-wp80", "mvaEleID-Fall17-noIso-V2-wp90", "mvaEleID-Fall17-noIso-V2-wpLoose" };
    virtual ~NtupleContainer();
    void SetRecoTree(TTree *tree);
    void SetGenTree(TTree *tree);
    void CreateTreeBranches();
    void ClearTreeBranches();

    // Trigger and event-level branches
    //first 64 trigger bits
    uint64_t fired0_;
    //last 8 trigger bits (but really only 2 bits used)
    uint8_t fired1_;
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
    vector<int8_t> genCharge_;
    vector<float> genPt_;
    vector<float> genEta_;
    vector<float> genPhi_;
    vector<float> genPz_;
    vector<float> genEn_;
    vector<float> genVxy_;
    vector<float> genVz_;
    vector<float> genMass_;
    
    // Reco branches
  
    // Reco electron branches
    uint8_t recoNGoodElectron_;
    vector<float> recoElectronPt_;
    vector<float> recoElectronEta_;
    vector<float> recoElectronPhi_;
    vector<float> recoElectronVxy_;
    vector<float> recoElectronVz_;
    vector<int8_t> recoElectronCharge_;
    vector<uint8_t> recoElectronIDResult_;
    
    // Reco photon branches
    uint8_t recoNGoodPhoton_;
    vector<float> recoPhotonPt_;
    vector<float> recoPhotonEta_;
    vector<float> recoPhotonPhi_;
    vector<uint8_t> recoPhotonIDResult_;

    // Reco muon branches
    uint8_t recoNGoodMuon_;
    vector<float> recoMuonPt_;
    vector<float> recoMuonEta_;
    vector<float> recoMuonPhi_;
    // muon ID is just 1 for all miniAOD muons...
    //vector<int> recoMuonIDResult_;
    vector<int8_t> recoMuonCharge_;
    vector<float> recoMuonVxy_;
    vector<float> recoMuonVz_;

    // Reco PackedCandidate track branches
    uint8_t recoNGoodTrk_;
    vector<float> recoTrkPt_;
    vector<float> recoTrkEta_;
    vector<float> recoTrkPhi_;
    vector<int8_t> recoTrkCharge_;

    // GsfElectron track branches
    //uint8_t gsfNGoodTrk_;
    vector<float> gsfTrkPt_;
    vector<float> gsfTrkEta_;
    vector<float> gsfTrkPhi_;
    vector<int8_t> gsfTrkCharge_;

    //mapping from the pos-only and neg-only list to the full list
    vector<uint8_t> mmeeTrxP;
    vector<uint8_t> mmeeTrxN;
    vector<uint8_t> gsfElsP;
    vector<uint8_t> gsfElsN;
    vector<uint8_t> muonsP;
    vector<uint8_t> muonsN;

    // Vertex branches
    float pvx_;
    float pvy_;
    float pvz_;

    //the maximum number of mmee vertices to save (if -1 no maximum)
    //int maxNmmee = -1;

    vector<float> mmeeVtxVxy_;
    vector<float> mmeeVtxVz_;
    vector<float> mmeeVtxSigmaVxy_;
    vector<float> mmeeVtxReducedChi2_;
    vector<float> mmeeVtxDr_;
    vector<float> mmeeVtxPt_;
    vector<float> mmeeVtxM_;
    vector<float> mmeeVtxM2_;
    vector<float> mmeeVtxPt2_;
    vector<uint8_t> mmeeVtxMuonP_;
    vector<uint8_t> mmeeVtxMuonN_;
    vector<uint8_t> mmeeVtxTrackP_;
    vector<uint8_t> mmeeVtxTrackN_;

    vector<float> mmelelVtxVxy_;
    vector<float> mmelelVtxVz_;
    vector<float> mmelelVtxSigmaVxy_;
    vector<float> mmelelVtxReducedChi2_;
    vector<float> mmelelVtxDr_;
    vector<float> mmelelVtxPt_;
    vector<float> mmelelVtxM_;
    vector<float> mmelelVtxM2_;
    vector<float> mmelelVtxPt2_;
    vector<uint8_t> mmelelVtxMuonP_;
    vector<uint8_t> mmelelVtxMuonN_;
    vector<uint8_t> mmelelVtxEleP_;
    vector<uint8_t> mmelelVtxEleN_;

    vector<float> elelVtxVxy_;
    vector<float> elelVtxVz_;
    vector<float> elelVtxSigmaVxy_;
    vector<float> elelVtxReducedChi2_;
    vector<float> elelVtxDr_;
    vector<float> elelVtxPt_;
    vector<float> elelVtxM_;
    vector<float> elelVtxM2_;
    vector<float> elelVtxPt2_;
    vector<uint8_t> elelVtxEleP_;
    vector<uint8_t> elelVtxEleN_;
    vector<uint8_t> elelVtxTrackP_;
    vector<uint8_t> elelVtxTrackN_;

    vector<float> pcpcVtxVxy_;
    vector<float> pcpcVtxVz_;
    vector<float> pcpcVtxSigmaVxy_;
    vector<float> pcpcVtxReducedChi2_;
    vector<float> pcpcVtxDr_;
    vector<float> pcpcVtxPt_;
    vector<float> pcpcVtxM_;
    vector<float> pcpcVtxM2_;
    vector<float> pcpcVtxPt2_;
    vector<uint8_t> pcpcVtxTrackP_;
    vector<uint8_t> pcpcVtxTrackN_;

    vector<float> mumuVtxVxy_;
    vector<float> mumuVtxVz_;
    vector<float> mumuVtxSigmaVxy_;
    vector<float> mumuVtxReducedChi2_;
    vector<float> mumuVtxDr_;
    vector<float> mumuVtxPt_;
    vector<float> mumuVtxM_;
    vector<float> mumuVtxM2_;
    vector<float> mumuVtxPt2_;
    vector<uint8_t> mumuVtxMuonP_;
    vector<uint8_t> mumuVtxMuonN_;

private:
    // Reco and gen TTrees
    TTree * recoT;
    TTree * genT;
    bool isData_;

};


#endif // NTUPLECONTAINER_HH
