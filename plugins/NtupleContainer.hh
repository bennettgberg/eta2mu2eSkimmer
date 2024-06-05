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
    //separate list for low pT electrons
    vector<std::string> lowPtElectronIDs = {"ID", "ptbiased", "unbiased" };
    //remove from consideration for vertexing purposes any electron matched to an Onia converted photon? or nah?
    bool removeOnia = false; // true;
    //list of LowPtElectrons to skip whilst doing the vertex computations (because they are too similar to an Onia converted photon)
    vector<int> skipListP;
    vector<int> skipListN;
    virtual ~NtupleContainer();
    void SetRecoTree(TTree *tree);
    void SetGenTree(TTree *tree);
    void CreateTreeBranches();
    void ClearTreeBranches();

    // Trigger and event-level branches
    //first 64 trigger bits
    uint32_t fired0_;
    //last 8 trigger bits (but really only 2 bits used)
    uint32_t fired1_;
    uint8_t fired2_;
    unsigned long long eventNum_;
    unsigned long long runNum_;
    unsigned long long lumiSec_;

    unsigned long long MClumiblock_;
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
  
    //photon conversion candidates
    uint8_t recoNOnia_;
    //info for the first track for the converted photon
    vector<float> recoOniaPt0_;
    vector<float> recoOniaEta0_;
    vector<float> recoOniaPhi0_;
    vector<float> recoOniaCharge0_;
    //infor for the second track for the converted photon
    vector<float> recoOniaPt1_;
    vector<float> recoOniaEta1_;
    vector<float> recoOniaPhi1_;
    vector<float> recoOniaCharge1_;
    //vector<float> recoOniaVtxVxy_;
    vector<float> recoOniaVtxVx_;
    vector<float> recoOniaVtxVy_;
    vector<float> recoOniaVtxVz_;

    // Reco electron branches
    uint8_t recoNGoodElectron_;
    vector<float> recoElectronPt_;
    vector<float> recoElectronEta_;
    vector<float> recoElectronPhi_;
    vector<float> recoElectronVxy_;
    vector<float> recoElectronVz_;
    vector<int8_t> recoElectronCharge_;
    vector<uint8_t> recoElectronIDResult_;
    //conversion veto
    vector<uint8_t> recoElectronConvVeto_;
    //number of missing layers in tracker hits
    vector<uint8_t> recoElectronNMhits_;

    // Reco electron branches
    uint8_t recoNGoodLowPtElectron_;
    vector<float> recoLowPtElectronPt_;
    vector<float> recoLowPtElectronEta_;
    vector<float> recoLowPtElectronPhi_;
    vector<float> recoLowPtElectronVxy_;
    vector<float> recoLowPtElectronVz_;
    vector<int8_t> recoLowPtElectronCharge_;
    vector<uint8_t> recoLowPtElectronIDResult_;
    
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
    vector<int8_t> recoMuonIDResult_;
    vector<int8_t> recoMuonCharge_;
    vector<float> recoMuonVxy_;
    vector<float> recoMuonVz_;

    // Reco PackedCandidate track branches
    uint8_t recoNGoodTrk_;
    vector<float> recoTrkPt_;
    vector<float> recoTrkEta_;
    vector<float> recoTrkPhi_;
    vector<int8_t> recoTrkCharge_;
    vector<float> recoTrkDxy_;
    vector<float> recoTrkDz_;

    // GsfElectron track branches
    //uint8_t gsfNGoodTrk_;
    vector<float> gsfTrkPt_;
    vector<float> gsfTrkEta_;
    vector<float> gsfTrkPhi_;
    vector<int8_t> gsfTrkCharge_;
    vector<float> gsfTrkDxy_;
    vector<float> gsfTrkDz_;

    vector<float> gsfLowPtTrkPt_;
    vector<float> gsfLowPtTrkEta_;
    vector<float> gsfLowPtTrkPhi_;
    vector<int8_t> gsfLowPtTrkCharge_;
    vector<float> gsfLowPtTrkDxy_;
    vector<float> gsfLowPtTrkDz_;

    //mapping from the pos-only and neg-only list to the full list
    vector<uint8_t> mmeeTrxP;
    vector<uint8_t> mmeeTrxN;
    vector<uint8_t> gsfElsP;
    vector<uint8_t> gsfElsN;
    vector<uint8_t> gsfLowPtElsP;
    vector<uint8_t> gsfLowPtElsN;
    vector<uint8_t> muonsP;
    vector<uint8_t> muonsN;

    // Vertex branches
    float pvx_;
    float pvy_;
    float pvz_;

    //the maximum number of mmee vertices to save (if -1 no maximum)
    //int maxNmmee = -1;

    //vector<float> mmeeVtxVxy_;
    vector<float> mmeeVtxVx_;
    vector<float> mmeeVtxVy_;
    vector<float> mmeeVtxVz_;
    vector<float> mmeeVtxSigmaVxy_;
    //vector<float> mmeeVtxReducedChi2_;
    vector<float> mmeeVtxChi2_;
    vector<int> mmeeVtxNdof_;
    vector<float> mmeeVtxDr_;
    vector<float> mmeeVtxPt_;
    vector<float> mmeeVtxM_;
    vector<uint8_t> mmeeVtxMuonP_;
    vector<uint8_t> mmeeVtxMuonN_;
    vector<uint8_t> mmeeVtxTrackP_;
    vector<uint8_t> mmeeVtxTrackN_;

    //vector<float> mmelelVtxVxy_;
    vector<float> mmelelVtxVx_;
    vector<float> mmelelVtxVy_;
    vector<float> mmelelVtxVz_;
    vector<float> mmelelVtxSigmaVxy_;
    //vector<float> mmelelVtxReducedChi2_;
    vector<float> mmelelVtxChi2_;
    vector<int> mmelelVtxNdof_;
    vector<float> mmelelVtxDr_;
    vector<float> mmelelVtxPt_;
    vector<float> mmelelVtxM_;
    vector<uint8_t> mmelelVtxMuonP_;
    vector<uint8_t> mmelelVtxMuonN_;
    vector<uint8_t> mmelelVtxEleP_;
    vector<uint8_t> mmelelVtxEleN_;

    //mu-mu-lowPtElectron-lowPtElectron vertices
    //vector<float> mmlplpVtxVxy_;
    vector<float> mmlplpVtxVx_;
    vector<float> mmlplpVtxVy_;
    vector<float> mmlplpVtxVz_;
    vector<float> mmlplpVtxSigmaVxy_;
    //vector<float> mmlplpVtxReducedChi2_;
    vector<float> mmlplpVtxChi2_;
    vector<int> mmlplpVtxNdof_;
    vector<float> mmlplpVtxDr_;
    vector<float> mmlplpVtxPt_;
    vector<float> mmlplpVtxM_;
    vector<uint8_t> mmlplpVtxMuonP_;
    vector<uint8_t> mmlplpVtxMuonN_;
    vector<uint8_t> mmlplpVtxEleP_;
    vector<uint8_t> mmlplpVtxEleN_;

    //lowPtElectron-lowPtElectron vertices
    //vector<float> lplpVtxVxy_;
    vector<float> lplpVtxVx_;
    vector<float> lplpVtxVy_;
    vector<float> lplpVtxVz_;
    vector<float> lplpVtxSigmaVxy_;
    //vector<float> lplpVtxReducedChi2_;
    vector<float> lplpVtxChi2_;
    vector<int> lplpVtxNdof_;
    vector<float> lplpVtxDr_;
    vector<float> lplpVtxPt_;
    vector<float> lplpVtxM_;
    vector<uint8_t> lplpVtxEleP_;
    vector<uint8_t> lplpVtxEleN_;
    vector<uint8_t> lplpVtxTrackP_;
    vector<uint8_t> lplpVtxTrackN_;

    //vector<float> elelVtxVxy_;
    vector<float> elelVtxVx_;
    vector<float> elelVtxVy_;
    vector<float> elelVtxVz_;
    vector<float> elelVtxSigmaVxy_;
    //vector<float> elelVtxReducedChi2_;
    vector<float> elelVtxChi2_;
    vector<int> elelVtxNdof_;
    vector<float> elelVtxDr_;
    vector<float> elelVtxPt_;
    vector<float> elelVtxM_;
    vector<uint8_t> elelVtxEleP_;
    vector<uint8_t> elelVtxEleN_;
    vector<uint8_t> elelVtxTrackP_;
    vector<uint8_t> elelVtxTrackN_;

    //vector<float> pcpcVtxVxy_;
    vector<float> pcpcVtxVx_;
    vector<float> pcpcVtxVy_;
    vector<float> pcpcVtxVz_;
    vector<float> pcpcVtxSigmaVxy_;
    //vector<float> pcpcVtxReducedChi2_;
    vector<float> pcpcVtxChi2_;
    vector<int> pcpcVtxNdof_;
    vector<float> pcpcVtxDr_;
    vector<float> pcpcVtxPt_;
    vector<float> pcpcVtxM_;
    vector<uint8_t> pcpcVtxTrackP_;
    vector<uint8_t> pcpcVtxTrackN_;

    //vector<float> mumuVtxVxy_;
    vector<float> mumuVtxVx_;
    vector<float> mumuVtxVy_;
    vector<float> mumuVtxVz_;
    vector<float> mumuVtxSigmaVxy_;
    //vector<float> mumuVtxReducedChi2_;
    vector<float> mumuVtxChi2_;
    vector<int> mumuVtxNdof_;
    vector<float> mumuVtxDr_;
    vector<float> mumuVtxPt_;
    vector<float> mumuVtxM_;
    vector<uint8_t> mumuVtxMuonP_;
    vector<uint8_t> mumuVtxMuonN_;

private:
    // Reco and gen TTrees
    TTree * recoT;
    TTree * genT;
    bool isData_;

};


#endif // NTUPLECONTAINER_HH
