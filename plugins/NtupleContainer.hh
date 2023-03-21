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
    NtupleContainer();
    //types of vertices
    // elec-elec, packedCand-packedCand, mu-mu-e-e (where e's are packedCandidates)
    vector<std::string> vtxTypes = { "elel", "pcpc", "mmee" };
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

    // Vertex branches
    float pvx_;
    float pvy_;
    float pvz_;
    std::map<std::string, vector<float>> recoVtxVxy_;
    std::map<std::string, vector<float>> recoVtxVz_;
    std::map<std::string, vector<float>> recoVtxSigmaVxy_;
    std::map<std::string, vector<float>> recoVtxReducedChi2_;
    std::map<std::string, vector<float>> recoVtxDr_;
    //std::map<std::string, vector<float>> recoVtxPt_;
    //std::map<std::string, vector<float>> recoVtxM_;
    std::map<std::string, vector<uint8_t>> recoVtxMuons_;
    std::map<std::string, vector<uint8_t>> recoVtxEles_;
    std::map<std::string, vector<uint8_t>> recoVtxTracks_;
    //which number in the FULL list of tracks (pos AND neg) is this (Pos OR Neg) track?
    vector<uint8_t> mmeeTrxP;
    vector<uint8_t> mmeeTrxN;
    vector<uint8_t> gsfElsP;
    vector<uint8_t> gsfElsN;

private:
    // Reco and gen TTrees
    TTree * recoT;
    TTree * genT;
    bool isData_;

};


#endif // NTUPLECONTAINER_HH
