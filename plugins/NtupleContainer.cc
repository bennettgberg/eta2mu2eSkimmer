#include "NtupleContainer.hh"

NtupleContainer::NtupleContainer() : isData_(true) {}

NtupleContainer::~NtupleContainer() {}

void NtupleContainer::SetRecoTree(TTree *tree) { recoT = tree; }
void NtupleContainer::SetGenTree(TTree *tree) { genT = tree; isData_ = false; }

void NtupleContainer::CreateTreeBranches() {

    //initialize vertex maps
    for ( std::string type : vtxTypes ) {
        recoVtxVxy_[type] = {};
        recoVtxVz_[type] = {};
        recoVtxSigmaVxy_[type] = {};
        recoVtxReducedChi2_[type] = {};
        recoVtxDr_[type] = {};
        recoVtxMuons_[type] = {};
        recoVtxEles_[type] = {};
        recoVtxTracks_[type] = {};
        mmeeTrxP = {};
        mmeeTrxN = {};
        gsfElsP = {};
        gsfElsN = {};
    }

    recoT->Branch("evt", &eventNum_);
    recoT->Branch("lumi_sec", &lumiSec_);
    recoT->Branch("run", &runNum_);
    recoT->Branch("nPV", &npv_);
    recoT->Branch("PV_vx", &pvx_);
    recoT->Branch("PV_vy", &pvy_);
    recoT->Branch("PV_vz", &pvz_);

    //vertices
    for( std::string type : vtxTypes ) {
        std::stringstream ssvxy;
        ssvxy << "Vertex_" << type << "_vxy";
        recoT->Branch(ssvxy.str().c_str(), &(recoVtxVxy_[type]));
        std::stringstream ssvz;
        ssvz << "Vertex_" << type << "_vz";
        recoT->Branch(ssvz.str().c_str(), &(recoVtxVz_[type]));
        std::stringstream sssigmavxy;
        sssigmavxy << "Vertex_" << type << "_sigmavxy";
        recoT->Branch(sssigmavxy.str().c_str(), &(recoVtxSigmaVxy_[type]));
        std::stringstream ssrchi2;
        ssrchi2 << "Vertex_" << type << "_reduced_chi2";
        recoT->Branch(ssrchi2.str().c_str(), &(recoVtxReducedChi2_[type]));
        std::stringstream ssdR;
        ssdR << "Vertex_" << type << "_dR";
        recoT->Branch(ssdR.str().c_str(), &(recoVtxDr_[type]));
        if( type == "mmee" ) {
            std::stringstream ssmuons;
            ssmuons << "Vertex_" << type << "_muons";
            recoT->Branch(ssmuons.str().c_str(), &recoVtxMuons_[type]);
        }
        if ( type == "elel" ) {
            std::stringstream ssels;
            ssels << "Vertex_" << type << "_eles";
            recoT->Branch(ssels.str().c_str(), &recoVtxEles_[type]);
        }
        if ( type == "mmee" || type == "pcpc" ) {
            std::stringstream sstrx;
            sstrx << "Vertex_" << type << "_tracks";
            recoT->Branch(sstrx.str().c_str(), &recoVtxTracks_[type]);
        }
    }

    recoT->Branch("nGoodElectron", &recoNGoodElectron_);
    recoT->Branch("Electron_pt",  &recoElectronPt_);
    recoT->Branch("Electron_eta", &recoElectronEta_);
    recoT->Branch("Electron_phi", &recoElectronPhi_);
    recoT->Branch("Electron_vxy", &recoElectronVxy_);
    recoT->Branch("Electron_vz",  &recoElectronVz_);
    recoT->Branch("Electron_charge", &recoElectronCharge_);
    recoT->Branch("Electron_id", &recoElectronIDResult_);

    recoT->Branch("nGoodPhoton", &recoNGoodPhoton_);
    recoT->Branch("Photon_pt",  &recoPhotonPt_);
    recoT->Branch("Photon_eta", &recoPhotonEta_);
    recoT->Branch("Photon_phi", &recoPhotonPhi_);
    recoT->Branch("Photon_id", &recoPhotonIDResult_);

    //muons :)
    recoT->Branch("nGoodMuon", &recoNGoodMuon_);
    recoT->Branch("Muon_pt",  &recoMuonPt_);
    recoT->Branch("Muon_eta", &recoMuonEta_);
    recoT->Branch("Muon_phi", &recoMuonPhi_);
    recoT->Branch("Muon_charge", &recoMuonCharge_);
    recoT->Branch("Muon_vxy", &recoMuonVxy_);
    recoT->Branch("Muon_vz", &recoMuonVz_);
    //recoT->Branch("Muon_id", &recoMuonIDResult_);

    //PackedCandidate tracks
    recoT->Branch("nGoodTrack", &recoNGoodTrk_);
    recoT->Branch("Track_pt", &recoTrkPt_);
    recoT->Branch("Track_eta", &recoTrkEta_);
    recoT->Branch("Track_phi", &recoTrkPhi_);
    recoT->Branch("Track_charge", &recoTrkCharge_);

    //gsfElectron tracks
    //there's exactly one GsfTrack for every electron so nGoodGsfTrack == nGoodElectron
    //recoT->Branch("nGoodGsfTrack", &gsfNGoodTrk_);
    recoT->Branch("GsfTrack_pt", &gsfTrkPt_);
    recoT->Branch("GsfTrack_eta", &gsfTrkEta_);
    recoT->Branch("GsfTrack_phi", &gsfTrkPhi_);
    recoT->Branch("GsfTrack_charge", &gsfTrkCharge_);

    if (!isData_) {
        genT->Branch("evt", &eventNum_);
        genT->Branch("nPU", &genpuobs_);
        genT->Branch("nPUtrue", &genputrue_);
        //genT->Branch("Generator_weight", &genwgt_);
        genT->Branch("nGenPart", &nGen_);
        genT->Branch("GenPart_pdgId", &genID_);
        genT->Branch("GenPart_hardProcess", &genHardProcess_);
        genT->Branch("GenPart_charge", &genCharge_);
        genT->Branch("GenPart_pt", &genPt_);
        genT->Branch("GenPart_eta", &genEta_);
        genT->Branch("GenPart_phi", &genPhi_);
        genT->Branch("GenPart_pz", &genPz_);
        genT->Branch("GenPart_energy", &genEn_);
        genT->Branch("GenPart_vxy", &genVxy_);
        genT->Branch("GenPart_vz", &genVz_);
        genT->Branch("GenPart_mass", &genMass_);
    }

}

void NtupleContainer::ClearTreeBranches() {

    for(std::string type : vtxTypes) {
        recoVtxVxy_[type].clear();
        recoVtxVz_[type].clear();
        recoVtxSigmaVxy_[type].clear();
        recoVtxReducedChi2_[type].clear();
        recoVtxDr_[type].clear();
        recoVtxMuons_[type].clear();
        recoVtxEles_[type].clear();
        recoVtxTracks_[type].clear();
    }

    recoElectronPt_.clear();
    recoElectronEta_.clear();
    recoElectronPhi_.clear();
    recoElectronVxy_.clear();
    recoElectronVz_.clear();
    recoElectronCharge_.clear();
    recoElectronIDResult_.clear();
    recoPhotonPt_.clear();
    recoPhotonEta_.clear();
    recoPhotonPhi_.clear();
    recoPhotonIDResult_.clear();
    recoMuonPt_.clear();
    recoMuonEta_.clear();
    recoMuonPhi_.clear();
    recoMuonCharge_.clear();
    recoMuonVxy_.clear();
    recoMuonVz_.clear();
    //recoMuonIDResult_.clear();

    //recoNElectron_ = 0;
    recoNGoodElectron_ = 0;
    //recoNPhoton_ = 0;
    recoNGoodPhoton_ = 0;

    //fired_ = 0;

    genID_.clear();
    genHardProcess_.clear();
    genCharge_.clear();
    genPt_.clear();
    genEta_.clear();
    genPhi_.clear();
    genPz_.clear();
    genEn_.clear();
    genVxy_.clear();
    genVz_.clear();
    genMass_.clear();

    // PackedCandidate tracks    
    recoTrkPt_.clear();
    recoTrkEta_.clear();
    recoTrkPhi_.clear();
    recoTrkCharge_.clear();

    // GsfElectron tracks    
    gsfTrkPt_.clear();
    gsfTrkEta_.clear();
    gsfTrkPhi_.clear();
    gsfTrkCharge_.clear();

    // Pile-up and event genweight
    genpuobs_ = -9999;
    genputrue_ = -9999;
    //genwgt_ = -9999;


}
