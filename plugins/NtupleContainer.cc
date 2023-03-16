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
        recoVtxPt_[type] = {};
        recoVtxM_[type] = {};
    }

    for( std::string id : photonIDs ) {
        recoPhotonIDResult_[id] = {};
    }

    for( std::string id : electronIDs ) {
        recoElectronIDResult_[id] = {};
    }

    recoT->Branch("evt", &eventNum_);
    recoT->Branch("lumi_sec", &lumiSec_);
    recoT->Branch("run", &runNum_);
    recoT->Branch("nPV", &npv_);
    //recoT->Branch("MET_filters_fail_bits", &METFiltersFailBits_);
    //recoT->Branch("trig_fired", &fired_);
    //recoT->Branch("reco_n_gm", &recoNGM_);
    //recoT->Branch("reco_n_good_gm", &recoNGoodGM_);
    //recoT->Branch("reco_gm_pt",  &recoGMPt_);
    //recoT->Branch("reco_gm_pt_err",  &recoGMPtError_);
    //recoT->Branch("reco_gm_eta", &recoGMEta_);
    //recoT->Branch("reco_gm_eta_err", &recoGMEtaError_);
    //recoT->Branch("reco_gm_phi", &recoGMPhi_);
    //recoT->Branch("reco_gm_phi_err", &recoGMPhiError_);
    //recoT->Branch("reco_gm_dxy", &recoGMDxy_);
    //recoT->Branch("reco_gm_dxy_err", &recoGMDxyError_);
    //recoT->Branch("reco_gm_dz",  &recoGMDz_);
    //recoT->Branch("reco_gm_dz_err",  &recoGMDzError_);
    //recoT->Branch("reco_gm_charge", &recoGMCharge_);
    //recoT->Branch("reco_gm_trk_chi2", &recoGMTrkChi2_);
    //recoT->Branch("reco_gm_trk_n_planes", &recoGMTrkNumPlanes_);
    //recoT->Branch("reco_gm_trk_n_hits", &recoGMTrkNumHits_);
    //recoT->Branch("reco_gm_trk_n_DT_hits", &recoGMTrkNumDTHits_);
    //recoT->Branch("reco_gm_trk_n_CSC_hits", &recoGMTrkNumCSCHits_);
    //recoT->Branch("reco_gm_isPF", &recoGMIsPF_);
    //recoT->Branch("reco_gm_PFIso", &recoGMPFIso_);
    //recoT->Branch("reco_gm_TrkIso", &recoGMTrkIso_);
    //recoT->Branch("reco_gm_trk_n_pix_hits", &recoGMTrkNumPixelHit_);
    //recoT->Branch("reco_gm_trk_n_trk_layers", &recoGMTrkNumTrkLayers_);
    //recoT->Branch("reco_Mmumu",  &recoMmumu_);
    recoT->Branch("PV_vx", &pvx_);
    recoT->Branch("PV_vy", &pvy_);
    recoT->Branch("PV_vz", &pvz_);
    //recoT->Branch("reco_vtx_gmgm_vxy", &gmgm_recoVtxVxy_);
    //recoT->Branch("reco_vtx_gmgm_vz",  &gmgm_recoVtxVz_);
    //recoT->Branch("reco_vtx_gmgm_sigmavxy", &gmgm_recoVtxSigmaVxy_);
    //recoT->Branch("reco_vtx_gmgm_reduced_chi2", &gmgm_recoVtxReducedChi2_);
    //recoT->Branch("reco_vtx_gmgm_dR",  &gmgm_recoVtxDr_);
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
        if( type == "mmee" || type == "mmpp" ) {
            //pt of the reco'd parent particle at this vertex
            std::stringstream sspT;
            sspT << "Vertex_" << type << "_pt";
            recoT->Branch(sspT.str().c_str(), &(recoVtxPt_[type]));
            //mass of the reco'd parent particle at this vertex
            std::stringstream ssm;
            ssm << "Vertex_" << type << "_m";
            recoT->Branch(ssm.str().c_str(), &(recoVtxM_[type]));
        }
    }
    ////packed candidate - packed candidate vertex
    //recoT->Branch("Vertex_pcpc_vxy", &pcpc_recoVtxVxy_);
    //recoT->Branch("Vertex_pcpc_vz",  &pcpc_recoVtxVz_);
    //recoT->Branch("Vertex_pcpc_sigmavxy", &pcpc_recoVtxSigmaVxy_);
    //recoT->Branch("Vertex_pcpc_reduced_chi2", &pcpc_recoVtxReducedChi2_);
    //recoT->Branch("Vertex_pcpc_dR",  &pcpc_recoVtxDr_);
    ////----
    recoT->Branch("nElectron", &recoNElectron_);
    recoT->Branch("nGoodElectron", &recoNGoodElectron_);
    recoT->Branch("Electron_pt",  &recoElectronPt_);
    recoT->Branch("Electron_eta", &recoElectronEta_);
    recoT->Branch("Electron_phi", &recoElectronPhi_);
    recoT->Branch("Electron_vxy", &recoElectronVxy_);
    recoT->Branch("Electron_vz",  &recoElectronVz_);
    recoT->Branch("Electron_charge", &recoElectronCharge_);
    //recoT->Branch("Electron_id", &recoElectronIDResult_);
    int idnum = 0;
    for(std::string id : electronIDs) {
        std::stringstream ss;
        ss << "Electron_id" << idnum;
        recoT->Branch(ss.str().c_str(), &recoElectronIDResult_[id]); 
        idnum++;
    }
    recoT->Branch("nPhoton", &recoNPhoton_);
    recoT->Branch("nGoodPhoton", &recoNGoodPhoton_);
    recoT->Branch("Photon_pt",  &recoPhotonPt_);
    recoT->Branch("Photon_eta", &recoPhotonEta_);
    recoT->Branch("Photon_phi", &recoPhotonPhi_);
    //recoT->Branch("Photon_id", &recoPhotonIDResult_);
    idnum = 0;
    for(std::string id : photonIDs) {
        std::stringstream ss;
        ss << "Photon_id" << idnum;
        recoT->Branch(ss.str().c_str(), &recoPhotonIDResult_[id]); 
        idnum++;
    }
    //muons :)
    recoT->Branch("nMuon", &recoNMuon_);
    recoT->Branch("nGoodMuon", &recoNGoodMuon_);
    recoT->Branch("Muon_pt",  &recoMuonPt_);
    recoT->Branch("Muon_eta", &recoMuonEta_);
    recoT->Branch("Muon_phi", &recoMuonPhi_);
    recoT->Branch("Muon_charge", &recoMuonCharge_);
    recoT->Branch("Muon_vxy", &recoMuonVxy_);
    recoT->Branch("Muon_vz", &recoMuonVz_);
    recoT->Branch("Muon_id", &recoMuonIDResult_);

    //PackedCandidate tracks
    recoT->Branch("nTrack", &recoNTrk_);
    recoT->Branch("nGoodTrack", &recoNGoodTrk_);
    recoT->Branch("Track_pt", &recoTrkPt_);
    recoT->Branch("Track_eta", &recoTrkEta_);
    recoT->Branch("Track_phi", &recoTrkPhi_);
    recoT->Branch("Track_charge", &recoTrkCharge_);

    //gsfElectron tracks
    recoT->Branch("nGsfTrack", &gsfNTrk_);
    recoT->Branch("nGoodGsfTrack", &gsfNGoodTrk_);
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

    //recoGMPt_.clear();
    //recoGMPtError_.clear();
    //recoGMEta_.clear();
    //recoGMEtaError_.clear();
    //recoGMPhi_.clear();
    //recoGMPhiError_.clear();
    //recoGMDxy_.clear();
    //recoGMDxyError_.clear();
    //recoGMDz_.clear();
    //recoGMDzError_.clear();
    //recoGMCharge_.clear();
    //recoGMTrkChi2_.clear();
    //recoGMTrkNumPlanes_.clear();
    //recoGMTrkNumHits_.clear();
    //recoGMTrkNumDTHits_.clear();
    //recoGMTrkNumCSCHits_.clear();
    //recoGMIsPF_.clear();
    //recoGMPFIso_.clear();
    //recoGMTrkIso_.clear();
    //recoGMTrkNumPixelHit_.clear();
    //recoGMTrkNumTrkLayers_.clear();
    //gmgm_recoVtxVxy_.clear();
    //gmgm_recoVtxVz_.clear();
    //gmgm_recoVtxSigmaVxy_.clear();
    //gmgm_recoVtxReducedChi2_.clear();
    //gmgm_recoVtxDr_.clear();
    for(std::string type : vtxTypes) {
        recoVtxVxy_[type].clear();
        recoVtxVz_[type].clear();
        recoVtxSigmaVxy_[type].clear();
        recoVtxReducedChi2_[type].clear();
        recoVtxDr_[type].clear();
        recoVtxPt_[type].clear();
        recoVtxM_[type].clear();
    }

    recoElectronPt_.clear();
    recoElectronEta_.clear();
    recoElectronPhi_.clear();
    recoElectronVxy_.clear();
    recoElectronVz_.clear();
    recoElectronCharge_.clear();
    //recoElectronIDResult_.clear();
    for(std::string id : electronIDs) {
        recoElectronIDResult_[id].clear();
    }
    recoPhotonPt_.clear();
    recoPhotonEta_.clear();
    recoPhotonPhi_.clear();
    //recoPhotonIDResult_.clear();
    for(std::string id : photonIDs) {
        recoPhotonIDResult_[id].clear();
    }
    recoMuonPt_.clear();
    recoMuonEta_.clear();
    recoMuonPhi_.clear();
    recoMuonCharge_.clear();
    recoMuonVxy_.clear();
    recoMuonVz_.clear();
    recoMuonIDResult_.clear();

    recoNElectron_ = 0;
    recoNGoodElectron_ = 0;
    recoNPhoton_ = 0;
    recoNGoodPhoton_ = 0;

    //fired_ = 0;

    //recoMmumu_ = -9999;

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
