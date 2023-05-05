#include "NtupleContainer.hh"
#include <iostream>

NtupleContainer::NtupleContainer() : isData_(true) {}

NtupleContainer::~NtupleContainer() {}

void NtupleContainer::SetRecoTree(TTree *tree) { recoT = tree; }
void NtupleContainer::SetGenTree(TTree *tree) { genT = tree; isData_ = false; }

void NtupleContainer::CreateTreeBranches() {

    mmeeTrxP = {};
    mmeeTrxN = {};
    gsfElsP = {};
    gsfElsN = {};

    recoT->Branch("evt", &eventNum_);
    recoT->Branch("lumi_sec", &lumiSec_);
    recoT->Branch("run", &runNum_);
    recoT->Branch("nPV", &npv_);
    recoT->Branch("PV_vx", &pvx_);
    recoT->Branch("PV_vy", &pvy_);
    recoT->Branch("PV_vz", &pvz_);

    //vertices
    recoT->Branch("Vertex_mmee_vxy", &mmeeVtxVxy_);
    recoT->Branch("Vertex_mmee_vz", &mmeeVtxVz_);
    recoT->Branch("Vertex_mmee_sigmaVxy", &mmeeVtxSigmaVxy_);
    recoT->Branch("Vertex_mmee_reduced_chi2", &mmeeVtxReducedChi2_);
    recoT->Branch("Vertex_mmee_dR", &mmeeVtxDr_);
    recoT->Branch("Vertex_mmee_muP", &mmeeVtxMuonP_);
    recoT->Branch("Vertex_mmee_muN", &mmeeVtxMuonN_);
    recoT->Branch("Vertex_mmee_trackP", &mmeeVtxTrackP_);
    recoT->Branch("Vertex_mmee_trackN", &mmeeVtxTrackN_);
    //for debugging only, do NOT include in final ntuple (waste of space)
    //recoT->Branch("Vertex_mmee_M", &mmeeVtxM_);
    //recoT->Branch("Vertex_mmee_Pt", &mmeeVtxPt_);
    
    recoT->Branch("Vertex_mmelel_vxy", &mmelelVtxVxy_);
    recoT->Branch("Vertex_mmelel_vz", &mmelelVtxVz_);
    recoT->Branch("Vertex_mmelel_sigmaVxy", &mmelelVtxSigmaVxy_);
    recoT->Branch("Vertex_mmelel_reduced_chi2", &mmelelVtxReducedChi2_);
    recoT->Branch("Vertex_mmelel_dR", &mmelelVtxDr_);
    recoT->Branch("Vertex_mmelel_muP", &mmelelVtxMuonP_);
    recoT->Branch("Vertex_mmelel_muN", &mmelelVtxMuonN_);
    recoT->Branch("Vertex_mmelel_eleP", &mmelelVtxEleP_);
    recoT->Branch("Vertex_mmelel_eleN", &mmelelVtxEleN_);

    recoT->Branch("Vertex_elel_vxy", &elelVtxVxy_);
    recoT->Branch("Vertex_elel_vz", &elelVtxVz_);
    recoT->Branch("Vertex_elel_sigmaVxy", &elelVtxSigmaVxy_);
    recoT->Branch("Vertex_elel_reduced_chi2", &elelVtxReducedChi2_);
    recoT->Branch("Vertex_elel_dR", &elelVtxDr_);
    recoT->Branch("Vertex_elel_eleP", &elelVtxEleP_);
    recoT->Branch("Vertex_elel_eleN", &elelVtxEleN_);
    //for debugging only, do NOT include in final ntuple (waste of space)
    //recoT->Branch("Vertex_elel_M", &elelVtxM_);
    //recoT->Branch("Vertex_elel_Pt", &elelVtxPt_);

    recoT->Branch("Vertex_pcpc_vxy", &pcpcVtxVxy_);
    recoT->Branch("Vertex_pcpc_vz", &pcpcVtxVz_);
    recoT->Branch("Vertex_pcpc_sigmaVxy", &pcpcVtxSigmaVxy_);
    recoT->Branch("Vertex_pcpc_reduced_chi2", &pcpcVtxReducedChi2_);
    recoT->Branch("Vertex_pcpc_dR", &pcpcVtxDr_);
    recoT->Branch("Vertex_pcpc_trackP", &pcpcVtxTrackP_);
    recoT->Branch("Vertex_pcpc_trackN", &pcpcVtxTrackN_);
    //for debugging only, do NOT include in final ntuple (waste of space)
    //recoT->Branch("Vertex_pcpc_M", &pcpcVtxM_);
    //recoT->Branch("Vertex_pcpc_Pt", &pcpcVtxPt_);

    recoT->Branch("Vertex_mumu_vxy", &mumuVtxVxy_);
    recoT->Branch("Vertex_mumu_vz", &mumuVtxVz_);
    recoT->Branch("Vertex_mumu_sigmaVxy", &mumuVtxSigmaVxy_);
    recoT->Branch("Vertex_mumu_reduced_chi2", &mumuVtxReducedChi2_);
    recoT->Branch("Vertex_mumu_dR", &mumuVtxDr_);
    recoT->Branch("Vertex_mumu_muP", &mumuVtxMuonP_);
    recoT->Branch("Vertex_mumu_muN", &mumuVtxMuonN_);
    //for debugging only, do NOT include in final ntuple (waste of space)
    //recoT->Branch("Vertex_mumu_M", &mumuVtxM_);
    //recoT->Branch("Vertex_mumu_Pt", &mumuVtxPt_);

    //good particles
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
    //which triggers fired?
    recoT->Branch("Triggers_fired0", &fired0_);
    recoT->Branch("Triggers_fired1", &fired1_);

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

    mmeeVtxVxy_.clear();
    mmeeVtxVz_.clear();
    mmeeVtxSigmaVxy_.clear();
    mmeeVtxReducedChi2_.clear();
    mmeeVtxDr_.clear();
    //mmeeVtxM_.clear();
    //mmeeVtxPt_.clear();
    //mmeeVtxM2_.clear();
    //mmeeVtxPt2_.clear();
    mmeeVtxMuonP_.clear();
    mmeeVtxMuonN_.clear();
    mmeeVtxTrackP_.clear();
    mmeeVtxTrackN_.clear();

    mmelelVtxVxy_.clear();
    mmelelVtxVz_.clear();
    mmelelVtxSigmaVxy_.clear();
    mmelelVtxReducedChi2_.clear();
    mmelelVtxDr_.clear();
    mmelelVtxM_.clear();
    mmelelVtxPt_.clear();
    mmelelVtxM2_.clear();
    mmelelVtxPt2_.clear();
    mmelelVtxMuonP_.clear();
    mmelelVtxMuonN_.clear();
    mmelelVtxEleP_.clear();
    mmelelVtxEleN_.clear();

    elelVtxVxy_.clear();
    elelVtxVz_.clear();
    elelVtxSigmaVxy_.clear();
    elelVtxReducedChi2_.clear();
    elelVtxDr_.clear();
    elelVtxM_.clear();
    elelVtxPt_.clear();
    elelVtxM2_.clear();
    elelVtxPt2_.clear();
    elelVtxEleP_.clear();
    elelVtxEleN_.clear();

    pcpcVtxVxy_.clear();
    pcpcVtxVz_.clear();
    pcpcVtxSigmaVxy_.clear();
    pcpcVtxReducedChi2_.clear();
    pcpcVtxDr_.clear();
    pcpcVtxM_.clear();
    pcpcVtxPt_.clear();
    pcpcVtxM2_.clear();
    pcpcVtxPt2_.clear();
    pcpcVtxTrackP_.clear();
    pcpcVtxTrackN_.clear();

    mumuVtxVxy_.clear();
    mumuVtxVz_.clear();
    mumuVtxSigmaVxy_.clear();
    mumuVtxReducedChi2_.clear();
    mumuVtxDr_.clear();
    mumuVtxM_.clear();
    mumuVtxPt_.clear();
    mumuVtxM2_.clear();
    mumuVtxPt2_.clear();
    mumuVtxMuonP_.clear();
    mumuVtxMuonN_.clear();

    gsfElsP.clear();
    gsfElsN.clear();
    mmeeTrxP.clear();
    mmeeTrxN.clear();

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

    fired0_ = 0;
    fired1_ = 0;

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
