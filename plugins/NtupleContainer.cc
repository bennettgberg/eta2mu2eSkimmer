#include "NtupleContainer.hh"
#include <iostream>

NtupleContainer::NtupleContainer() : isData_(true) {}

NtupleContainer::~NtupleContainer() {}

void NtupleContainer::SetRecoTree(TTree *tree) { recoT = tree; }
void NtupleContainer::SetGenTree(TTree *tree) { genT = tree; isData_ = false; }

void NtupleContainer::CreateTreeBranches() {

    mmeeTrxP = {};
    mmeeTrxN = {};
    muonsP = {};
    muonsN = {};
    gsfElsP = {};
    gsfElsN = {};
    gsfLowPtElsP = {};
    gsfLowPtElsN = {};

    recoT->Branch("evt", &eventNum_);
    recoT->Branch("lumi_sec", &lumiSec_);
    recoT->Branch("run", &runNum_);
    recoT->Branch("nPV", &npv_);
    recoT->Branch("PV_vx", &pvx_);
    recoT->Branch("PV_vy", &pvy_);
    recoT->Branch("PV_vz", &pvz_);

    //vertices
//    //recoT->Branch("Vertex_mmee_vxy", &mmeeVtxVxy_);
//    recoT->Branch("Vertex_mmee_vx", &mmeeVtxVx_);
//    recoT->Branch("Vertex_mmee_vy", &mmeeVtxVy_);
//    recoT->Branch("Vertex_mmee_vz", &mmeeVtxVz_);
//    recoT->Branch("Vertex_mmee_sigmaVxy", &mmeeVtxSigmaVxy_);
//    recoT->Branch("Vertex_mmee_reduced_chi2", &mmeeVtxReducedChi2_);
//    recoT->Branch("Vertex_mmee_dR", &mmeeVtxDr_);
//    recoT->Branch("Vertex_mmee_muP", &mmeeVtxMuonP_);
//    recoT->Branch("Vertex_mmee_muN", &mmeeVtxMuonN_);
//    recoT->Branch("Vertex_mmee_trackP", &mmeeVtxTrackP_);
//    recoT->Branch("Vertex_mmee_trackN", &mmeeVtxTrackN_);
//    //for debugging only, do NOT include in final ntuple (waste of space)
//    //  actually include only for MC
//    if(!isData_) {
//        recoT->Branch("Vertex_mmee_M", &mmeeVtxM_);
//        recoT->Branch("Vertex_mmee_Pt", &mmeeVtxPt_);
//    }
//    
//    //recoT->Branch("Vertex_mmelel_vxy", &mmelelVtxVxy_);
//    recoT->Branch("Vertex_mmelel_vx", &mmelelVtxVx_);
//    recoT->Branch("Vertex_mmelel_vy", &mmelelVtxVy_);
//    recoT->Branch("Vertex_mmelel_vz", &mmelelVtxVz_);
//    recoT->Branch("Vertex_mmelel_sigmaVxy", &mmelelVtxSigmaVxy_);
//    recoT->Branch("Vertex_mmelel_reduced_chi2", &mmelelVtxReducedChi2_);
//    recoT->Branch("Vertex_mmelel_dR", &mmelelVtxDr_);
//    recoT->Branch("Vertex_mmelel_muP", &mmelelVtxMuonP_);
//    recoT->Branch("Vertex_mmelel_muN", &mmelelVtxMuonN_);
//    recoT->Branch("Vertex_mmelel_eleP", &mmelelVtxEleP_);
//    recoT->Branch("Vertex_mmelel_eleN", &mmelelVtxEleN_);
//    //for debugging only, do NOT include in final ntuple (waste of space)
//    // except for MC
//    if(!isData_) {
//        recoT->Branch("Vertex_mmelel_M", &mmelelVtxM_);
//        recoT->Branch("Vertex_mmelel_Pt", &mmelelVtxPt_);
//    }
//
//    //recoT->Branch("Vertex_elel_vxy", &elelVtxVxy_);
//    recoT->Branch("Vertex_elel_vx", &elelVtxVx_);
//    recoT->Branch("Vertex_elel_vy", &elelVtxVy_);
//    recoT->Branch("Vertex_elel_vz", &elelVtxVz_);
//    recoT->Branch("Vertex_elel_sigmaVxy", &elelVtxSigmaVxy_);
//    recoT->Branch("Vertex_elel_reduced_chi2", &elelVtxReducedChi2_);
//    recoT->Branch("Vertex_elel_dR", &elelVtxDr_);
//    recoT->Branch("Vertex_elel_eleP", &elelVtxEleP_);
//    recoT->Branch("Vertex_elel_eleN", &elelVtxEleN_);
//    //for debugging only, do NOT include in final ntuple (waste of space)
//    // except for MC
//    if(!isData_) {
//        recoT->Branch("Vertex_mmlplp_M", &mmlplpVtxM_);
//        recoT->Branch("Vertex_mmlplp_Pt", &mmlplpVtxPt_);
//    }

    //recoT->Branch("Vertex_mmlplp_vxy", &mmlplpVtxVxy_);
    recoT->Branch("Vertex_mmlplp_vx", &mmlplpVtxVx_);
    recoT->Branch("Vertex_mmlplp_vy", &mmlplpVtxVy_);
    recoT->Branch("Vertex_mmlplp_vz", &mmlplpVtxVz_);
    recoT->Branch("Vertex_mmlplp_sigmaVxy", &mmlplpVtxSigmaVxy_);
    recoT->Branch("Vertex_mmlplp_reduced_chi2", &mmlplpVtxReducedChi2_);
    recoT->Branch("Vertex_mmlplp_dR", &mmlplpVtxDr_);
    recoT->Branch("Vertex_mmlplp_muP", &mmlplpVtxMuonP_);
    recoT->Branch("Vertex_mmlplp_muN", &mmlplpVtxMuonN_);
    recoT->Branch("Vertex_mmlplp_eleP", &mmlplpVtxEleP_);
    recoT->Branch("Vertex_mmlplp_eleN", &mmlplpVtxEleN_);
    if(!isData_) {
        recoT->Branch("Vertex_lplp_M", &lplpVtxM_);
        recoT->Branch("Vertex_lplp_Pt", &lplpVtxPt_);
    }

    //recoT->Branch("Vertex_lplp_vxy", &lplpVtxVxy_);
    recoT->Branch("Vertex_lplp_vx", &lplpVtxVx_);
    recoT->Branch("Vertex_lplp_vy", &lplpVtxVy_);
    recoT->Branch("Vertex_lplp_vz", &lplpVtxVz_);
    recoT->Branch("Vertex_lplp_sigmaVxy", &lplpVtxSigmaVxy_);
    recoT->Branch("Vertex_lplp_reduced_chi2", &lplpVtxReducedChi2_);
    recoT->Branch("Vertex_lplp_dR", &lplpVtxDr_);
    recoT->Branch("Vertex_lplp_eleP", &lplpVtxEleP_);
    recoT->Branch("Vertex_lplp_eleN", &lplpVtxEleN_);
    //for debugging only, do NOT include in final ntuple (waste of space)
    //  except for MC!
    if(!isData_) {
        recoT->Branch("Vertex_lplp_M", &lplpVtxM_);
        recoT->Branch("Vertex_lplp_Pt", &lplpVtxPt_);
    }

//    //recoT->Branch("Vertex_pcpc_vxy", &pcpcVtxVxy_);
//    recoT->Branch("Vertex_pcpc_vx", &pcpcVtxVx_);
//    recoT->Branch("Vertex_pcpc_vy", &pcpcVtxVy_);
//    recoT->Branch("Vertex_pcpc_vz", &pcpcVtxVz_);
//    recoT->Branch("Vertex_pcpc_sigmaVxy", &pcpcVtxSigmaVxy_);
//    recoT->Branch("Vertex_pcpc_reduced_chi2", &pcpcVtxReducedChi2_);
//    recoT->Branch("Vertex_pcpc_dR", &pcpcVtxDr_);
//    recoT->Branch("Vertex_pcpc_trackP", &pcpcVtxTrackP_);
//    recoT->Branch("Vertex_pcpc_trackN", &pcpcVtxTrackN_);
//    //for debugging only, do NOT include in final ntuple (waste of space)
//    // except for MC!
//    if(!isData_) {
//        recoT->Branch("Vertex_pcpc_M", &pcpcVtxM_);
//        recoT->Branch("Vertex_pcpc_Pt", &pcpcVtxPt_);
//    }

    //recoT->Branch("Vertex_mumu_vxy", &mumuVtxVxy_);
    recoT->Branch("Vertex_mumu_vx", &mumuVtxVx_);
    recoT->Branch("Vertex_mumu_vy", &mumuVtxVy_);
    recoT->Branch("Vertex_mumu_vz", &mumuVtxVz_);
    recoT->Branch("Vertex_mumu_sigmaVxy", &mumuVtxSigmaVxy_);
    recoT->Branch("Vertex_mumu_reduced_chi2", &mumuVtxReducedChi2_);
    recoT->Branch("Vertex_mumu_dR", &mumuVtxDr_);
    recoT->Branch("Vertex_mumu_muP", &mumuVtxMuonP_);
    recoT->Branch("Vertex_mumu_muN", &mumuVtxMuonN_);
    //for debugging only, do NOT include in final ntuple (waste of space)
    // except for MC
    if(!isData_) {
        recoT->Branch("Vertex_mumu_M", &mumuVtxM_);
        recoT->Branch("Vertex_mumu_Pt", &mumuVtxPt_);
    }

    //photon conversions
    recoT->Branch("nOnia", &recoNOnia_);
    recoT->Branch("Onia_pt0",  &recoOniaPt0_);
    recoT->Branch("Onia_eta0", &recoOniaEta0_);
    recoT->Branch("Onia_phi0", &recoOniaPhi0_);
    recoT->Branch("Onia_charge0", &recoOniaCharge0_);
    recoT->Branch("Onia_pt1",  &recoOniaPt1_);
    recoT->Branch("Onia_eta1", &recoOniaEta1_);
    recoT->Branch("Onia_phi1", &recoOniaPhi1_);
    recoT->Branch("Onia_charge1", &recoOniaCharge1_);
    //recoT->Branch("Onia_vxy", &recoOniaVtxVxy_);
    recoT->Branch("Onia_vx", &recoOniaVtxVx_);
    recoT->Branch("Onia_vy", &recoOniaVtxVy_);
    recoT->Branch("Onia_vz", &recoOniaVtxVz_);

    //good particles
//    recoT->Branch("nGoodElectron", &recoNGoodElectron_);
//    recoT->Branch("Electron_pt",  &recoElectronPt_);
//    recoT->Branch("Electron_eta", &recoElectronEta_);
//    recoT->Branch("Electron_phi", &recoElectronPhi_);
//    recoT->Branch("Electron_vxy", &recoElectronVxy_);
//    recoT->Branch("Electron_vz",  &recoElectronVz_);
//    recoT->Branch("Electron_charge", &recoElectronCharge_);
//    recoT->Branch("Electron_id", &recoElectronIDResult_);

    recoT->Branch("nGoodLowPtElectron", &recoNGoodLowPtElectron_);
    recoT->Branch("LowPtElectron_pt",  &recoLowPtElectronPt_);
    recoT->Branch("LowPtElectron_eta", &recoLowPtElectronEta_);
    recoT->Branch("LowPtElectron_phi", &recoLowPtElectronPhi_);
    recoT->Branch("LowPtElectron_vxy", &recoLowPtElectronVxy_);
    recoT->Branch("LowPtElectron_vz",  &recoLowPtElectronVz_);
    recoT->Branch("LowPtElectron_charge", &recoLowPtElectronCharge_);
    recoT->Branch("LowPtElectron_id", &recoLowPtElectronIDResult_);

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
//    recoT->Branch("nGoodTrack", &recoNGoodTrk_);
//    recoT->Branch("Track_pt", &recoTrkPt_);
//    recoT->Branch("Track_eta", &recoTrkEta_);
//    recoT->Branch("Track_phi", &recoTrkPhi_);
//    recoT->Branch("Track_charge", &recoTrkCharge_);
//    recoT->Branch("Track_dxy", &recoTrkDxy_);
//    recoT->Branch("Track_dz", &recoTrkDz_);

    //gsfElectron tracks
    //there's exactly one GsfTrack for every electron so nGoodGsfTrack == nGoodElectron
    //recoT->Branch("nGoodGsfTrack", &gsfNGoodTrk_);
//    recoT->Branch("GsfTrack_pt", &gsfTrkPt_);
//    recoT->Branch("GsfTrack_eta", &gsfTrkEta_);
//    recoT->Branch("GsfTrack_phi", &gsfTrkPhi_);
//    recoT->Branch("GsfTrack_charge", &gsfTrkCharge_);
//    recoT->Branch("GsfTrack_dxy", &gsfTrkDxy_);
//    recoT->Branch("GsfTrack_dz", &gsfTrkDz_);

    recoT->Branch("GsfLowPtTrack_pt", &gsfLowPtTrkPt_);
    recoT->Branch("GsfLowPtTrack_eta", &gsfLowPtTrkEta_);
    recoT->Branch("GsfLowPtTrack_phi", &gsfLowPtTrkPhi_);
    recoT->Branch("GsfLowPtTrack_charge", &gsfLowPtTrkCharge_);
    recoT->Branch("GsfLowPtTrack_dxy", &gsfLowPtTrkDxy_);
    recoT->Branch("GsfLowPtTrack_dz", &gsfLowPtTrkDz_);
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

    //mmeeVtxVxy_.clear();
    mmeeVtxVx_.clear();
    mmeeVtxVy_.clear();
    mmeeVtxVz_.clear();
    mmeeVtxSigmaVxy_.clear();
    mmeeVtxReducedChi2_.clear();
    mmeeVtxDr_.clear();
    mmeeVtxM_.clear();
    mmeeVtxPt_.clear();
    mmeeVtxMuonP_.clear();
    mmeeVtxMuonN_.clear();
    mmeeVtxTrackP_.clear();
    mmeeVtxTrackN_.clear();

    //mmelelVtxVxy_.clear();
    mmelelVtxVx_.clear();
    mmelelVtxVy_.clear();
    mmelelVtxVz_.clear();
    mmelelVtxSigmaVxy_.clear();
    mmelelVtxReducedChi2_.clear();
    mmelelVtxDr_.clear();
    mmelelVtxM_.clear();
    mmelelVtxPt_.clear();
    mmelelVtxMuonP_.clear();
    mmelelVtxMuonN_.clear();
    mmelelVtxEleP_.clear();
    mmelelVtxEleN_.clear();

    //mmlplpVtxVxy_.clear();
    mmlplpVtxVx_.clear();
    mmlplpVtxVy_.clear();
    mmlplpVtxVz_.clear();
    mmlplpVtxSigmaVxy_.clear();
    mmlplpVtxReducedChi2_.clear();
    mmlplpVtxDr_.clear();
    mmlplpVtxM_.clear();
    mmlplpVtxPt_.clear();
    mmlplpVtxMuonP_.clear();
    mmlplpVtxMuonN_.clear();
    mmlplpVtxEleP_.clear();
    mmlplpVtxEleN_.clear();

    //lplpVtxVxy_.clear();
    lplpVtxVx_.clear();
    lplpVtxVy_.clear();
    lplpVtxVz_.clear();
    lplpVtxSigmaVxy_.clear();
    lplpVtxReducedChi2_.clear();
    lplpVtxDr_.clear();
    lplpVtxM_.clear();
    lplpVtxPt_.clear();
    lplpVtxEleP_.clear();
    lplpVtxEleN_.clear();

    //elelVtxVxy_.clear();
    elelVtxVx_.clear();
    elelVtxVy_.clear();
    elelVtxVz_.clear();
    elelVtxSigmaVxy_.clear();
    elelVtxReducedChi2_.clear();
    elelVtxDr_.clear();
    elelVtxM_.clear();
    elelVtxPt_.clear();
    elelVtxEleP_.clear();
    elelVtxEleN_.clear();

    //pcpcVtxVxy_.clear();
    pcpcVtxVx_.clear();
    pcpcVtxVy_.clear();
    pcpcVtxVz_.clear();
    pcpcVtxSigmaVxy_.clear();
    pcpcVtxReducedChi2_.clear();
    pcpcVtxDr_.clear();
    pcpcVtxM_.clear();
    pcpcVtxPt_.clear();
    pcpcVtxTrackP_.clear();
    pcpcVtxTrackN_.clear();

    //mumuVtxVxy_.clear();
    mumuVtxVx_.clear();
    mumuVtxVy_.clear();
    mumuVtxVz_.clear();
    mumuVtxSigmaVxy_.clear();
    mumuVtxReducedChi2_.clear();
    mumuVtxDr_.clear();
    mumuVtxM_.clear();
    mumuVtxPt_.clear();
    mumuVtxMuonP_.clear();
    mumuVtxMuonN_.clear();

    gsfElsP.clear();
    gsfElsN.clear();
    gsfLowPtElsP.clear();
    gsfLowPtElsN.clear();
    mmeeTrxP.clear();
    mmeeTrxN.clear();
    muonsP.clear();
    muonsN.clear();

    recoOniaPt0_.clear();
    recoOniaEta0_.clear();
    recoOniaPhi0_.clear();
    recoOniaCharge0_.clear();
    recoOniaPt1_.clear();
    recoOniaEta1_.clear();
    recoOniaPhi1_.clear();
    recoOniaCharge1_.clear();
    //recoOniaVtxVxy_.clear();
    recoOniaVtxVx_.clear();
    recoOniaVtxVy_.clear();
    recoOniaVtxVz_.clear();

    recoElectronPt_.clear();
    recoElectronEta_.clear();
    recoElectronPhi_.clear();
    recoElectronVxy_.clear();
    recoElectronVz_.clear();
    recoElectronCharge_.clear();
    recoElectronIDResult_.clear();
    recoLowPtElectronPt_.clear();
    recoLowPtElectronEta_.clear();
    recoLowPtElectronPhi_.clear();
    recoLowPtElectronVxy_.clear();
    recoLowPtElectronVz_.clear();
    recoLowPtElectronCharge_.clear();
    recoLowPtElectronIDResult_.clear();
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
    
    skipListP.clear();
    skipListN.clear();

    recoNOnia_ = 0;
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
    gsfTrkDxy_.clear();
    gsfTrkDz_.clear();

    gsfLowPtTrkPt_.clear();
    gsfLowPtTrkEta_.clear();
    gsfLowPtTrkPhi_.clear();
    gsfLowPtTrkCharge_.clear();
    gsfLowPtTrkDxy_.clear();
    gsfLowPtTrkDz_.clear();

    // Pile-up and event genweight
    genpuobs_ = -9999;
    genputrue_ = -9999;
    //genwgt_ = -9999;


}
