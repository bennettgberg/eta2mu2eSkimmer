#include <algorithm>
#include <cmath> 
#include <memory>
#include <random>
#include <vector>
#include <boost/format.hpp>
#include <boost/any.hpp>
#include <limits>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "CommonTools/Utils/interface/InvariantMass.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

//new miniAOD includes
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "CommonTools/Egamma/interface/ConversionTools.h" 

#include "TTree.h"
#include "TMath.h"
#include <TLorentzVector.h>

#include "NtupleContainer.hh"
#include "utils.hh"

#include "TVectorD.h"    // for fixing tracks
#include "TMatrixDSym.h" // for fixing tracks

//include for random number generator
#include <cstdlib>

class eta2mu2eAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {

public:
    explicit eta2mu2eAnalyzer(const edm::ParameterSet&);
    ~eta2mu2eAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions&);
    
    typedef std::pair<std::string, float> IdPairf;
    typedef std::pair<std::string, bool> IdPairb;

private:
    bool getCollections(const edm::Event&);

    void beginJob() override;
    void beginRun(edm::Run const&, edm::EventSetup const&) override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endRun(edm::Run const&, edm::EventSetup const&) override;
    void endJob() override;

    TTree *recoT, *genT;
    //for MC only
    TH1F *allGenPt, *matchedGenPt, *alldR, *matchedGenPtE, *alldRE, *alldRMu, *allGenPtMu, *allGenPtEl,
        *matchedGenPtMu, *allGenPtEta, *matchedGenPtEta, *matchedGenPtEtaE; //,
        //*hdRP, *hdRN;
    //max dR between gen and reco particles for a successful gen-match to be declared
    const float drCut = 0.01;
    //max dR betwixt Converted photon (electron) track and electron for a successful match to be declared (and hence the electron is deletted)
    const float drOniaCut = 0.05;
    //min, max invariant mass values for reco'd eta meson for a successful reco to be declared
    const float etaMassMin = 0.0;
    const float etaMassMax = 9999.0;

    NtupleContainer nt;
    edm::Service<TFileService> fs;

    std::mt19937 m_random_generator; 

    bool isData;
    const std::string triggerProcessName_;

    // Tokens
    const edm::EDGetTokenT<pat::MuonCollection> pfRecoMuToken_;
    const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
    const edm::EDGetTokenT<edm::TriggerResults> trigResultsToken_;
    //const edm::EDGetTokenT<trigger::TriggerEvent> trigEventToken_;
    const edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfosToken_;
    const edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken_;
    const edm::EDGetTokenT<reco::VertexCollection>primaryVertexToken_;
    const edm::EDGetTokenT<pat::ElectronCollection> recoElectronToken_;
    //lowpT electrons
    const edm::EDGetTokenT<pat::ElectronCollection> recoLowPtElectronToken_;
    const edm::EDGetTokenT<pat::PhotonCollection> recoPhotonToken_;
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> esToken_;
    const edm::EDGetTokenT<double> rhoToken_;
    const edm::EDGetTokenT<pat::PackedCandidateCollection> trkToken_;
    //photon conversions?
    //const edm::EDGetTokenT<pat::CompositeCandidateCollection> conToken_;
    const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    const edm::EDGetTokenT<reco::ConversionCollection> convsToken_;

    // Handles
    edm::Handle<pat::MuonCollection> recoMuonHandle_;
    edm::Handle<reco::VertexCollection> primaryVertexHandle_;
    edm::Handle<reco::GenParticleCollection> genParticleHandle_;
    edm::Handle<edm::TriggerResults> trigResultsHandle_;
    //edm::Handle<trigger::TriggerEvent> trigEventHandle_;
    edm::Handle<std::vector<PileupSummaryInfo>> pileupInfosHandle_;
    edm::Handle<GenEventInfoProduct> genEvtInfoHandle_;
    edm::Handle<pat::ElectronCollection> recoElectronHandle_;
    //low pT
    edm::Handle<pat::ElectronCollection> recoLowPtElectronHandle_;
    edm::Handle<pat::PhotonCollection> recoPhotonHandle_;
    edm::Handle<double> rhoHandle_;
    edm::Handle<pat::PackedCandidateCollection> trkHandle_;
    //for photon conversions
    //edm::Handle<pat::CompositeCandidateCollection> conHandle_;
    edm::Handle<reco::BeamSpot> bsHandle_;
    edm::Handle<reco::ConversionCollection> convsHandle_;
    
    std::vector<std::string> triggerPathsWithoutVersionNum_;
    std::vector<std::string> triggerPathsWithVersionNum_;
    std::vector<bool> trigExist_;
    HLTConfigProvider hltConfig_;

};


eta2mu2eAnalyzer::eta2mu2eAnalyzer(const edm::ParameterSet& ps):
    isData(ps.getParameter<bool>("isData")),
    triggerProcessName_(ps.getParameter<std::string>("triggerProcessName")),
    
    pfRecoMuToken_(consumes<pat::MuonCollection>(ps.getParameter<edm::InputTag>("muon_collection"))),
    genParticleToken_(consumes<reco::GenParticleCollection>(ps.getParameter<edm::InputTag>("genParticle"))),
    trigResultsToken_(consumes<edm::TriggerResults>(ps.getParameter<edm::InputTag>("trigResult"))),
    //trigEventToken_(consumes<trigger::TriggerEvent>(ps.getParameter<edm::InputTag>("trigEvent"))),
    pileupInfosToken_(consumes<std::vector<PileupSummaryInfo>>(ps.getParameter<edm::InputTag>("pileups"))),
    genEvtInfoToken_(consumes<GenEventInfoProduct>(ps.getParameter<edm::InputTag>("genEvt"))),
    primaryVertexToken_(consumes<reco::VertexCollection>(ps.getParameter<edm::InputTag>("primary_vertices"))),
    recoElectronToken_(consumes<pat::ElectronCollection>(ps.getParameter<edm::InputTag>("electron_collection"))),
    recoLowPtElectronToken_(consumes<pat::ElectronCollection>(ps.getParameter<edm::InputTag>("lowpt_electron_collection"))),
    recoPhotonToken_(consumes<pat::PhotonCollection>(ps.getParameter<edm::InputTag>("photon_collection"))),
    esToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
    rhoToken_(consumes<double>(ps.getParameter<edm::InputTag>("rho"))),
    trkToken_(consumes<pat::PackedCandidateCollection>(ps.getParameter<edm::InputTag>("packed_candidate"))),
    //conToken_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("composite_candidate")))
    bsToken_(consumes<reco::BeamSpot>(ps.getParameter<edm::InputTag>("beamspot"))),
    convsToken_(consumes<reco::ConversionCollection>(ps.getParameter<edm::InputTag>("conversions")))
{
    usesResource("TFileService");
    m_random_generator = std::mt19937(37428479);
}

eta2mu2eAnalyzer::~eta2mu2eAnalyzer() = default;

void eta2mu2eAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    // Only specify tags with reasonable defaults -- check out cfg for others

    edm::ParameterSetDescription desc;
    desc.add<bool>("isData", 0);
    desc.add<std::string>("triggerProcessName", "HLT");

    desc.add<edm::InputTag>("muon_collection", edm::InputTag("slimmedMuons"));
    desc.add<edm::InputTag>("electron_collection", edm::InputTag("slimmedElectrons"));
    desc.add<edm::InputTag>("lowpt_electron_collection", edm::InputTag("slimmedLowPtElectrons"));
    desc.add<edm::InputTag>("photon_collection", edm::InputTag("slimmedPhotons"));
    desc.add<edm::InputTag>("primary_vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
    desc.add<edm::InputTag>("genParticle", edm::InputTag("prunedGenParticles"));
    desc.add<edm::InputTag>("trigResult", edm::InputTag("TriggerResults", "", "HLT"));
    //desc.add<edm::InputTag>("trigEvent", edm::InputTag("hltTriggerSummaryAOD", "", "HLT"));
    desc.add<edm::InputTag>("pileups", edm::InputTag("slimmedAddPileupInfo"));
    desc.add<edm::InputTag>("genEvt", edm::InputTag("generator"));
    desc.add<edm::InputTag>("rho", edm::InputTag("fixedGridRhoFastjetAll"));
    desc.add<edm::InputTag>("packed_candidate", edm::InputTag("packedPFCandidates"));
    //desc.add<edm::InputTag>("composite_candidate", edm::InputTag("oniaPhotonCandidates", "conversions"));
    desc.add<edm::InputTag>("beamspot", edm::InputTag("offlineBeamSpot"));
    desc.add<edm::InputTag>("conversions", edm::InputTag("reducedEgamma", "reducedConversions"));
    
    descriptions.add("eta2mu2eAnalyzer", desc);
}

void eta2mu2eAnalyzer::beginJob()
{
    recoT = fs->make<TTree>("recoT", "recoT");
    nt.SetRecoTree(recoT);
    if (!isData) {
        genT = fs->make<TTree>("genT", "genT");
        nt.SetGenTree(genT);
        allGenPt = new TH1F("allGenPt", "allGenPt", 10000, 0., 100.);
        matchedGenPt = new TH1F("matchedGenPt", "matchedGenPt", 10000, 0., 100.);
        matchedGenPtE = new TH1F("matchedGenPtE", "matchedGenPtE", 10000, 0., 100.);
        alldR = new TH1F("alldR", "alldR", 10000, 0., 1.);
        alldRE = new TH1F("alldRE", "alldRE", 10000, 0., 1.);
        alldRMu = new TH1F("alldRMu", "alldRMu", 10000, 0., 1.);
        allGenPtMu = new TH1F("allGenPtMu", "allGenPtMu", 10000, 0., 100.);
        allGenPtEl = new TH1F("allGenPtEl", "allGenPtEl", 10000, 0., 100.);
        matchedGenPtMu = new TH1F("matchedGenPtMu", "matchedGenPtMu", 10000, 0., 100.);
        allGenPtEta = new TH1F("allGenPtEta", "allGenPtEta", 10000, 0., 100.);
        matchedGenPtEta = new TH1F("matchedGenPtEta", "matchedGenPtEta", 10000, 0., 100.);
        matchedGenPtEtaE = new TH1F("matchedGenPtEtaE", "matchedGenPtEtaE", 10000, 0., 100.);

        srand(time(NULL));

        // Generate a random number.
        int randomInt = rand();

        // Print the random integer
        std::cout << "Random Integer: " << randomInt << std::endl;
        nt.MClumiblock_ = randomInt;
    }

    //hdRP = new TH1F("hdRP", "hdRP", 10000, 0.0, 10.0);
    //hdRN = new TH1F("hdRN", "hdRN", 10000, 0.0, 10.0);
    nt.CreateTreeBranches();
}


void eta2mu2eAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    using namespace edm;

    bool changed = true;
    if (hltConfig_.init(iRun, iSetup, triggerProcessName_, changed)) {
        if (changed) {
            LogInfo("HLTConfig") << "eta2mu2eAnalyzer::beginRun: " << "hltConfig init for Run" << iRun.run();
            hltConfig_.dump("ProcessName");
            hltConfig_.dump("GlobalTag");
            hltConfig_.dump("TableName");
        }
    } 
    else {
        LogError("HLTConfig") << "eta2mu2eAnalyzer::beginRun: config extraction failure with triggerProcessName -> " << triggerProcessName_;
        return;
    }

    // Add trigger paths if they exist
    triggerPathsWithoutVersionNum_.clear();
    triggerPathsWithVersionNum_.clear();
    trigExist_.clear();

    vector<std::string> triggerNames { "HLT_Dimuon0_Jpsi3p5_Muon2",
        "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R",
        "HLT_Dimuon0_Jpsi_L1_NoOS",
        "HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R",
        "HLT_Dimuon0_Jpsi_NoVertexing_NoOS",
        "HLT_Dimuon0_Jpsi_NoVertexing",
        "HLT_Dimuon0_Jpsi",
        "HLT_Dimuon0_LowMass_L1_0er1p5R",
        "HLT_Dimuon0_LowMass_L1_0er1p5",
        "HLT_Dimuon0_LowMass_L1_4R",
        "HLT_Dimuon0_LowMass_L1_4",
        "HLT_Dimuon0_LowMass_L1_TM530",
        "HLT_Dimuon0_LowMass",
        "HLT_Dimuon0_Upsilon_L1_4p5NoOS",
        "HLT_Dimuon0_Upsilon_L1_4p5",
        "HLT_Dimuon0_Upsilon_L1_4p5er2p0M",
        "HLT_Dimuon0_Upsilon_L1_4p5er2p0",
        "HLT_Dimuon0_Upsilon_L1_5M",
        "HLT_Dimuon0_Upsilon_L1_5",
        "HLT_Dimuon0_Upsilon_Muon_L1_TM0",
        "HLT_Dimuon0_Upsilon_Muon_NoL1Mass",
        "HLT_Dimuon0_Upsilon_NoVertexing",
        "HLT_Dimuon10_PsiPrime_Barrel_Seagulls",
        "HLT_Dimuon10_Upsilon_y1p4",
        "HLT_Dimuon12_Upsilon_y1p4",
        "HLT_Dimuon14_Phi_Barrel_Seagulls",
        "HLT_Dimuon14_PsiPrime_noCorrL1",
        "HLT_Dimuon14_PsiPrime",
        "HLT_Dimuon18_PsiPrime_noCorrL1",
        "HLT_Dimuon18_PsiPrime",
        "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
        "HLT_Dimuon24_Phi_noCorrL1",
        "HLT_Dimuon24_Upsilon_noCorrL1",
        "HLT_Dimuon25_Jpsi_noCorrL1",
        "HLT_Dimuon25_Jpsi",
        "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05",
        "HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon",
        "HLT_DoubleMu3_TkMu_DsTau3Mu",
        "HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass",
        "HLT_DoubleMu3_Trk_Tau3mu",
        "HLT_DoubleMu4_3_Bs",
        "HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG",
        "HLT_DoubleMu4_3_Jpsi",
        "HLT_DoubleMu4_3_LowMass",
        "HLT_DoubleMu4_3_Photon4_BsToMMG",
        "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
        "HLT_DoubleMu4_JpsiTrk_Bc",
        "HLT_DoubleMu4_Jpsi_Displaced",
        "HLT_DoubleMu4_Jpsi_NoVertexing",
        "HLT_DoubleMu4_LowMass_Displaced",
        "HLT_DoubleMu4_MuMuTrk_Displaced",
        "HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL",
        "HLT_Mu20_TkMu0_Phi",
        "HLT_Mu25_TkMu0_Onia",
        "HLT_Mu25_TkMu0_Phi",
        "HLT_Mu30_TkMu0_Psi",
        "HLT_Mu30_TkMu0_Upsilon",
        "HLT_Mu4_L1DoubleMu",
        "HLT_Mu7p5_L2Mu2_Jpsi",
        "HLT_Mu7p5_L2Mu2_Upsilon",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15",
        "HLT_Trimuon5_3p5_2_Upsilon_Muon",
        "HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon" };

    for(std::string tName : triggerNames) {
        triggerPathsWithoutVersionNum_.emplace_back(tName);
    }
    
    const std::vector<std::string>& pathNames = hltConfig_.triggerNames();
    for (auto trigPathNoVersion : triggerPathsWithoutVersionNum_) {
        auto matchedPaths(hltConfig_.restoreVersion(pathNames, trigPathNoVersion));
        if (matchedPaths.size() == 0) {
            LogWarning("TriggerNotFound") << "Could not find matched full trigger path with --> " << trigPathNoVersion;
            triggerPathsWithVersionNum_.push_back("None");
            trigExist_.push_back(false);
        }
        else {
            trigExist_.push_back(true);
            triggerPathsWithVersionNum_.push_back(matchedPaths[0]);
            if (hltConfig_.triggerIndex(matchedPaths[0]) >= hltConfig_.size()) {
                LogError("TriggerError") << "Cannot find trigger path --> " << matchedPaths[0];
                return;
            }
            //else {
            //    std::cout << "Seems to be no problem with trigger path " << trigPathNoVersion << " : " << matchedPaths[0] << std::endl;
            //}
        }
    }
}


bool eta2mu2eAnalyzer::getCollections(const edm::Event& iEvent) {
    using namespace edm;

    char error_msg[] = "eta2mu2eAnalyzer::GetCollections: Error in getting product %s from Event!";

    bool ret = true;
    auto getHandle = [&]<typename T>(const EDGetTokenT<T> &token, Handle<T> &handle, std::string name) {
        iEvent.getByToken(token, handle);
        if (!handle.isValid()) {
            LogError("HandleError") << boost::str(boost::format(error_msg) % name);
            ret = false;
        }
    };

    getHandle(pfRecoMuToken_, recoMuonHandle_, "muon_collection");
    getHandle(primaryVertexToken_, primaryVertexHandle_, "primary_vertices");
    getHandle(trigResultsToken_, trigResultsHandle_, "trigResults");
    //getHandle(trigEventToken_, trigEventHandle_, "trigEvent");
    getHandle(recoElectronToken_, recoElectronHandle_, "electron_collection");
    getHandle(recoLowPtElectronToken_, recoLowPtElectronHandle_, "lowpt_electron_collection");
    getHandle(recoPhotonToken_, recoPhotonHandle_, "photon_collection");
    getHandle(rhoToken_, rhoHandle_, "rho");
    getHandle(trkToken_, trkHandle_, "packed_candidate");
    //getHandle(conToken_, conHandle_, "composite_candidate");
    getHandle(bsToken_, bsHandle_, "beamspot");
    getHandle(convsToken_, convsHandle_, "conversions");
    if (!isData) {
        getHandle(genEvtInfoToken_, genEvtInfoHandle_, "genEventInfo");
        getHandle(genParticleToken_, genParticleHandle_, "genParticle");
        getHandle(pileupInfosToken_, pileupInfosHandle_, "pileupInfos");
    }
    
    return ret;
}

void eta2mu2eAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using std::cout, std::vector, std::endl;

    if (!getCollections(iEvent))
        return;

    // Clear branches before filling
    nt.ClearTreeBranches();

    // Start filling

    // Event information
    nt.eventNum_ = iEvent.id().event();
    if(!isData) {
        nt.lumiSec_ = nt.MClumiblock_;
    }
    else {
        nt.lumiSec_ = iEvent.luminosityBlock();
    }
    nt.runNum_ = iEvent.id().run();
    if(nt.runNum_ == 357271) {
        std::cout << "Run " << (int)nt.runNum_ << ", LS " << (int)nt.lumiSec_ << ", Event " << (int)nt.eventNum_  << std::endl;
    }
    //nt.npv_ = *primaryVertexFilterHandle_;
    nt.npv_ = primaryVertexHandle_->size();

    // Primary vertex
    reco::Vertex pv = *primaryVertexHandle_->begin();
    nt.pvx_ = pv.x();
    nt.pvy_ = pv.y();
    nt.pvz_ = pv.z();

//    // Add all electrons to ntuple, regardless of ID
    // don't care about NElectron, only care about NGoodElectron
    //nt.recoNElectron_ = recoElectronHandle_->size();
    nt.recoNGoodElectron_ = 0;
    //get the tracks from the electron objects
    vector<reco::GsfTrackRef> elTracksP{};    
    vector<reco::GsfTrackRef> elTracksN{};    

    //nt.gsfNTrk_ = recoElectronHandle_->size();
    //nt.gsfNGoodTrk_ = 0;
//    std::cout << "**Regular Electrons**" << std::endl;
    for (size_t i = 0; i < recoElectronHandle_->size(); i++) {
        pat::ElectronRef electronRef(recoElectronHandle_, i);
        nt.recoElectronPt_.push_back(electronRef->pt());
        nt.recoElectronEta_.push_back(electronRef->eta());
        nt.recoElectronPhi_.push_back(electronRef->phi());
        nt.recoElectronVxy_.push_back(electronRef->trackPositionAtVtx().rho()); //?????
        nt.recoElectronVz_.push_back(electronRef->trackPositionAtVtx().z());
        nt.recoElectronCharge_.push_back(electronRef->charge());
        //std::cout << "Electron " << (int)nt.recoNGoodElectron_ << ": pT=" << electronRef->pt() << ", eta=" << electronRef->eta() << ", phi=" << electronRef->phi() 
        //    << ", charge=" << (int)electronRef->charge() << std::endl;

        //try getting the electron ID
        //if ( i == 0 ) {
        //    vector<IdPairf> IDsList = electronRef->electronIDs();
        //    std::cout << "List of available electron IDs: " << std::endl;
        //    for( IdPairf x : IDsList ) {
        //        std::cout << x.first << " : " << x.second << std::endl;
        //    }
        //}

        //Directly from CMSSW code:
        //// ---- methods for electron ID ----
        /// Returns a specific electron ID associated to the pat::Electron given its name
        // For cut-based IDs, the value map has the following meaning:
        // 0: fails,
        // 1: passes electron ID only,
        // 2: passes electron Isolation only,
        // 3: passes electron ID and Isolation only,
        // 4: passes conversion rejection,
        // 5: passes conversion rejection and ID,
        // 6: passes conversion rejection and Isolation,
        // 7: passes the whole selection.
        // For more details have a look at:
        // https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID
        // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCategoryBasedElectronID
        // Note: an exception is thrown if the specified ID is not available

        //format of my Electron ID (8 bits):
        // 0: wp80 ID, 1: wp80 conversion rejection, 2: wp90 ID, 3: wp90 conversion rejection, 4: loose ID, 5: loose conv. rej. (6,7: 0)
        uint8_t full_id = 0;
        int idnum = 0;
        //loop over the 3 different elec IDs (wp80, wp90, loose)
        for(std::string id : nt.electronIDs) {
            uint8_t thisid = electronRef->electronID(id);
            //std::cout << id << " : " << (int)thisid << "; ";
            uint8_t add_id = 0;
            if ( thisid % 2 == 1 ) {
                add_id += 1;
            }
            if ( thisid > 3 ) {
                add_id += 2;
            }
            full_id = full_id + (add_id << (2*idnum));
            idnum++;
        }
        //std::cout << "full_id : " << (int)full_id << std::endl;
        nt.recoElectronIDResult_.push_back( full_id );

        reco::GsfTrackRef eltrack = electronRef->gsfTrack();
        constexpr auto missingHitType = reco::HitPattern::MISSING_INNER_HITS;
        uint8_t mHits = eltrack->hitPattern().numberOfLostHits(missingHitType);
        //std::cout << "mHits: " << (int)mHits << std::endl;
        nt.recoElectronNMhits_.push_back(mHits);

        //const reco::GsfElectron *gsfElectron = static_cast<const reco::GsfElectron*>(*electronRef);

        //if (gsfElectron) {
        //    // Successfully casted to reco::GsfElectron
        //    // Now you can use the gsfElectron object
        //} else {
        //    // Failed to cast, the pat::Electron is not a reco::GsfElectron
        //    std::cout << "Error: not a good gsfElectron :/" << std::endl;
        //}

        ////see if there's a matched conversion candidate?
        bool conversion = ConversionTools::hasMatchedConversion( *electronRef, *convsHandle_, bsHandle_->position());
        nt.recoElectronConvVeto_.push_back(conversion);

        nt.gsfTrkPt_.push_back(eltrack->pt());
        nt.gsfTrkEta_.push_back(eltrack->eta());
        nt.gsfTrkPhi_.push_back(eltrack->phi());
        nt.gsfTrkCharge_.push_back(eltrack->charge());

        nt.gsfTrkDxy_.push_back(eltrack->dxy());
        nt.gsfTrkDz_.push_back(eltrack->dz());

        if ( electronRef->charge() > 0 ) {
            elTracksP.push_back(eltrack);
            nt.gsfElsP.push_back(nt.recoNGoodElectron_);
        }
        else {
            elTracksN.push_back(eltrack);
            nt.gsfElsN.push_back(nt.recoNGoodElectron_);
        }
        //cout << "about to increment; just pushed back: " << (int)nt.recoNGoodElectron_ << "/" << (int)(recoElectronHandle_->size()) << " good electrons." << endl;
        nt.recoNGoodElectron_++;
        //cout<<"just incremented NGoodElectron: " << (int)nt.recoNGoodElectron_ << "; i: " << (int)i << "; handle size: " << (int)(recoElectronHandle_->size()) << endl;
    }

    nt.recoNGoodLowPtElectron_ = 0;
    //get the tracks from the electron objects
    vector<reco::GsfTrackRef> elLowPtTracksP{};    
    vector<reco::GsfTrackRef> elLowPtTracksN{};    
    //now add also the Low pT electrons!!
    //std::cout << "**Low pT Electrons**" << std::endl;
    for (size_t i = 0; i < recoLowPtElectronHandle_->size(); i++) {
        pat::ElectronRef electronRef(recoLowPtElectronHandle_, i);
        nt.recoLowPtElectronPt_.push_back(electronRef->pt());
        nt.recoLowPtElectronEta_.push_back(electronRef->eta());
        nt.recoLowPtElectronPhi_.push_back(electronRef->phi());
        nt.recoLowPtElectronVxy_.push_back(electronRef->trackPositionAtVtx().rho()); //?????
        nt.recoLowPtElectronVz_.push_back(electronRef->trackPositionAtVtx().z());
        nt.recoLowPtElectronCharge_.push_back(electronRef->charge());

        //std::cout << "LowPtElectron " << (int)nt.recoNGoodElectron_ <<": pT="<<electronRef->pt()<<", eta="<<electronRef->eta()<< ", phi=" <<electronRef->phi() 
        //    << ", charge=" << (int)electronRef->charge() << std::endl;

        //format of my Electron ID (8 bits):
        // 0: wp80 ID, 1: wp80 conversion rejection, 2: wp90 ID, 3: wp90 conversion rejection, 4: loose ID, 5: loose conv. rej. (6,7: 0)
        uint8_t full_id = 0;
        int idnum = 0;
        //loop over the 3 different elec IDs (wp80, wp90, loose)
        //for(std::string id : nt.electronIDs) {
        for(std::string id : nt.lowPtElectronIDs) {
            uint8_t thisid = electronRef->electronID(id);
            uint8_t add_id = 0;
            if ( thisid % 2 == 1 ) {
                add_id += 1;
            }
            if ( thisid > 3 ) {
                add_id += 2;
            }
            full_id = full_id + (add_id << (2*idnum));
            idnum++;
        }
        nt.recoLowPtElectronIDResult_.push_back( full_id );
        reco::GsfTrackRef eltrack = electronRef->gsfTrack();

        nt.gsfLowPtTrkPt_.push_back(eltrack->pt());
        nt.gsfLowPtTrkEta_.push_back(eltrack->eta());
        nt.gsfLowPtTrkPhi_.push_back(eltrack->phi());
        nt.gsfLowPtTrkCharge_.push_back(eltrack->charge());
        nt.gsfLowPtTrkDxy_.push_back(eltrack->dxy());
        nt.gsfLowPtTrkDz_.push_back(eltrack->dz());


        if ( electronRef->charge() > 0 ) {
            elLowPtTracksP.push_back(eltrack);
            nt.gsfLowPtElsP.push_back(nt.recoNGoodLowPtElectron_);
        }
        else {
            elLowPtTracksN.push_back(eltrack);
            nt.gsfLowPtElsN.push_back(nt.recoNGoodLowPtElectron_);
        }
        nt.recoNGoodLowPtElectron_++;
    }

    //cout << "Event " << (int)nt.eventNum_ << ": " << (int)nt.recoNGoodElectron_ << " good electrons total. " << (int)elTracksP.size() << " positive and " << (int)elTracksN.size() << " negative." << std::endl;
    // Also add all photons to ntuple, regardless of ID
    // Photon ID only produces 1 or 0
    //nt.recoNPhoton_ = recoPhotonHandle_->size();
    nt.recoNGoodPhoton_ = 0;
    for (size_t i = 0; i < recoPhotonHandle_->size(); i++) {
        pat::PhotonRef photonRef(recoPhotonHandle_, i);
        nt.recoPhotonPt_.push_back(photonRef->pt());
        nt.recoPhotonEta_.push_back(photonRef->eta());
        nt.recoPhotonPhi_.push_back(photonRef->phi());
        //for(std::string id : nt.photonIDs) {
        //    nt.recoPhotonIDResult_[id].push_back( (float) photonRef->photonID(id) );
        //}
        nt.recoNGoodPhoton_++;
        //my photon ID format:
        // bit 0: wp80
        // bit 1: wp90
        uint8_t full_id = 0;
        int idnum = 0;
        for(std::string id : nt.photonIDs) {
            uint8_t thisid = photonRef->photonID(id);
            //std::cout << id << " : " << (int)thisid << "; ";
            full_id += (thisid << idnum);
            idnum++;
        } 
        //std::cout << "full_id : " << (int)full_id << std::endl;
        nt.recoPhotonIDResult_.push_back( full_id );
    }
    
    //vector of all Muons (Positive and Negative separately)
    vector<pat::Muon> muonsP {};
    vector<pat::Muon> muonsN {};
    nt.recoNGoodMuon_ = 0;
    for (size_t i = 0; i < recoMuonHandle_->size(); i++) {
        pat::MuonRef muonRef(recoMuonHandle_, i);
        //muon info will be added later, once we're sure this is a useful muon
        //nt.recoMuonPt_.push_back(muonRef->pt());
        //nt.recoMuonEta_.push_back(muonRef->eta());
        //nt.recoMuonPhi_.push_back(muonRef->phi());
        //nt.recoMuonCharge_.push_back(muonRef->charge());
        if ( muonRef->charge() > 0 ) {
            muonsP.push_back(* muonRef );
        }
        else {
            muonsN.push_back(* muonRef );
        }
        ////nt.recoMuonIDResult_.push_back( (float) (muonRef->muonID("All")) );
        //NGoodMuon is set later
        ////nt.recoNGoodMuon_++;
    }

    // Pick pair of muons with smallest vertex chi square fit for all collection combos
    edm::ESHandle<TransientTrackBuilder> theB = iSetup.getHandle(esToken_);
    //deprecated!
    //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
    KalmanVertexFitter kvf(true);

    //all tracks (will be mostly pions methinks) -- separated betwixt positive and negative charge
    vector<reco::Track> allTracksP{};    //positive only
    vector<reco::Track> allTracksN{};    //negative only

    //now get the PackedCandidate tracks
    //nt.recoNTrk_ = trkHandle_->size();
    nt.recoNGoodTrk_ = 0;
    for (pat::PackedCandidateCollection::const_iterator iTra1 = trkHandle_->begin(); iTra1 != trkHandle_->end(); iTra1++) {
        //try setting track mass to electron mass to see what happens
        //iTra1->setMass(ele_mass); //doesn't work
        if (!(iTra1->hasTrackDetails())) continue;

        //electrons only
        //float mass = iTra1->mass();
        //if( fabs(mass - ele_mass) > (ele_mass/10.0) ) continue;
        if( iTra1->charge() == 0 ) continue;
        //std::cout << "mass: " << mass << std::endl;

        auto iTrack1 = iTra1->bestTrack();
        if (iTrack1 == nullptr) continue;
        //
        if (iTrack1->pt()<0.5) continue;
        //if (iTrack1->charge() < 0.5) continue; // positive
        if (!(iTrack1->quality(reco::TrackBase::highPurity))) continue;
        //
        if( fabs(iTrack1->eta()) > 2.9 ) continue;
        //
        //TLorentzVector p4el1;
        //p4el1.SetPtEtaPhiM(iTrack1->pt(),iTrack1->eta(),iTrack1->phi(), ele_mass);
        //
        reco::TransientTrack elec1TT = theB->build( fix_track(iTra1->bestTrack()) );
        //reco::TransientTrack elec1TT = theB->build( fix_track(iTra1->bestTrack()) );
        if (!elec1TT.isValid()) continue;
 
        if ( iTrack1->charge() > 0.5) {
            allTracksP.push_back(*iTrack1);
        }
        else if ( iTrack1->charge() < -0.5 ) {
            allTracksN.push_back(*iTrack1); 
        }
        else {
            continue;
        }
        //nt.recoNGoodTrk_++;
        //nt.recoTrkPt_.push_back(iTrack1->pt());
        //nt.recoTrkEta_.push_back(iTrack1->eta());
        //nt.recoTrkPhi_.push_back(iTrack1->phi());
        //nt.recoTrkCharge_.push_back(iTrack1->charge());
    }


    // Assign each trigger result to a different bit
    nt.fired0_ = 0;
    nt.fired1_ = 0;
    for (size_t i = 0; i < triggerPathsWithVersionNum_.size(); i++) {
        if (trigExist_.at(i)) {
            std::string trigPath = triggerPathsWithVersionNum_[i];
            //first 64 triggers belong to the first trigger word, next few to the second one.
            if(i < 64) {
                //bool firstfired = nt.fired0_&(1<<11);
                nt.fired0_ |= (trigResultsHandle_->accept(hltConfig_.triggerIndex(trigPath)) << i);
                //std::cout << "Trigger bit " << (int)i << ": " << trigPath << ". New fired0: " << (int)nt.fired0_ << std::endl;
                //if(!firstfired && (nt.fired0_&(1<<11))) std::cout << "trig bit 11 found!!!: " << trigPath << std::endl;
            }
            else {
                nt.fired1_ |= (trigResultsHandle_->accept(hltConfig_.triggerIndex(trigPath)) << (i-64));
            }
        }
        else {
            if(i < 64) {
                nt.fired0_ |= (0 <<i);
            }
            else {
                nt.fired1_ |= (0 <<(i-64));
            }
        }
    }
    
    //std::cout << "Number of onia ConvertedPhoton Candidates: " << (int)conHandle_->size() << std::endl;
    //std::cout << "Printing onia photons:" << std::endl;
    ////get converted photons, see how many there are, if their electron tracks can be gen-matched to regular electrons?
//    int ctr = 0;
//    for (pat::CompositeCandidateCollection::const_iterator iCon1 = conHandle_->begin(); iCon1 != conHandle_->end(); iCon1++) {
//        //const reco::Track* trk0 = iCon1->bestTrack(); //track0;
//        //if(!trk0) {
//        //    std::cout << "ConPhot " << ctr << ": bestTrack is null!!" << std::endl;
//        //    continue;
//        //}
//        bool hasTrack0 = iCon1->hasUserData("track0");
//        if(!hasTrack0) {
//            std::cout << "No track0 :(" << std::endl;
//            continue;
//        }
//        //reco::Vertex conVert = iCon1->vertex();
//        const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> conVert = iCon1->vertex();
//        float vx = conVert.X();
//        float vy = conVert.Y();
//        //float vxy = (std::sqrt(conVert.X()*conVert.X() + conVert.Y()*conVert.Y()));
//        float vz = conVert.Z();
//        //nt.recoOniaVtxVxy_.push_back(vxy);
//        nt.recoOniaVtxVx_.push_back(vx);
//        nt.recoOniaVtxVy_.push_back(vy);
//        nt.recoOniaVtxVz_.push_back(vz);
//        //else {
//        //    std::cout << "Has track0!!" << std::endl;
//        //}
//        const reco::Track* trk0 = iCon1->userData<reco::Track>("track0");
//        const reco::Track* trk1 = iCon1->userData<reco::Track>("track1");
//        //std::cout << "ConPhot " << ctr << " trk0: pT=" << trk0->pt() << ", eta=" << trk0->eta() << ", phi=" << trk0->phi() 
//        //   << ", charge=" << trk0->charge() << std::endl;
//        //std::cout << "ConPhot " << ctr << " trk1: pT=" << trk1->pt() << ", eta=" << trk1->eta() << ", phi=" << trk1->phi() 
//        //   << ", charge=" << trk1->charge() << std::endl;
//
//        nt.recoOniaPt0_.push_back(trk0->pt());
//        nt.recoOniaEta0_.push_back(trk0->eta());
//        nt.recoOniaPhi0_.push_back(trk0->phi());
//        nt.recoOniaCharge0_.push_back((uint8_t)trk0->charge());
//        nt.recoOniaPt1_.push_back(trk1->pt());
//        nt.recoOniaEta1_.push_back(trk1->eta());
//        nt.recoOniaPhi1_.push_back(trk1->phi());
//        nt.recoOniaCharge1_.push_back((uint8_t)trk1->charge());
//        nt.recoNOnia_++;
//
//        const reco::Track* trkP = (trk0->charge()>0 ? trk0 : trk1);
//        const reco::Track* trkN = (trk1->charge()<0 ? trk1 : trk0);
//        if(trkP == trkN) {
//            std::cout << "Error! trk0 and trk1 have same charge!! Cand " << ctr << "; Qtrk0: " << trk0->charge() << ", Qtrk1: " << trk1->charge() << std::endl;
//        }
//        //see if can genmatch with one of the electrons
//        float mindRP = 9999.0;
//        float mindRN = 9999.0;
//        //which index has the lowest dR betwixt the pos and neg track respectively?
//        int argp = -1;
//        int argn = -1;
//        //for(reco::GsfTrackRef tp : elTracksP) {
//        int ii = 0;
//        for(reco::GsfTrackRef tp : elLowPtTracksP) {
//            float dR = reco::deltaR(*trkP, *tp);
//            //std::cout << "dR: " << dR << std::endl;
//            if(dR < mindRP) {
//                mindRP = dR;
//                argp = ii;
//            }
//            ii++;
//        }
//        ii = 0;
//        //for(reco::GsfTrackRef tn : elTracksN) {
//        for(reco::GsfTrackRef tn : elLowPtTracksN) {
//            float dR = reco::deltaR(*trkN, *tn);
//            //std::cout << "dR: " << dR << std::endl;
//            if(dR < mindRN) {
//                mindRN = dR;
//                argn = ii;
//            }
//        }
//        //if(mindRP<9999 && mindRN<9999) {
//        //    //std::cout << "mindR for Conv. phot. tracks vs. GsfElectron tracks, pos/neg: " << mindRP << "/" << mindRN << std::endl;
//        //    hdRP->Fill(mindRP);
//        //    hdRN->Fill(mindRN);
//        //}
//        //now if the dR values are small enough, remove the matched electron/positron from consideration
//        if(nt.removeOnia && mindRP<drOniaCut){
//            nt.skipListP.push_back(argp);
//            //std::cout << "SKIP elecP " << argp << " out of " << (int)elLowPtTracksP.size() << ": dR=" << mindRP << std::endl;
//        }
//        if(nt.removeOnia && mindRN<drOniaCut){
//            nt.skipListN.push_back(argn);
//            //std::cout << "SKIP elecN " << argn << " out of " << (int)elLowPtTracksN.size() << ": dR=" << mindRN << std::endl;
//        }
//        ctr++;
//    } 

    /****** GEN INFO *******/

    //std::cout << "done computing vertices." << std::endl;
    if (!isData) {

        nt.nGen_ = (int)genParticleHandle_->size();
        
        //TODO: add gen weight info 
        // Gen weight
        //nt.genwgt_ = genEvtInfoHandle_->weight();

        // Pile-up
        for (const auto & pileupInfo : *pileupInfosHandle_) {
            if (pileupInfo.getBunchCrossing() == 0) {
                nt.genpuobs_ = pileupInfo.getPU_NumInteractions();
                nt.genputrue_ = pileupInfo.getTrueNumInteractions();
                break;
            }
        }

        //gen level eta meson
        TLorentzVector genEtaVec;
        //reco'd eta meson using PC tracks
        TLorentzVector recoEtaVec;
        //reco'd eta meson using GsfElectrons
        TLorentzVector recoEtaVecE;
        //start out true, but if any leptons fail to be genmatched then it becomes false
        bool genmatched = true;
        //same but for GsfElectrons instead of PC tracks
        bool genmatchedE = true;
        for (size_t i = 0; i < genParticleHandle_->size(); i++) {
            reco::GenParticleRef genParticle(genParticleHandle_, i);
            // ?? what is this for??
            //if (!genParticle->isHardProcess()) continue;
            nt.genID_.push_back(genParticle->pdgId());
            nt.genHardProcess_.push_back(genParticle->isHardProcess());
            nt.genCharge_.push_back(genParticle->charge());
            nt.genPt_.push_back(genParticle->pt());
            nt.genEta_.push_back(genParticle->eta());
            nt.genPhi_.push_back(genParticle->phi());
            nt.genPz_.push_back(genParticle->pz());
            nt.genEn_.push_back(genParticle->energy());
            nt.genVxy_.push_back(genParticle->vertex().rho());
            //nt.genVxy_.push_back(genParticle->vxy()); //?????
            nt.genVz_.push_back(genParticle->vz());
            nt.genMass_.push_back(genParticle->mass());

            TLorentzVector genLepVec;
            genLepVec.SetPtEtaPhiM(genParticle->pt(), genParticle->eta(), genParticle->phi(), genParticle->mass());
            //skip particles that aren't supposed to be there
            if(genParticle->pdgId() == 990 || genParticle->pdgId() == 221) continue;
            if(i == 0) genEtaVec = genLepVec;
            else genEtaVec = genEtaVec + genLepVec;
            //4vector for the reco'd particle
            TLorentzVector recoLepVec;
            if(fabs(genParticle->pdgId()) == 11) {
                allGenPtEl->Fill(genParticle->pt());
                //for MC for electrons, do genmatching for allTracksP to find the best track corr. to the gen particle; fill the histogram
                    // so that later can calculate reco efficiency
                //TLorentzVector genpart;
                //genpart.SetPtEtaPhiM( genParticle->pt(), genParticle->eta(), genParticle->phi(), genParticle->mass() );
                //loop thru each track and find the closest dR match
                vector<reco::Track> alltracks = (genParticle->charge() > 0 ? allTracksP : allTracksN);
                float mindR = 9999.0;
                for(unsigned int k = 0; k < alltracks.size(); k++) {
                    reco::Track trk = alltracks[k];
                    float gendR = reco::deltaR(trk, *genParticle);
                    if( gendR < mindR ) {
                        mindR = gendR;
                        recoLepVec.SetPtEtaPhiM(trk.pt(), trk.eta(), trk.phi(), ele_mass);
                    }
                } //end loop over alltracks
                //now fill the histograms
                allGenPt->Fill(genParticle->pt());
                if(i == 0) recoEtaVec = recoLepVec; else recoEtaVec += recoLepVec;
                if(mindR < drCut) {
                    matchedGenPt->Fill(genParticle->pt());
                }
                else {
                    genmatched = false;
                }
                alldR->Fill(mindR);

                //now do the same for pat::Electrons instead of PC tracks.
                TLorentzVector recoLepVecE;
                mindR = 9999.0;
                //for (size_t i = 0; i < recoElectronHandle_->size(); i++) {
                for (size_t i = 0; i < recoLowPtElectronHandle_->size(); i++) {
                    //pat::ElectronRef ele(recoElectronHandle_, i);
                    pat::ElectronRef ele(recoLowPtElectronHandle_, i);
                    //make sure the Electron has the right charge!!
                    if(ele->charge() * (genParticle->charge()) < 0) continue;
                    float gendR = reco::deltaR(*ele, *genParticle);
                    if( gendR < mindR ) {
                        mindR = gendR;
                        recoLepVecE.SetPtEtaPhiM(ele->pt(), ele->eta(), ele->phi(), ele->mass());
                    }
                } //end loop over alltracks
                if(i == 0) recoEtaVecE = recoLepVecE; else recoEtaVecE += recoLepVecE;
                //now fill the histograms
                if(mindR < drCut) {
                    matchedGenPtE->Fill(genParticle->pt());
                    //std::cout << "Event " << (int)nt.eventNum_ << " : gen electron " << (int)i << " genmatched to gsfel (pdgId " << (int)genParticle->pdgId() << ")." << endl;
                }
                else {
                    genmatchedE = false;
                    //std::cout << "Event " << (int)nt.eventNum_ << " : gen electron " << (int)i << " NOT genmatched (pdgId " << (int)genParticle->pdgId() << ")." << endl;
                }
                alldRE->Fill(mindR);

            } //end is gen electron/positron
            //otherwise could be a muon
            else if(fabs(genParticle->pdgId()) == 13) {
                //for MC for electrons, do genmatching for allTracksP to find the best track corr. to the gen particle; fill the histogram
                    // so that later can calculate reco efficiency
                //loop thru each track and find the closest dR match
                float mindR = 9999.0;
                allGenPtMu->Fill(genParticle->pt());
                for (size_t i = 0; i < recoMuonHandle_->size(); i++) {
                    pat::MuonRef muo(recoMuonHandle_, i);
                    float gendR = reco::deltaR(*muo, *genParticle);
                    if( gendR < mindR ) {
                        mindR = gendR;
                        recoLepVec.SetPtEtaPhiM(muo->pt(), muo->eta(), muo->phi(), muo->mass());
                    }
                } //end loop over all muons
                if(i==0) recoEtaVec = recoLepVec; else recoEtaVec += recoLepVec;
                if(i==0) recoEtaVecE = recoLepVec; else recoEtaVecE += recoLepVec;
                //now fill the histograms
                if(mindR < drCut) {
                    matchedGenPtMu->Fill(genParticle->pt());
                }
                else {
                    genmatched = false;
                    genmatchedE = false;
                }
                alldRMu->Fill(mindR);

            } //end is gen muon
            
        }
        //now get the pt of the full eta meson; see if all 4 particles genmatched successfully AND reco invar mass in good range!
        allGenPtEta->Fill(genEtaVec.Pt());
        if(genmatched && recoEtaVec.M() > etaMassMin && recoEtaVec.M() < etaMassMax) {
            matchedGenPtEta->Fill(genEtaVec.Pt());
        }
        if(genmatchedE && recoEtaVecE.M() > etaMassMin && recoEtaVecE.M() < etaMassMax) {
            matchedGenPtEtaE->Fill(genEtaVec.Pt());
            //std::cout << "Event " << (int)nt.eventNum_ << " : eta meson fully genmatched!" << endl;
        }

        genT->Fill();
    }

    //std::cout << "computing vertices 0" << std::endl;
    //first get the 4-particle vertices, then get rid of all the particles/tracks that aren't involved in those.
    // EL-EL-MU-MU
    //VertexTracks primVertTrx = computeVertices(allTracksP, allTracksN, muTracksP, muTracksN, "mmee", theB, kvf);
    //VertexTracks primVertTrx = computeVertices(allTracksP, allTracksN, muonsP, muonsN, "mmee", theB, kvf, nt);
    //debugging
    ////std::cout << "AFTER computeVertices recoVtxTrackP_: ";
    ////for(uint8_t rvtp : nt.recoVtxTrackP_) cout << (int)rvtp << ", " ;
    ////std::cout << "AFTER computeVertices recoVtxTrackN_: ";
    ////for(uint8_t rvtn : nt.recoVtxTrackN_) cout << (int)rvtn << ", " ;
    ////std::cout << std::endl;
    ////std::cout << std::endl << "AFTER computeVertices recoVtxDr_: ";
    ////for(float rvdr : nt.recoVtxDr_) cout << rvdr << ", ";
    ////std::cout << std::endl;
    ////std::cout << "AFTER computeVertices recoVtxMuonP_: ";
    ////for(uint8_t rvtn : nt.recoVtxMuonP_) cout << (int)rvtn << ", " ;
    ////std::cout << std::endl;
    ////std::cout << "AFTER computeVertices recoVtxMuonN_: ";
    ////for(uint8_t rvtn : nt.recoVtxMuonN_) cout << (int)rvtn << ", " ;
    ////std::cout << std::endl;
    //// if the mmee vertex is no good, then no need to save the event!
    ////std::cout << "reassigning allTracksP and N" << std::endl;
    //std::cout << "Event " << (int)nt.eventNum_ << " nMuonsP: " << (int)nt.muonsP.size() << "; nMuonsN: " << (int)nt.muonsN.size() << std::endl;
    //std::cout << "     muonsP: "; for(auto mp : nt.muonsP) std::cout << (int)mp << ", "; std::cout << std::endl;
    //std::cout << "     muonsN: "; for(auto mn : nt.muonsN) std::cout << (int)mn << ", "; std::cout << std::endl;
    //allTracksP = primVertTrx.tracksP;
    //allTracksN = primVertTrx.tracksN;
    //*** commenting out these changes so that can use ALL muons for 2-lep vertexing!! ****//
    //muonsN = primVertTrx.muonsN;
    //muonsP = primVertTrx.muonsP;
    //****                             *************                                   **                ***//
    //std::cout << "Event " << (int)nt.eventNum_ << " nGoodMuonsP: " << (int)nt.muonsP.size() << "; nGoodMuonsN: " << (int)nt.muonsN.size() << std::endl;
    //std::cout << "     muonsP: "; for(auto mp : nt.muonsP) std::cout << (int)mp << ", "; std::cout << std::endl;
    //std::cout << "     muonsN: "; for(auto mn : nt.muonsN) std::cout << (int)mn << ", "; std::cout << std::endl;

    //maybe implement this later? but not needed for now.
    ////if there's a max number of vertices to save, delet the rest of them!
    //if(nt.maxNmmee > 0) {
    //    //new vectors for mmee vertex
    //    for(int nv=0; nv < nt.maxNmmee; nv++) {
    //        //find the vertex with the highest prob, add it to the new list
    //    }
    //}
    //std::cout << "computing vertices 1" << std::endl;
    //now get the vertices for just 2 GsfElectrons
    // EL-EL 
    computeVertices(elTracksP, elTracksN, "elel", theB, kvf, nt);
    //LowPtElectron-LowPtElectron
    computeVertices(elLowPtTracksP, elLowPtTracksN, "lplp", theB, kvf, nt);
    //computeKinematicVertices(elLowPtTracksP, elLowPtTracksN, "lplp", theB, kvf, nt);

    // MU-MU 
    VertexTracks primVertTrx = computeVertices(muonsP, muonsN, "mumu", theB, kvf, nt, pv);
    //save ONLY the muons that form good muon-muon vertices (otherwise will be too many)
    muonsN = primVertTrx.muonsN;
    muonsP = primVertTrx.muonsP;

    //mu-mu-el-el
    computeVertices(elTracksP, elTracksN, muonsP, muonsN, "mmelel", theB, kvf, nt);
    //mu-mu-LowPtElectron-LowPtElectron
    computeVertices(elLowPtTracksP, elLowPtTracksN, muonsP, muonsN, "mmlplp", theB, kvf, nt);

    //std::cout << "computing vertices 2" << std::endl;
    //lastly get the vertices for just 2 packed candidate tracks (electrons or pions)
    // PC-PC
    //cout << "all tracks" << std::endl;
    //computeVertices(allTracksP, allTracksN, "pcpc", theB, kvf, nt);


    //std::cout << "Event " << (int)nt.eventNum_ << " filled!" << std::endl;
    recoT->Fill();

    return;
}

void eta2mu2eAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void eta2mu2eAnalyzer::endJob() {
    fs->cd();
    if(!isData) {
        allGenPt->Write();
        matchedGenPt->Write();
        alldR->Write();
        matchedGenPtE->Write();
        alldRE->Write();
        allGenPtMu->Write();
        allGenPtEl->Write();
        matchedGenPtMu->Write();
        alldRMu->Write();
        allGenPtEta->Write();
        matchedGenPtEta->Write();
        matchedGenPtEtaE->Write();
    }
    //hdRP->Write();
    //hdRN->Write();
}

// define this as a plug-in
DEFINE_FWK_MODULE(eta2mu2eAnalyzer);
