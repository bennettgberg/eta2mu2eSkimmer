#include <algorithm>
#include <cmath> 
#include <memory>
#include <random>
#include <vector>
#include <boost/format.hpp>
#include <boost/any.hpp>

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
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
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

#include "TTree.h"
#include "TMath.h"
#include <TLorentzVector.h>

#include "NtupleContainer.hh"
#include "utils.hh"

class TMAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {

public:
    explicit TMAnalyzer(const edm::ParameterSet&);
    ~TMAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions&);
    
    typedef std::pair<std::string, float> IdPairf;
    typedef std::pair<std::string, bool> IdPairb;

private:
    bool getCollections(const edm::Event&);

    void beginJob() override;
    void beginRun(edm::Run const&, edm::EventSetup const&) override;
    float calcVertices(vector<reco::TransientTrack>, TransientVertex, std::string);
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endRun(edm::Run const&, edm::EventSetup const&) override;
    void endJob() override;

    TTree *recoT, *genT;
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
    const edm::EDGetTokenT<pat::PhotonCollection> recoPhotonToken_;
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> esToken_;
    const edm::EDGetTokenT<double> rhoToken_;
    const edm::EDGetTokenT<pat::PackedCandidateCollection> trkToken_;

    // Handles
    edm::Handle<pat::MuonCollection> recoMuonHandle_;
    edm::Handle<reco::VertexCollection> primaryVertexHandle_;
    edm::Handle<reco::GenParticleCollection> genParticleHandle_;
    edm::Handle<edm::TriggerResults> trigResultsHandle_;
    //edm::Handle<trigger::TriggerEvent> trigEventHandle_;
    edm::Handle<std::vector<PileupSummaryInfo>> pileupInfosHandle_;
    edm::Handle<GenEventInfoProduct> genEvtInfoHandle_;
    edm::Handle<pat::ElectronCollection> recoElectronHandle_;
    edm::Handle<pat::PhotonCollection> recoPhotonHandle_;
    edm::Handle<double> rhoHandle_;
    edm::Handle<pat::PackedCandidateCollection> trkHandle_;
    
    std::vector<std::string> triggerPathsWithoutVersionNum_;
    std::vector<std::string> triggerPathsWithVersionNum_;
    std::vector<bool> trigExist_;
    HLTConfigProvider hltConfig_;

};


TMAnalyzer::TMAnalyzer(const edm::ParameterSet& ps):
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
    recoPhotonToken_(consumes<pat::PhotonCollection>(ps.getParameter<edm::InputTag>("photon_collection"))),
    esToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
    rhoToken_(consumes<double>(ps.getParameter<edm::InputTag>("rho"))),
    trkToken_(consumes<pat::PackedCandidateCollection>(ps.getParameter<edm::InputTag>("packed_candidate")))
{
    usesResource("TFileService");
    m_random_generator = std::mt19937(37428479);
}

TMAnalyzer::~TMAnalyzer() = default;

void TMAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    // Only specify tags with reasonable defaults -- check out cfg for others

    edm::ParameterSetDescription desc;
    desc.add<bool>("isData", 0);
    desc.add<std::string>("triggerProcessName", "HLT");

    desc.add<edm::InputTag>("muon_collection", edm::InputTag("slimmedMuons"));
    desc.add<edm::InputTag>("electron_collection", edm::InputTag("slimmedElectrons"));
    desc.add<edm::InputTag>("photon_collection", edm::InputTag("slimmedPhotons"));
    desc.add<edm::InputTag>("primary_vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
    desc.add<edm::InputTag>("genParticle", edm::InputTag("prunedGenParticles"));
    desc.add<edm::InputTag>("trigResult", edm::InputTag("TriggerResults", "", "HLT"));
    //desc.add<edm::InputTag>("trigEvent", edm::InputTag("hltTriggerSummaryAOD", "", "HLT"));
    desc.add<edm::InputTag>("pileups", edm::InputTag("slimmedAddPileupInfo"));
    desc.add<edm::InputTag>("genEvt", edm::InputTag("generator"));
    desc.add<edm::InputTag>("rho", edm::InputTag("fixedGridRhoFastjetAll"));
    desc.add<edm::InputTag>("packed_candidate", edm::InputTag("packedPFCandidates"));
    
    descriptions.add("TMAnalyzer", desc);
}

void TMAnalyzer::beginJob()
{
    recoT = fs->make<TTree>("recoT", "recoT");
    nt.SetRecoTree(recoT);
    if (!isData) {
        //genT = fs->make<TTree>("genT", "genT");
        nt.SetGenTree(recoT);
    }
    nt.CreateTreeBranches();
}


void TMAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    using namespace edm;

    bool changed = true;
    if (hltConfig_.init(iRun, iSetup, triggerProcessName_, changed)) {
        if (changed) {
            LogInfo("HLTConfig") << "TMAnalyzer::beginRun: " << "hltConfig init for Run" << iRun.run();
            hltConfig_.dump("ProcessName");
            hltConfig_.dump("GlobalTag");
            hltConfig_.dump("TableName");
        }
    } 
    else {
        LogError("HLTConfig") << "TMAnalyzer::beginRun: config extraction failure with triggerProcessName -> " << triggerProcessName_;
        return;
    }

    // Add trigger paths if they exist
    triggerPathsWithoutVersionNum_.clear();
    triggerPathsWithVersionNum_.clear();
    trigExist_.clear();

    triggerPathsWithoutVersionNum_.emplace_back("HLT_IsoMu27"); // For MET trigger eff. studies in data
    triggerPathsWithoutVersionNum_.emplace_back("HLT_L2Mu10_NoVertex_NoBPTX");    // For dSA eff. studies in data
    triggerPathsWithoutVersionNum_.emplace_back("HLT_L2Mu10_NoVertex_NoBPTX3BX"); // For dSA eff. studies in data
    
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
        }
    }
}


bool TMAnalyzer::getCollections(const edm::Event& iEvent) {
    using namespace edm;

    char error_msg[] = "TMAnalyzer::GetCollections: Error in getting product %s from Event!";

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
    getHandle(recoPhotonToken_, recoPhotonHandle_, "photon_collection");
    getHandle(rhoToken_, rhoHandle_, "rho");
    getHandle(trkToken_, trkHandle_, "packed_candidate");
    if (!isData) {
        getHandle(genEvtInfoToken_, genEvtInfoHandle_, "genEventInfo");
        getHandle(genParticleToken_, genParticleHandle_, "genParticle");
        getHandle(pileupInfosToken_, pileupInfosHandle_, "pileupInfos");
    }
    
    return ret;
}

void TMAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using std::cout, std::vector, std::endl;

    if (!getCollections(iEvent))
        return;

    // Clear branches before filling
    nt.ClearTreeBranches();

    // Start filling

    // Event information
    nt.eventNum_ = iEvent.id().event();
    nt.lumiSec_ = iEvent.luminosityBlock();
    nt.runNum_ = iEvent.id().run();
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
    for (size_t i = 0; i < recoElectronHandle_->size(); i++) {
        pat::ElectronRef electronRef(recoElectronHandle_, i);
        nt.recoElectronPt_.push_back(electronRef->pt());
        nt.recoElectronEta_.push_back(electronRef->eta());
        nt.recoElectronPhi_.push_back(electronRef->phi());
        nt.recoElectronVxy_.push_back(electronRef->trackPositionAtVtx().rho()); //?????
        nt.recoElectronVz_.push_back(electronRef->trackPositionAtVtx().z());
        nt.recoElectronCharge_.push_back(electronRef->charge());
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

        nt.gsfTrkPt_.push_back(eltrack->pt());
        nt.gsfTrkEta_.push_back(eltrack->eta());
        nt.gsfTrkPhi_.push_back(eltrack->phi());
        nt.gsfTrkCharge_.push_back(eltrack->charge());


        if ( electronRef->charge() > 0 ) {
            elTracksP.push_back(eltrack);
            nt.gsfElsP.push_back(nt.recoNGoodElectron_);
        }
        else {
            elTracksN.push_back(eltrack);
            nt.gsfElsN.push_back(nt.recoNGoodElectron_);
        }
        nt.recoNGoodElectron_++;
    }

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
        //nt.recoMuonIDResult_.push_back( (float) (muonRef->muonID("All")) );
        //nt.recoNGoodMuon_++;
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
    // andre: don't need tracks for true muonium
    //
    //for (pat::PackedCandidateCollection::const_iterator iTra1 = trkHandle_->begin(); iTra1 != trkHandle_->end(); iTra1++) {
    //    //try setting track mass to electron mass to see what happens
    //    //iTra1->setMass(ele_mass); //doesn't work
    //    if (!(iTra1->hasTrackDetails())) continue;

    //    //electrons only
    //    //float mass = iTra1->mass();
    //    //if( fabs(mass - ele_mass) > (ele_mass/10.0) ) continue;
    //    if( iTra1->charge() == 0 ) continue;
    //    //std::cout << "mass: " << mass << std::endl;

    //    auto iTrack1 = iTra1->bestTrack();
    //    if (iTrack1 == nullptr) continue;
    //    //
    //    if (iTrack1->pt()<0.5) continue;
    //    //if (iTrack1->charge() < 0.5) continue; // positive
    //    if (!(iTrack1->quality(reco::TrackBase::highPurity))) continue;
    //    //
    //    if( fabs(iTrack1->eta()) > 2.9 ) continue;
    //    //
    //    //TLorentzVector p4el1;
    //    //p4el1.SetPtEtaPhiM(iTrack1->pt(),iTrack1->eta(),iTrack1->phi(), ele_mass);
    //    //
    //    reco::TransientTrack elec1TT = theB->build( fix_track(iTra1->bestTrack()) );
    //    //reco::TransientTrack elec1TT = theB->build( fix_track(iTra1->bestTrack()) );
    //    if (!elec1TT.isValid()) continue;
 
    //    if ( iTrack1->charge() > 0.5) {
    //        allTracksP.push_back(*iTrack1);
    //    }
    //    else if ( iTrack1->charge() < -0.5 ) {
    //        allTracksN.push_back(*iTrack1); 
    //    }
    //    else {
    //        continue;
    //    }
    //    //nt.recoNGoodTrk_++;
    //    //nt.recoTrkPt_.push_back(iTrack1->pt());
    //    //nt.recoTrkEta_.push_back(iTrack1->eta());
    //    //nt.recoTrkPhi_.push_back(iTrack1->phi());
    //    //nt.recoTrkCharge_.push_back(iTrack1->charge());
    //}

// TODO: implement triggers
//
//    // Assign each trigger result to a different bit
//    nt.fired_ = 0;
//    for (size_t i = 0; i < triggerPathsWithVersionNum_.size(); i++) {
//        if (trigExist_.at(i)) {
//            std::string trigPath = triggerPathsWithVersionNum_[i];
//            nt.fired_ |= (trigResultsHandle_->accept(hltConfig_.triggerIndex(trigPath)) << i);
//        }
//        else {
//            nt.fired_ |= (0 <<i);
//        }
//    }
//    

    //first get the 4-particle vertices, then get rid of all the particles/tracks that aren't involved in those.
    // EL-EL-MU-MU
    //VertexTracks primVertTrx = computeVertices(allTracksP, allTracksN, muTracksP, muTracksN, "mmee", theB, kvf);
    //VertexTracks primVertTrx = computeVertices(allTracksP, allTracksN, muonsP, muonsN, "mmee", theB, kvf, nt);
    // if the mmee vertex is no good, then no need to save the event!
    //allTracksP = primVertTrx.tracksP;
    //allTracksN = primVertTrx.tracksN;
    //muonsN = primVertTrx.muonsN;
    //muonsP = primVertTrx.muonsP;
    for(auto muonRef : muonsP) {
        nt.recoMuonPt_.push_back(muonRef.pt());
        nt.recoMuonEta_.push_back(muonRef.eta());
        nt.recoMuonPhi_.push_back(muonRef.phi());
        nt.recoMuonCharge_.push_back(muonRef.charge());
        nt.recoNGoodMuon_++;
    }
    for(auto muonRef : muonsN) {
        nt.recoMuonPt_.push_back(muonRef.pt());
        nt.recoMuonEta_.push_back(muonRef.eta());
        nt.recoMuonPhi_.push_back(muonRef.phi());
        nt.recoMuonCharge_.push_back(muonRef.charge());
        nt.recoNGoodMuon_++;
    }
  
    //for(auto iTrack1 : allTracksP) {
    //    nt.recoNGoodTrk_++;
    //    nt.recoTrkPt_.push_back(iTrack1.pt());
    //    nt.recoTrkEta_.push_back(iTrack1.eta());
    //    nt.recoTrkPhi_.push_back(iTrack1.phi());
    //    nt.recoTrkCharge_.push_back(iTrack1.charge());
    //}
    //for(auto iTrack1 : allTracksN) {
    //    nt.recoNGoodTrk_++;
    //    nt.recoTrkPt_.push_back(iTrack1.pt());
    //    nt.recoTrkEta_.push_back(iTrack1.eta());
    //    nt.recoTrkPhi_.push_back(iTrack1.phi());
    //    nt.recoTrkCharge_.push_back(iTrack1.charge());
    //}

    //now get the vertices for just 2 muons
    // MU-MU 
    computeVertices(muonsP, muonsN, "mumu", theB, kvf, nt);
    //lastly get the vertices for just 2 packed candidate tracks (electrons or pions)
    // PC-PC
    //computeVertices(allTracksP, allTracksN, "pcpc", theB, kvf, nt);

    /****** GEN INFO *******/

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
        }

        //genT->Fill();
    }

    // only add event if at least one combination of mumugamma has a mass close to the eta
    bool saveEvent = false;
    for (auto muonRefP : muonsP) {
      if (saveEvent) break;
      for (auto muonRefN : muonsN) {
        if (saveEvent) break;
        for (size_t i = 0; i < recoPhotonHandle_->size(); i++) {
          pat::PhotonRef photonRef(recoPhotonHandle_, i);
          TLorentzVector p_mu_plus(muonRefP.pt(), muonRefP.eta(), muonRefP.phi(), mu_mass);
          TLorentzVector p_mu_minus(muonRefN.pt(), muonRefN.eta(), muonRefN.phi(), mu_mass);
          TLorentzVector p_photon(photonRef->pt(), photonRef->eta(), photonRef->phi(), 0.);
          TLorentzVector sum = p_mu_plus + p_mu_minus + p_photon;
          if (sum.M() > 0.3 && sum.M() < 0.7) {
            saveEvent = true;
            break;
          }
        }
      }
    }

    if (saveEvent)
      recoT->Fill();

    return;
}

void TMAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void TMAnalyzer::endJob() {}

// define this as a plug-in
DEFINE_FWK_MODULE(TMAnalyzer);
