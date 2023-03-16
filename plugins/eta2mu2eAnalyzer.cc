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

#include "TVectorD.h"    // for fixing tracks
#include "TMatrixDSym.h" // for fixing tracks

const float ele_mass = 0.000511; //GeV
const float mu_mass = 0.106; //GeV
const float pi_mass = 0.140; //GeV

struct VertexTracks {
    vector<reco::Track> tracksP;
    vector<reco::Track> tracksN;
    vector<reco::Track> muonsP;
    vector<reco::Track> muonsN;
};
reco::Track fix_track(const reco::Track *tk, double delta=1e-8);

reco::Track fix_track(const reco::TrackRef& tk)
{
    reco::Track t = reco::Track(*tk);
    return fix_track(&t);
}

/* Check for a not positive definite covariance matrix. If the covariance matrix is not positive definite, we force it to be positive definite by
 * adding the minimum eigenvalue to the diagonal of the covariance matrix plus `delta`.
 * See https://nhigham.com/2020/12/22/what-is-a-modified-cholesky-factorization/ */
reco::Track fix_track(const reco::Track *tk, double delta)
{
    unsigned int i, j;
    double min_eig = 1;

    /* Get the original covariance matrix. */
    reco::TrackBase::CovarianceMatrix cov = tk->covariance();

    /* Convert it from an SMatrix to a TMatrixD so we can get the eigenvalues. */
    TMatrixDSym new_cov(cov.kRows);
    for (i = 0; i < cov.kRows; i++) {
        for (j = 0; j < cov.kRows; j++) {
            /* Need to check for nan or inf, because for some reason these
             * cause a segfault when calling Eigenvectors().
             *
             * No idea what to do here or why this happens. */
            if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
                cov(i,j) = 1e-6;
            new_cov(i,j) = cov(i,j);
        }
    }

    /* Get the eigenvalues. */
    TVectorD eig(cov.kRows);
    new_cov.EigenVectors(eig);
    for (i = 0; i < cov.kRows; i++)
        if (eig(i) < min_eig)
            min_eig = eig(i);

    /* If the minimum eigenvalue is less than zero, then subtract it from the
     * diagonal and add `delta`. */
    if (min_eig < 0) {
        for (i = 0; i < cov.kRows; i++)
            cov(i,i) -= min_eig - delta;
    }

    return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
}
class eta2mu2eAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {

public:
    explicit eta2mu2eAnalyzer(const edm::ParameterSet&);
    ~eta2mu2eAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions&);
    
    typedef std::pair<std::string, float> IdPairf;
    typedef std::pair<std::string, bool> IdPairb;

private:
    //structure to hold the muons (p+n) and tracks (p+n) forming a good vertex

    bool getCollections(const edm::Event&);

    //virtual void beginJob() override;
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    ////void computeVertices(vector<T> &, vector<T> &, std::string type) override;
    //virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endJob() override;
    void beginJob() override;
    void beginRun(edm::Run const&, edm::EventSetup const&) override;
    float calcVertices(vector<reco::TransientTrack>, TransientVertex, std::string);
    //type T defined as template (so can take either GsfElectronRef or regular TrackRef)
    //template<typename T>
    VertexTracks computeVertices(vector<reco::Track> &, vector<reco::Track> &, std::string type, edm::ESHandle<TransientTrackBuilder>, KalmanVertexFitter);
    VertexTracks computeVertices(vector<reco::GsfTrackRef> &, vector<reco::GsfTrackRef> &, std::string type, edm::ESHandle<TransientTrackBuilder>, KalmanVertexFitter);
    VertexTracks computeVertices(vector<reco::Track> &, vector<reco::Track> &, vector<reco::Track> &, vector<reco::Track> &, std::string type, edm::ESHandle<TransientTrackBuilder>, KalmanVertexFitter);
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
    recoPhotonToken_(consumes<pat::PhotonCollection>(ps.getParameter<edm::InputTag>("photon_collection"))),
    esToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
    rhoToken_(consumes<double>(ps.getParameter<edm::InputTag>("rho"))),
    trkToken_(consumes<pat::PackedCandidateCollection>(ps.getParameter<edm::InputTag>("packed_candidate")))
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
    desc.add<edm::InputTag>("photon_collection", edm::InputTag("slimmedPhotons"));
    desc.add<edm::InputTag>("primary_vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
    desc.add<edm::InputTag>("genParticle", edm::InputTag("prunedGenParticles"));
    desc.add<edm::InputTag>("trigResult", edm::InputTag("TriggerResults", "", "HLT"));
    //desc.add<edm::InputTag>("trigEvent", edm::InputTag("hltTriggerSummaryAOD", "", "HLT"));
    desc.add<edm::InputTag>("pileups", edm::InputTag("slimmedAddPileupInfo"));
    desc.add<edm::InputTag>("genEvt", edm::InputTag("generator"));
    //desc.add<edm::InputTag>("transient_track_builder", edm::InputTag("TransientTrackBuilder"));
    desc.add<edm::InputTag>("rho", edm::InputTag("fixedGridRhoFastjetAll"));
    desc.add<edm::InputTag>("packed_candidate", edm::InputTag("packedPFCandidates"));
    
    descriptions.add("eta2mu2eAnalyzer", desc);
}

void eta2mu2eAnalyzer::beginJob()
{
    recoT = fs->make<TTree>("recoT", "recoT");
    nt.SetRecoTree(recoT);
    if (!isData) {
        genT = fs->make<TTree>("genT", "genT");
        nt.SetGenTree(genT);
    }
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

    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET120_PFMHT120_IDTight");
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET130_PFMHT130_IDTight");
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET140_PFMHT140_IDTight");
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight");
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight");
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET120_PFMHT120_IDTight_PFHT60"); // 2017+2018
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"); // 2017+2018
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET200_HBHECleaned"); // 2017+2018
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET200_HBHE_BeamHaloCleaned"); // 2017+2018
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned"); // 2017+2018
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFHT500_PFMET100_PFMHT100_IDTight"); // 2017+2018
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFHT700_PFMET85_PFMHT85_IDTight"); // 2017+2018
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFHT800_PFMET75_PFMHT75_IDTight"); // 2017+2018
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET170_HBHECleaned"); // 2016
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET300"); // 2016
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_MET200"); // 2016
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_PFHT300_PFMET110"); // 2016
    triggerPathsWithoutVersionNum_.emplace_back("HLT_IsoMu27"); // For MET trigger eff. studies in data
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight"); // Alternative triggers
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_DoubleMu3_DCA_PFMET50_PFMHT60"); // Alternative triggers
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_DoubleMu3_DZ_PFMET50_PFMHT60");  // Alternative triggers
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_DoubleMu3_DZ_PFMET70_PFMHT70");  // Alternative triggers
    //triggerPathsWithoutVersionNum_.emplace_back("HLT_DoubleMu3_DZ_PFMET90_PFMHT90");  // Alternative triggers
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

//helper function to calculate the vertices just given the transient tracks list
// returns probability of the vertex fit (-1 if the fit failed).
float eta2mu2eAnalyzer::calcVertices(vector<reco::TransientTrack> transient_tracks, TransientVertex tv, std::string type) {
    float vxy = -9999;
    float sigma_vxy = -9999;
    float vtx_chi2 = 999999;
    float vz = -9999;
    float prob = -1.0;

    if (tv.isValid()) {
        reco::Vertex vertex = reco::Vertex(tv);
        vxy = sqrt(vertex.x()*vertex.x() + vertex.y()*vertex.y());
        sigma_vxy = (1/vxy)*sqrt(vertex.x()*vertex.x()*vertex.xError()*vertex.xError() +
                vertex.y()*vertex.y()*vertex.yError()*vertex.yError());
        //sigma_vxy = (1/vxy)*(vertex.x()*vertex.xError() + vertex.y()*vertex.yError());
        vtx_chi2 = vertex.normalizedChi2();
        vz = vertex.z();

        //get probability for this vertex
        prob = TMath::Prob( vertex.chi2() , vertex.ndof() );
    }
    else {
        return prob;
    }

    //only save if prob > .1 ??
    if ( prob > 0.1 ) {
        nt.recoVtxReducedChi2_[type].push_back(vtx_chi2);
        nt.recoVtxVxy_[type].push_back(vxy);
        nt.recoVtxVz_[type].push_back(vz);
        nt.recoVtxSigmaVxy_[type].push_back(sigma_vxy);
        //std::cout << "probability: " << prob << std::endl;
    }
    return prob;
    
} 


//auto computeVertices = [&](vector<reco::TrackRef> coll_1, vector<reco::TrackRef> coll_2, std::string type) {
//auto computeVertices = [&](vector<reco::GsfTrackRef> coll_1, vector<reco::GsfTrackRef> coll_2, std::string type) {
// T is a template class that can be either a GsfTrackRef, or a reco::Track, or a reco::TrackRef, etc.
//        (anything that can go into the TransientTrackBuilder::build() method)
//template<typename T>
VertexTracks eta2mu2eAnalyzer::computeVertices(vector<reco::Track> & coll_1, vector<reco::Track> & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf) {
    VertexTracks myVertTracks;
    myVertTracks.tracksP = {};
    myVertTracks.tracksN = {};
    myVertTracks.muonsP = {};
    myVertTracks.muonsN = {};
    for (size_t i = 0; i < coll_1.size(); i++) {
        for (size_t j = 0; j < coll_2.size(); j++) {
            //if ( j > 15 || i > 15 ) continue;
            reco::Track part_i, part_j;
            part_i = coll_1[i];
            part_j = coll_2[j];

            float dr = -9999;
            TransientVertex tv;

            vector<reco::TransientTrack> transient_tracks{};
            transient_tracks.push_back(theB->build(fix_track(&part_i)));
            transient_tracks.push_back(theB->build(fix_track(&part_j)));
            tv = kvf.vertex(transient_tracks);
            float probVtx = calcVertices(transient_tracks, tv, type);
            //if ( probVtx > 0 ) {
            if ( probVtx > 0.1 ) {
                dr = reco::deltaR(part_i, part_j);
                nt.recoVtxDr_[type].push_back(dr);
                myVertTracks.tracksP.push_back(part_i);
                myVertTracks.tracksN.push_back(part_j);
            }
            
        } // j loop
    } // i loop
    return myVertTracks;
} // computeVertices

//overloaded function
VertexTracks eta2mu2eAnalyzer::computeVertices(vector<reco::GsfTrackRef> & coll_1, vector<reco::GsfTrackRef> & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf) {
    VertexTracks myVertTracks;
    myVertTracks.tracksP = {};
    myVertTracks.tracksN = {};
    myVertTracks.muonsP = {};
    myVertTracks.muonsN = {};
    for (size_t i = 0; i < coll_1.size(); i++) {
        for (size_t j = 0; j < coll_2.size(); j++) {
            //reco::TrackRef part_i, part_j;
            //reco::GsfTrackRef part_i, part_j;
            reco::GsfTrackRef part_i, part_j;
            //if (i < coll_1.size())
                part_i = coll_1[i];
            //if (j < coll_2.size())
                part_j = coll_2[j];

            float dr = -9999;
            TransientVertex tv;
            if (part_i.isNonnull() && part_j.isNonnull() ) { // && i != j) {
                vector<reco::TransientTrack> transient_tracks{};
                transient_tracks.push_back(theB->build(part_i));
                transient_tracks.push_back(theB->build(part_j));
                tv = kvf.vertex(transient_tracks);
                float probVtx = calcVertices(transient_tracks, tv, type);
                //if ( probVtx > 0 ) {
                if ( probVtx > 0.1 ) {
                    dr = reco::deltaR(*part_i, *part_j);
                    nt.recoVtxDr_[type].push_back(dr);
                }
            } 

        } // j loop
    } // i loop
    return myVertTracks;
} // computeVertices

//this version for fitting the 2 electron tracks and the 2 muons all together!
VertexTracks eta2mu2eAnalyzer::computeVertices(vector<reco::Track> & coll_1, vector<reco::Track> & coll_2, vector<reco::Track> & coll_3, vector<reco::Track> & coll_4, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf) {
    VertexTracks myVertTracks;
    myVertTracks.tracksP = {};
    myVertTracks.tracksN = {};
    myVertTracks.muonsP = {};
    myVertTracks.muonsN = {};
    //vectors of bools to save what particles have already been added to the collection to return
    vector<bool> igood {}; for(size_t i = 0; i < coll_1.size(); i++) igood.push_back(false);
    vector<bool> jgood {}; for(size_t i = 0; i < coll_2.size(); i++) jgood.push_back(false);
    vector<bool> kgood {}; for(size_t i = 0; i < coll_3.size(); i++) kgood.push_back(false);
    vector<bool> lgood {}; for(size_t i = 0; i < coll_4.size(); i++) lgood.push_back(false);
    unsigned int printevery = 100000;
    for (size_t i = 0; i < coll_1.size(); i++) {
        reco::Track part_i;
        part_i = coll_1[i];
        //suppose this is an electron
        TLorentzVector ele_i;
        ele_i.SetPtEtaPhiM(part_i.pt(),part_i.eta(),part_i.phi(), ele_mass);
        //now suppose this is a pion
        TLorentzVector pi_i;
        pi_i.SetPtEtaPhiM(part_i.pt(),part_i.eta(),part_i.phi(), pi_mass);
        //if ( i > 15 ) break; // ??
        for (size_t j = 0; j < coll_2.size(); j++) {
            //if ( j > 15 ) break; // ??
            //print a message once in a while to let us know it's still working
            //if ( (i*j) % printevery == (printevery-1) ) std::cout << "i=" << i << "/" << coll_1.size() << "; j=" << j << "/" << coll_2.size() << std::endl;
            reco::Track part_j;
            part_j = coll_2[j];
            TLorentzVector ele_j;
            ele_j.SetPtEtaPhiM(part_j.pt(),part_j.eta(),part_j.phi(), ele_mass);
            TLorentzVector pi_j;
            pi_j.SetPtEtaPhiM(part_j.pt(),part_j.eta(),part_j.phi(), pi_mass);
            for(size_t k = 0; k < coll_3.size(); k++) {
                reco::Track part_k;
                part_k = coll_3[k];
                TLorentzVector mu_k;
                mu_k.SetPtEtaPhiM(part_k.pt(),part_k.eta(),part_k.phi(), mu_mass);
                //if (! (part_k.isNonnull()) ) {
                //    continue;
                //}
                for(size_t l = 0; l < coll_4.size(); l++) {
                    reco::Track part_l;
                    part_l = coll_4[l];
                    TLorentzVector mu_l;
                    mu_l.SetPtEtaPhiM(part_l.pt(),part_l.eta(),part_l.phi(), mu_mass);
                    //if (! (part_l.isNonnull()) ) {
                    //    continue;
                    //}

                    //first compute the pT/M of the 4 particle system to make sure it's consistent with eta/eta' meson
                    TLorentzVector eta_mmee = ele_i + ele_j + mu_k + mu_l; 
                    TLorentzVector eta_mmpp = pi_i + pi_j + mu_k + mu_l; 
                    float eta_Mmmee = eta_mmee.M();
                    float eta_Mmmpp = eta_mmpp.M();
                    float eta_Ptmmee = eta_mmee.Pt();
                    float eta_Ptmmpp = eta_mmpp.Pt();
                    //if the invar. mass of the 4vector isn't in the range of interest, continue
                    bool passed_mmee = false;
                    bool passed_mmpp = false;
                    if( eta_Mmmee > 0.4 && eta_Mmmee < 1.1 && eta_Ptmmee > 12 )  {
                        passed_mmee = true;
                    }
                    if( eta_Mmmpp > 0.4 && eta_Mmmpp < 1.1 && eta_Ptmmpp > 12 )  {
                        passed_mmpp = true;
                    }
                    if( ! (passed_mmee || passed_mmpp) ) {
                        continue;
                    }

                    float dr = -9999;
                    TransientVertex tv;
                    vector<reco::TransientTrack> transient_tracks{};

                    // build all 4 transient tracks ( 2 el + 2 mu )
                    transient_tracks.push_back(theB->build(fix_track(&part_i)));
                    transient_tracks.push_back(theB->build(fix_track(&part_j)));
                    transient_tracks.push_back(theB->build(fix_track(&part_k)));
                    transient_tracks.push_back(theB->build(fix_track(&part_l)));

                    tv = kvf.vertex(transient_tracks);
                    float probVtx = calcVertices(transient_tracks, tv, type);
                    //if ( probVtx > 0 ) {
                    if ( probVtx > 0.1 ) {
                        //dr of the electron pair
                        dr = reco::deltaR(part_i, part_j);
                        nt.recoVtxDr_[type].push_back(dr);
                        if( passed_mmee ) {
                            nt.recoVtxPt_["mmee"].push_back(eta_Ptmmee);
                            nt.recoVtxM_["mmee"].push_back(eta_Mmmee);
                        }
                        if( passed_mmpp ) {
                            nt.recoVtxPt_["mmpp"].push_back(eta_Ptmmpp);
                            nt.recoVtxM_["mmpp"].push_back(eta_Mmmpp);
                        }
                        //check if the particle is saved yet or nah before adding it
                        if( !igood[i] ) {
                            myVertTracks.tracksP.push_back(part_i);
                            igood[i] = true;
                        }
                        if( !jgood[j] ) {
                            myVertTracks.tracksN.push_back(part_j);
                            jgood[j] = true;
                        }
                        if ( !kgood[k] ) {
                            myVertTracks.muonsP.push_back(part_k);
                            kgood[k] = true;
                        }
                        if( !lgood[l] ) {
                            myVertTracks.muonsN.push_back(part_l);
                            lgood[l] = true;
                        } 
                    }
                } //l loop
            } //k loop

        } // j loop
    } // i loop
    return myVertTracks;
} // computeVertices
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
//    // for now "good" electron means only ID is passed
//    // i.e. (IDmap % 2 == 1)
    nt.recoNElectron_ = recoElectronHandle_->size();
    nt.recoNGoodElectron_ = 0;
    //get the tracks from the electron objects
    vector<reco::GsfTrackRef> elTracksP{};    
    vector<reco::GsfTrackRef> elTracksN{};    

    nt.gsfNTrk_ = recoElectronHandle_->size();
    nt.gsfNGoodTrk_ = 0;
    //const edm::ValueMap<float> & eIDmap = *recoElectronIDHandle_;
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

        // or use cut-based ID instead of mva??: "cutBasedElectronID-Fall17-94X-V2-loose"
        //nt.recoElectronIDResult_.push_back( electronRef->electronID("cutBasedElectronID-Fall17-94X-V2-loose") );
        //nt.recoElectronIDResult_.push_back( electronRef->electronID("mvaEleID-Fall17-noIso-V2-wpLoose") );
        for(std::string id : nt.electronIDs) {
            nt.recoElectronIDResult_[id].push_back( electronRef->electronID(id) );
        }

        reco::GsfTrackRef eltrack = electronRef->gsfTrack();

        nt.gsfTrkPt_.push_back(eltrack->pt());
        nt.gsfTrkEta_.push_back(eltrack->eta());
        nt.gsfTrkPhi_.push_back(eltrack->phi());
        nt.gsfTrkCharge_.push_back(eltrack->charge());

        nt.gsfNGoodTrk_++;

        if ( electronRef->charge() > 0 ) {
            elTracksP.push_back(eltrack);
        }
        else {
            elTracksN.push_back(eltrack);
        }
    }
//
    // Also add all photons to ntuple, regardless of ID
    // Photon ID only produces 1 or 0
    nt.recoNPhoton_ = recoPhotonHandle_->size();
    nt.recoNGoodPhoton_ = 0;
    for (size_t i = 0; i < recoPhotonHandle_->size(); i++) {
        pat::PhotonRef photonRef(recoPhotonHandle_, i);
        nt.recoPhotonPt_.push_back(photonRef->pt());
        nt.recoPhotonEta_.push_back(photonRef->eta());
        nt.recoPhotonPhi_.push_back(photonRef->phi());
        for(std::string id : nt.photonIDs) {
            nt.recoPhotonIDResult_[id].push_back( (float) photonRef->photonID(id) );
        }
    }
    
    //vector of all Muon tracks (Positive and Negative separately)
    vector<reco::Track> muTracksP {};
    vector<reco::Track> muTracksN {};
    // Also add all muons to ntuple, regardless of ID
    nt.recoNMuon_ = recoMuonHandle_->size();
    nt.recoNGoodMuon_ = 0;
    for (size_t i = 0; i < recoMuonHandle_->size(); i++) {
        pat::MuonRef muonRef(recoMuonHandle_, i);
        nt.recoMuonPt_.push_back(muonRef->pt());
        nt.recoMuonEta_.push_back(muonRef->eta());
        nt.recoMuonPhi_.push_back(muonRef->phi());
        nt.recoMuonCharge_.push_back(muonRef->charge());
        if ( muonRef->charge() > 0 ) {
            muTracksP.push_back(* (muonRef->bestTrack()) );
        }
        else {
            muTracksN.push_back(* (muonRef->bestTrack()) );
        }
        //nt.recoMuonIDResult_.push_back(phIDmap[photonRef]);
        nt.recoMuonIDResult_.push_back( (float) (muonRef->muonID("All")) );
    }

    // Pick pair of muons with smallest vertex chi square fit for all collection combos
    edm::ESHandle<TransientTrackBuilder> theB = iSetup.getHandle(esToken_);
    //deprecated!
    //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
    KalmanVertexFitter kvf(true);

    //all tracks (will be mostly pions methinks)
    vector<reco::Track> allTracksP{};    //positive only
    vector<reco::Track> allTracksN{};    //negative only
    //for(size_t i = 0; i < trkHandle_->size(); i++) {
    //    reco::PackedCandidateRef packedRef(trkHandle_, i);
    //    reco::TrackRef trackRef = packedRef->Track();
    //    allTracks.push_back(trackRef);
    //}
    //now get the PackedCandidate tracks
    nt.recoNTrk_ = trkHandle_->size();
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
        nt.recoNGoodTrk_++;
        nt.recoTrkPt_.push_back(iTrack1->pt());
        nt.recoTrkEta_.push_back(iTrack1->eta());
        nt.recoTrkPhi_.push_back(iTrack1->phi());
        nt.recoTrkCharge_.push_back(iTrack1->charge());
    }

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

    // EL-EL-MU-MU
    VertexTracks primVertTrx = computeVertices(allTracksP, allTracksN, muTracksP, muTracksN, "mmee", theB, kvf);
    // if the mmee vertex is no good, then no need to save the event!
    allTracksP = primVertTrx.tracksP;
    allTracksN = primVertTrx.tracksN;
    nt.recoNGoodTrk_ = allTracksP.size() + allTracksN.size();

    // EL-EL 
    computeVertices(elTracksP, elTracksN, "elel", theB, kvf);
    // PC-PC
    computeVertices(allTracksP, allTracksN, "pcpc", theB, kvf);

    /****** GEN INFO *******/

    if (!isData) {

        nt.nGen_ = (int)genParticleHandle_->size();
        
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

        genT->Fill();
    }

    recoT->Fill();

    return;
}

void eta2mu2eAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void eta2mu2eAnalyzer::endJob() {}

// define this as a plug-in
DEFINE_FWK_MODULE(eta2mu2eAnalyzer);
