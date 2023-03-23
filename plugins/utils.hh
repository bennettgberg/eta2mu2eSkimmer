#ifndef __UTILS__
#define __UTILS__

#include <vector>
using std::vector;
#include <TLorentzVector.h>

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TVectorD.h"    // for fixing tracks
#include "TMatrixDSym.h" // for fixing tracks

#include "NtupleContainer.hh"

const float ele_mass = 0.000511; //GeV
const float mu_mass = 0.106; //GeV
const float pi_mass = 0.140; //GeV


float calcVertices(vector<reco::TransientTrack> transient_tracks, TransientVertex tv, std::string type, NtupleContainer & nt);

//struct for saving the used tracks and muons
// -- used as return type for computeVertices functions
struct VertexTracks {
    vector<reco::Track> tracksP;
    vector<reco::Track> tracksN;
    vector<pat::Muon> muonsP;
    vector<pat::Muon> muonsN;
};

// compute vertices for two track vectors coll_1 and coll_2, and add them to the ntuple.
VertexTracks computeVertices(vector<reco::Track> & coll_1, vector<reco::Track> & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt);

// compute vertex for two sets of GsfElectron vectors coll_1 and coll_2.
VertexTracks computeVertices(vector<reco::GsfTrackRef> & coll_1, vector<reco::GsfTrackRef> & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt);

// this version for fitting the 2 electron tracks and the 2 muons all together!
VertexTracks computeVertices(vector<reco::Track> & coll_1, vector<reco::Track> & coll_2, vector<pat::Muon> & coll_3, vector<pat::Muon> & coll_4, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt);

//fix tracks with non-pos-def covariance matrices -- needed to prevent crashing
reco::Track fix_track(const reco::Track *tk, double delta=1e-8);

//got this code from Sergey Polikarpov
reco::Track fix_track(const reco::TrackRef& tk);

#endif // __UTILS__
