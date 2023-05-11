#include "utils.hh"

//helper function to calculate the vertices just given the transient tracks list
// returns probability of the vertex fit (-1 if the fit failed).
float calcVertices(vector<reco::TransientTrack> transient_tracks, TransientVertex tv, std::string type, NtupleContainer & nt) {
    float vxy = -9999;
    float sigma_vxy = -9999;
    float vtx_chi2 = 999999;
    float vz = -9999;
    float prob = -1.0;

    //only process valid transient vertices
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

    //only save if prob > .1 
    //if ( prob > 0.1 ) {
    //    nt.recoVtxReducedChi2_[type].push_back(vtx_chi2);
    //    nt.recoVtxVxy_[type].push_back(vxy);
    //    nt.recoVtxVz_[type].push_back(vz);
    //    nt.recoVtxSigmaVxy_[type].push_back(sigma_vxy);
    //    //std::cout << "probability: " << prob << std::endl;
    //}
    //TODO: keep in this || elel?? decide.
    if ( prob > 0.1 ) {
        if(type == "mmee") {
            nt.mmeeVtxReducedChi2_.push_back(vtx_chi2);
            nt.mmeeVtxVxy_.push_back(vxy);
            nt.mmeeVtxVz_.push_back(vz);
            nt.mmeeVtxSigmaVxy_.push_back(sigma_vxy);
        }
        else if(type == "elel") {
            nt.elelVtxReducedChi2_.push_back(vtx_chi2);
            nt.elelVtxVxy_.push_back(vxy);
            nt.elelVtxVz_.push_back(vz);
            nt.elelVtxSigmaVxy_.push_back(sigma_vxy);
        }
        else if(type == "lplp") {
            nt.lplpVtxReducedChi2_.push_back(vtx_chi2);
            nt.lplpVtxVxy_.push_back(vxy);
            nt.lplpVtxVz_.push_back(vz);
            nt.lplpVtxSigmaVxy_.push_back(sigma_vxy);
        }
        else if(type == "pcpc") {
            nt.pcpcVtxReducedChi2_.push_back(vtx_chi2);
            nt.pcpcVtxVxy_.push_back(vxy);
            nt.pcpcVtxVz_.push_back(vz);
            nt.pcpcVtxSigmaVxy_.push_back(sigma_vxy);
        }
        else if(type == "mumu") {
            nt.mumuVtxReducedChi2_.push_back(vtx_chi2);
            nt.mumuVtxVxy_.push_back(vxy);
            nt.mumuVtxVz_.push_back(vz);
            nt.mumuVtxSigmaVxy_.push_back(sigma_vxy);
        }
        else if(type == "mmelel") {
            nt.mmelelVtxReducedChi2_.push_back(vtx_chi2);
            nt.mmelelVtxVxy_.push_back(vxy);
            nt.mmelelVtxVz_.push_back(vz);
            nt.mmelelVtxSigmaVxy_.push_back(sigma_vxy);
        }
        else if(type == "mmlplp") {
            nt.mmlplpVtxReducedChi2_.push_back(vtx_chi2);
            nt.mmlplpVtxVxy_.push_back(vxy);
            nt.mmlplpVtxVz_.push_back(vz);
            nt.mmlplpVtxSigmaVxy_.push_back(sigma_vxy);
        }
        //std::cout << "probability: " << prob << std::endl;
    }
    return prob;
    
} 

void computeVertex(pat::Muon & coll_1, pat::Muon & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt) {
  reco::Track part_1, part_2;
  part_1 = *(coll_1.bestTrack());
  part_2 = *(coll_2.bestTrack());
  //first build the transient vertex and transient tracks.
  float dr = -9999;
  TransientVertex tv;

  vector<reco::TransientTrack> transient_tracks{};
  transient_tracks.push_back(theB->build(fix_track(&part_1)));
  transient_tracks.push_back(theB->build(fix_track(&part_2)));
  tv = kvf.vertex(transient_tracks);
  float probVtx = calcVertices(transient_tracks, tv, type, nt);
  if ( probVtx > 0.1 ) {
    dr = reco::deltaR(part_1, part_2);
    nt.mumuVtxDr_.push_back(dr);
  }
}


//compute vertices for two muons vectors coll_1 and coll_2, and add them to the ntuple.
// (overloaded below)
// Return type: struct VertexTracks (defined above)
void computeVertices(vector<pat::Muon> & coll_1, vector<pat::Muon> & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt) {
    for (size_t i = 0; i < coll_1.size(); i++) {
        for (size_t j = 0; j < coll_2.size(); j++) {
            //if ( j > 15 || i > 15 ) continue;
            reco::Track part_i, part_j;
            part_i = *(coll_1[i].bestTrack());
            part_j = *(coll_2[j].bestTrack());

            //first build the transient vertex and transient tracks.
            float dr = -9999;
            TransientVertex tv;

            vector<reco::TransientTrack> transient_tracks{};
            transient_tracks.push_back(theB->build(fix_track(&part_i)));
            transient_tracks.push_back(theB->build(fix_track(&part_j)));
            tv = kvf.vertex(transient_tracks);
            float probVtx = calcVertices(transient_tracks, tv, type, nt);
            //only fill the ntuple if the (chi2) prob is > .1
            // the rest of the ntuple info is filled in the calcVertices function (above)
            //if ( probVtx > 0 ) {
            //if ( probVtx > 0.1 ) {
            //    dr = reco::deltaR(part_i, part_j);
            //    nt.recoVtxDr_[type].push_back(dr);
            //}
            if ( probVtx > 0.1 ) {
                dr = reco::deltaR(part_i, part_j);
                nt.mumuVtxDr_.push_back(dr);
            }
            
        } // j loop
    } // i loop
} // computeVertices

//compute vertices for two track vectors coll_1 and coll_2, and add them to the ntuple.
// (overloaded below)
// Return type: struct VertexTracks (defined above)
void computeVertices(vector<reco::Track> & coll_1, vector<reco::Track> & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt) {
    for (size_t i = 0; i < coll_1.size(); i++) {
        for (size_t j = 0; j < coll_2.size(); j++) {
            //if ( j > 15 || i > 15 ) continue;
            reco::Track part_i, part_j;
            part_i = coll_1[i];
            part_j = coll_2[j];

            //first build the transient vertex and transient tracks.
            float dr = -9999;
            TransientVertex tv;

            vector<reco::TransientTrack> transient_tracks{};
            transient_tracks.push_back(theB->build(fix_track(&part_i)));
            transient_tracks.push_back(theB->build(fix_track(&part_j)));
            tv = kvf.vertex(transient_tracks);
            float probVtx = calcVertices(transient_tracks, tv, type, nt);
            //only fill the ntuple if the (chi2) prob is > .1
            // the rest of the ntuple info is filled in the calcVertices function (above)
            //if ( probVtx > 0 ) {
            if ( probVtx > 0.1 ) {
            //if ( 0.2 > 0.1 ) {
                dr = reco::deltaR(part_i, part_j);
                //nt.recoVtxDr_[type].push_back(dr);
                nt.pcpcVtxDr_.push_back(dr);
                //there are positive AND negative tracks
                // number pushed back is 4 bits for the first track # (in the good list), followed by 4 bits for the 2nd track #
                //  so max we can handle is 16 good tracks...
                //uint8_t full_val = nt.mmeeTrxP[i] + (nt.mmeeTrxN[j] << 4);
                uint8_t trackP = nt.mmeeTrxP[i];
                uint8_t trackN = nt.mmeeTrxN[j];
                //nt.recoVtxTracks_[type].push_back(full_val);
                //nt.recoVtxTrackP_[type].push_back(trackP);
                //nt.recoVtxTrackN_[type].push_back(trackN);
                nt.pcpcVtxTrackP_.push_back(trackP);
                nt.pcpcVtxTrackN_.push_back(trackN);
            }
            
        } // j loop
    } // i loop
} // computeVertices

//overloaded function
// now compute vertex for two sets of GsfElectron vectors coll_1 and coll_2.
void computeVertices(vector<reco::GsfTrackRef> & coll_1, vector<reco::GsfTrackRef> & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt) {
    for (size_t i = 0; i < coll_1.size(); i++) {
        for (size_t j = 0; j < coll_2.size(); j++) {
            //std::cout << "Vertex cand " << (int)i << " " << (int)j << std::endl;
            //reco::TrackRef part_i, part_j;
            //reco::GsfTrackRef part_i, part_j;
            reco::GsfTrackRef part_i, part_j;
            part_i = coll_1[i];
            part_j = coll_2[j];

            float dr = -9999;
            TransientVertex tv;
            if (part_i.isNonnull() && part_j.isNonnull() ) { // && i != j) {
                vector<reco::TransientTrack> transient_tracks{};
                transient_tracks.push_back(theB->build(part_i));
                transient_tracks.push_back(theB->build(part_j));
                tv = kvf.vertex(transient_tracks);
                float probVtx = calcVertices(transient_tracks, tv, type, nt);
                //if ( probVtx > 0 ) {
                if ( probVtx > 0.1 ) {
                //if ( 0.2 > 0.1 ) {
                    dr = reco::deltaR(*part_i, *part_j);
                    if(type == "elel") {
                        nt.elelVtxDr_.push_back(dr);
                        //first 4 bits for first electron, last 4 for second electron
                        //uint8_t full_val = nt.gsfElsP[i] + (nt.gsfElsN[j] << 4);
                        uint8_t eleP = nt.gsfElsP[i];
                        uint8_t eleN = nt.gsfElsN[j];
                        //std::cout << "eleP: " << (int)eleP << "; eleN: " << (int)eleN << std::endl;
                        //std::cout << "full vector gsfElsP: ";
                        //for(uint8_t elp : nt.gsfElsP) std::cout << (int)elp << ", ";
                        //std::cout << std::endl;
                        //std::cout << "full vector gsfElsN: ";
                        //for(uint8_t eln : nt.gsfElsN) std::cout << (int)eln << ", ";
                        //std::cout << std::endl;
                        //std::cout << "nGoodElectron: " << (int)nt.recoNGoodElectron_ << std::endl;
                        //nt.recoVtxEleP_[type].push_back(eleP);
                        //nt.recoVtxEleN_[type].push_back(eleN);
                        nt.elelVtxEleP_.push_back(eleP);
                        nt.elelVtxEleN_.push_back(eleN);
                        //std::cout << "eleP: " << (int)eleP << "; eleN: " << (int)eleN;
                    } //end elel type
                    else if(type == "lplp") {
                        nt.lplpVtxDr_.push_back(dr);
                        uint8_t eleP = nt.gsfLowPtElsP[i];
                        uint8_t eleN = nt.gsfLowPtElsN[j];
                        nt.lplpVtxEleP_.push_back(eleP);
                        nt.lplpVtxEleN_.push_back(eleN);
                    } //end lplp type
                    else {
                        std::cout << "Error!! Unrecognized vertex type " << type << std::endl;
                        throw std::runtime_error("Error: unrecognized type for 2-lepton vertex.");
                    }
                }
            } //end is Nonnull 
            else { //is null ....
                std::cout << "Null GsfTrack! Event" << (int)nt.eventNum_ << " ";
                if(!part_i.isNonnull()) std::cout << "part i: " << (int)i << " "; 
                if(!part_j.isNonnull()) std::cout << "part j: " << (int)j << " "; 
                std::cout << std::endl;
                if(type == "elel") {
                    nt.elelVtxDr_.push_back(9999);
                    uint8_t eleP = nt.gsfElsP[i];
                    uint8_t eleN = nt.gsfElsN[j];
                    nt.elelVtxEleP_.push_back(eleP);
                    nt.elelVtxEleN_.push_back(eleN);
                    nt.elelVtxReducedChi2_.push_back(9999);
                    nt.elelVtxVxy_.push_back(9999);
                    nt.elelVtxVz_.push_back(9999);
                    nt.elelVtxSigmaVxy_.push_back(9999);
                } //elel type
                else if(type == "lplp") {
                    nt.lplpVtxDr_.push_back(9999);
                    uint8_t eleP = nt.gsfLowPtElsP[i];
                    uint8_t eleN = nt.gsfLowPtElsN[j];
                    nt.lplpVtxEleP_.push_back(eleP);
                    nt.lplpVtxEleN_.push_back(eleN);
                    nt.lplpVtxReducedChi2_.push_back(9999);
                    nt.lplpVtxVxy_.push_back(9999);
                    nt.lplpVtxVz_.push_back(9999);
                    nt.lplpVtxSigmaVxy_.push_back(9999);
                }
            }
        } // j loop
    } // i loop
} // computeVertices

//add a track to the ntuple
NtupleContainer addTrack(NtupleContainer nt, reco::Track iTrack1) {
    //can't add another if we're at the maximum
    if(nt.recoNGoodTrk_ == 255) {
        throw std::runtime_error("Too many good tracks!");
    }
    nt.recoTrkPt_.push_back(iTrack1.pt());
    nt.recoTrkEta_.push_back(iTrack1.eta());
    nt.recoTrkPhi_.push_back(iTrack1.phi());
    nt.recoTrkCharge_.push_back(iTrack1.charge());
    nt.recoNGoodTrk_++; 
    //std::cout << "just incremented recoNGoodTrk: " << (int)nt.recoNGoodTrk_ << std::endl;
    return nt;
}

//add a muon to the ntuple
NtupleContainer addMuon(NtupleContainer nt, pat::Muon iMuon1) {
    if(nt.recoNGoodMuon_ == 255) {
        throw std::runtime_error("Too many good muons!");
    }
    nt.recoMuonPt_.push_back(iMuon1.pt());
    nt.recoMuonEta_.push_back(iMuon1.eta());
    nt.recoMuonPhi_.push_back(iMuon1.phi());
    nt.recoMuonCharge_.push_back(iMuon1.charge());
    nt.recoNGoodMuon_++; 
    return nt;
}

//this version for fitting the 2 electron tracks and the 2 muons all together!
// this will be the first version called -- will also filter out the muons and tracks not used in any good vertex
VertexTracks computeVertices(vector<reco::Track> & coll_1, vector<reco::Track> & coll_2, vector<pat::Muon> & coll_3, vector<pat::Muon> & coll_4, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt) {
    VertexTracks myVertTracks;
    myVertTracks.tracksP = {};
    myVertTracks.tracksN = {};
    myVertTracks.muonsP = {};
    myVertTracks.muonsN = {};
    //vectors of bools to save the order each particle is added to the vector of good particles
    // this will be needed to keep track of what particles are in each vertex!
    vector<int> igood {}; for(size_t i = 0; i < coll_1.size(); i++) igood.push_back(-1);
    vector<int> jgood {}; for(size_t i = 0; i < coll_2.size(); i++) jgood.push_back(-1);
    vector<int> kgood {}; for(size_t i = 0; i < coll_3.size(); i++) kgood.push_back(-1);
    vector<int> lgood {}; for(size_t i = 0; i < coll_4.size(); i++) lgood.push_back(-1);
    //unsigned int printevery = 100000;
    //first loop over positive tracks
    for (size_t i = 0; i < coll_1.size(); i++) {
        reco::Track part_i;
        part_i = coll_1[i];
        //suppose this is an electron
        TLorentzVector ele_i;
        ele_i.SetPtEtaPhiM(part_i.pt(),part_i.eta(),part_i.phi(), ele_mass);
        //now suppose this is a pion
        //TLorentzVector pi_i;
        //pi_i.SetPtEtaPhiM(part_i.pt(),part_i.eta(),part_i.phi(), pi_mass);
        //if ( i > 15 ) break; // ??
        //then loop over negative tracks
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
            //now loop over pos muons
            for(size_t k = 0; k < coll_3.size(); k++) {
                reco::Track part_k;
                part_k = *(coll_3[k].bestTrack());
                //this particle is clearly a muon.
                TLorentzVector mu_k;
                mu_k.SetPtEtaPhiM(part_k.pt(),part_k.eta(),part_k.phi(), mu_mass);
                //mu_k.SetPtEtaPhiM(coll_3[k].pt(),coll_3[k].eta(),coll_3[k].phi(), mu_mass);
                //if (! (part_k.isNonnull()) ) {
                //    continue;
                //}
                //finally loop over neg muons
                for(size_t l = 0; l < coll_4.size(); l++) {
                    reco::Track part_l;
                    part_l = *(coll_4[l].bestTrack());
                    //first make a 4 vector of the 4 particles, see if their invariant mass is in the regime of interest or nah
                    TLorentzVector mu_l;
                    mu_l.SetPtEtaPhiM(part_l.pt(),part_l.eta(),part_l.phi(), mu_mass);
                    //mu_l.SetPtEtaPhiM(coll_4[l].pt(),coll_4[l].eta(),coll_4[l].phi(), mu_mass);

                    //first compute the pT/M of the 4 particle system to make sure it's consistent with eta/eta' meson
                    //first 4-vector assumes tracks are electrons
                    TLorentzVector eta_mmee = ele_i + ele_j + mu_k + mu_l; 
                    ////this 4-vector assumes tracks are pions
                    //TLorentzVector eta_mmpp = pi_i + pi_j + mu_k + mu_l; 
                    float eta_Mmmee = eta_mmee.M();
                    //float eta_Mmmpp = eta_mmpp.M();
                    float eta_Ptmmee = eta_mmee.Pt();
                    //float eta_Ptmmpp = eta_mmpp.Pt();
                    //if the invar. mass of the 4vector isn't in the range of interest, continue
                    bool passed_mmee = false;
                    bool passed_mmpp = false;
                    if( eta_Mmmee > 0.4 && eta_Mmmee < 1.1 && eta_Ptmmee > 12 )  {
                        passed_mmee = true;
                    }
                    //if( eta_Mmmpp > 0.4 && eta_Mmmpp < 1.1 && eta_Ptmmpp > 12 )  {
                    //    passed_mmpp = true;
                    //}
                    if( ! (passed_mmee || passed_mmpp) ) {
                        continue;
                    }

                    float dr = -9999;
                    TransientVertex tv;
                    vector<reco::TransientTrack> transient_tracks{};

                    // build all 4 transient tracks ( 2 el(or pi) + 2 mu )
                    transient_tracks.push_back(theB->build(fix_track(&part_i)));
                    transient_tracks.push_back(theB->build(fix_track(&part_j)));
                    transient_tracks.push_back(theB->build(fix_track(&part_k)));
                    transient_tracks.push_back(theB->build(fix_track(&part_l)));

                    tv = kvf.vertex(transient_tracks);
                    float probVtx = calcVertices(transient_tracks, tv, type, nt);
                    //choose arbitrary very loose cutoff of prob .1
                    if ( probVtx > 0.1 ) {
                    //if ( probVtx > 0 ) {
                    //if ( 0.2 > 0.1 ) {
                        //std::cout << "good vertex! Tracks: " << (int)i << " positive; " << (int)j << " negative; Muons: " << (int)k << " positive; " << (int)l << " negative" << std::endl;
                        //dr of the electron pair
                        dr = reco::deltaR(part_i, part_j);
                        //nt.recoVtxDr_[type].push_back(dr);
                        nt.mmeeVtxDr_.push_back(dr);
                        //save these temporarily, just for debugging!!
                        //nt.recoVtxM_[type].push_back(eta_Mmmee);    /////
                        //nt.recoVtxPt_[type].push_back(eta_Ptmmee);  /////
                        nt.mmeeVtxM_.push_back(eta_Mmmee);    /////
                        nt.mmeeVtxPt_.push_back(eta_Ptmmee);  /////
                        /////////
                        //check if the particle is saved yet or nah before adding it
                            //make sure you also keep track of what particle numbers are involved with each vertex.
                        //if the positive track has not already been added to the list of all good tracks, add it.
                        if( igood[i] < 0 ) {
                            igood[i] = static_cast<int>(nt.recoNGoodTrk_);
                            nt.mmeeTrxP.push_back(nt.recoNGoodTrk_);
                            //std::cout << (int)i << " positive track is the " << (int)nt.recoNGoodTrk_ << " track overall." << std::endl;
                            //std::cout << "OUTSIDE boutta increment recoNGoodTrk: " << (int)nt.recoNGoodTrk_ << std::endl;
                            nt = addTrack(nt, part_i);
                            //std::cout << "OUTSIDE just incremented recoNGoodTrk: " << (int)nt.recoNGoodTrk_ << std::endl;
                            myVertTracks.tracksP.push_back(part_i);
                        }
                        //do the same for negative tracks, etc.
                        if( jgood[j] < 0 ) {
                            jgood[j] = static_cast<int>(nt.recoNGoodTrk_);
                            nt.mmeeTrxN.push_back(nt.recoNGoodTrk_);
                            //std::cout << (int)j << " negative track is the " << (int)nt.recoNGoodTrk_ << " track overall." << std::endl;
                            nt = addTrack(nt, part_j);
                            myVertTracks.tracksN.push_back(part_j);
                        }
                        if ( kgood[k] < 0 ) {
                            kgood[k] = static_cast<int>(nt.recoNGoodMuon_);
                            //std::cout << (int)k << " positive muon is the " << (int)nt.recoNGoodMuon_ << " muon overall." << std::endl;
                            nt.muonsP.push_back(nt.recoNGoodMuon_);
                            nt = addMuon(nt, coll_3[k]);
                            myVertTracks.muonsP.push_back(coll_3[k]);
                        }
                        if( lgood[l] < 0 ) {
                            lgood[l] = static_cast<int>(nt.recoNGoodMuon_);
                            //std::cout << (int)l << " negative muon is the " << (int)nt.recoNGoodMuon_ << " muon overall." << std::endl;
                            nt.muonsN.push_back(nt.recoNGoodMuon_);
                            nt = addMuon(nt, coll_4[l]);
                            myVertTracks.muonsN.push_back(coll_4[l]);
                        } 
                        //for this vertex, two different pointers to constituent particles
                        // first for the tracks (electrons or pions), then for muons
                        //  for each of the two, first 4 bits are for the positive particle, last 4 bits for the neg particle
                        //uint8_t tracks = igood[i] + (jgood[j] << 4); 
                        uint8_t trackP = static_cast<uint8_t>(igood[i]); 
                        uint8_t trackN = static_cast<uint8_t>(jgood[j]); 
                        uint8_t muonP = static_cast<uint8_t>(kgood[k]);
                        uint8_t muonN = static_cast<uint8_t>(lgood[l]);

                        //std::cout << "pushing back tracks " << (int)trackP << ", " << (int)trackN << std::endl;
                        //std::cout << "pushing back muons " << (int)muonP << ", " << (int)muonN << std::endl;
                        //nt.recoVtxTrackP_[type].push_back(trackP);
                        //nt.recoVtxTrackN_[type].push_back(trackN);
                        //nt.recoVtxMuonP_[type].push_back(muonP);
                        //nt.recoVtxMuonN_[type].push_back(muonN);
                        //if(nt.recoVtxTrackP_[type].size() != nt.recoVtxDr_[type].size()){
                        //    std::cout << "Error!!!!!!!!!!" << (int)nt.recoVtxTrackP_[type].size() << " vs. " << (int)nt.recoVtxDr_[type].size() << std::endl;
                        //    throw std::runtime_error("what the fuck");
                        //}
                        ////quadruple checking that the assignment is done correctly
                        //int lastadd = nt.recoVtxMuonN_[type].size()-1;
                        //TLorentzVector trackp2;
                        //int lastTrackp = static_cast<int>(nt.recoVtxTrackP_[type][lastadd]);
                        //trackp2.SetPtEtaPhiM(nt.recoTrkPt_[lastTrackp],nt.recoTrkEta_[lastTrackp],nt.recoTrkPhi_[lastTrackp], ele_mass);
                        //TLorentzVector trackn2;
                        //trackn2.SetPtEtaPhiM(nt.recoTrkPt_[static_cast<int>(nt.recoVtxTrackN_[type][lastadd])],nt.recoTrkEta_[static_cast<int>(nt.recoVtxTrackN_[type][lastadd])],nt.recoTrkPhi_[static_cast<int>(nt.recoVtxTrackN_[type][lastadd])], ele_mass);
                        //TLorentzVector muonp2;
                        //muonp2.SetPtEtaPhiM(nt.recoMuonPt_[static_cast<int>(nt.recoVtxMuonP_[type][lastadd])],nt.recoMuonEta_[static_cast<int>(nt.recoVtxMuonP_[type][lastadd])],nt.recoMuonPhi_[static_cast<int>(nt.recoVtxMuonP_[type][lastadd])], mu_mass);
                        //TLorentzVector muonn2;
                        //muonn2.SetPtEtaPhiM(nt.recoMuonPt_[static_cast<int>(nt.recoVtxMuonN_[type][lastadd])],nt.recoMuonEta_[static_cast<int>(nt.recoVtxMuonN_[type][lastadd])],nt.recoMuonPhi_[static_cast<int>(nt.recoVtxMuonN_[type][lastadd])], mu_mass);
                        //TLorentzVector eta2 = trackp2 + trackn2 + muonp2 + muonn2;
                        ////better be very close to M()...
                        //    //tho not exactly because this is using the Muon info instead of the Muon bestTrack info like before
                        //float mass2 = eta2.M();
                        //nt.recoVtxM2_[type].push_back(mass2);
                        //float pt2 = eta2.Pt();
                        //nt.recoVtxPt2_[type].push_back(pt2);
                        nt.mmeeVtxTrackP_.push_back(trackP);
                        nt.mmeeVtxTrackN_.push_back(trackN);
                        nt.mmeeVtxMuonP_.push_back(muonP);
                        nt.mmeeVtxMuonN_.push_back(muonN);
                        if(nt.mmeeVtxTrackP_.size() != nt.mmeeVtxDr_.size()){
                            std::cout << "Error!!!!!!!!!!" << (int)nt.mmeeVtxTrackP_.size() << " vs. " << (int)nt.mmeeVtxDr_.size() << std::endl;
                            throw std::runtime_error("what the h*ck TrackP size diff from Dr size???");
                        }
                        //quadruple checking that the assignment is done correctly
                        //int lastadd = nt.mmeeVtxMuonN_.size()-1;
                        //TLorentzVector trackp2;
                        //int lastTrackp = static_cast<int>(nt.mmeeVtxTrackP_[lastadd]);
                        //trackp2.SetPtEtaPhiM(nt.recoTrkPt_[lastTrackp],nt.recoTrkEta_[lastTrackp],nt.recoTrkPhi_[lastTrackp], ele_mass);
                        //std::cout << "trackp2: " << lastTrackp << " pT: " << nt.recoTrkPt_[lastTrackp] << std::endl;
                        //TLorentzVector trackn2;
                        //trackn2.SetPtEtaPhiM(nt.recoTrkPt_[static_cast<int>(nt.mmeeVtxTrackN_[lastadd])],nt.recoTrkEta_[static_cast<int>(nt.mmeeVtxTrackN_[lastadd])],nt.recoTrkPhi_[static_cast<int>(nt.mmeeVtxTrackN_[lastadd])], ele_mass);
                        //TLorentzVector muonp2;
                        //muonp2.SetPtEtaPhiM(nt.recoMuonPt_[static_cast<int>(nt.mmeeVtxMuonP_[lastadd])],nt.recoMuonEta_[static_cast<int>(nt.mmeeVtxMuonP_[lastadd])],nt.recoMuonPhi_[static_cast<int>(nt.mmeeVtxMuonP_[lastadd])], mu_mass);
                        //TLorentzVector muonn2;
                        //muonn2.SetPtEtaPhiM(nt.recoMuonPt_[static_cast<int>(nt.mmeeVtxMuonN_[lastadd])],nt.recoMuonEta_[static_cast<int>(nt.mmeeVtxMuonN_[lastadd])],nt.recoMuonPhi_[static_cast<int>(nt.mmeeVtxMuonN_[lastadd])], mu_mass);
                        //TLorentzVector eta2 = trackp2 + trackn2 + muonp2 + muonn2;
                        ////better be very close to M()...
                        //    //tho not exactly because this is using the Muon info instead of the Muon bestTrack info like before
                        //float mass2 = eta2.M();
                        //nt.mmeeVtxM2_.push_back(mass2);
                        //float pt2 = eta2.Pt();
                        //nt.mmeeVtxPt2_.push_back(pt2);
                    }
                } //l loop
            } //k loop

        } // j loop
    } // i loop
    return myVertTracks;
} // computeVertices

//this version for fitting the 2 GsfElectron tracks and the 2 muons all together!
void computeVertices(vector<reco::GsfTrackRef> & coll_1, vector<reco::GsfTrackRef> & coll_2, vector<pat::Muon> & coll_3, vector<pat::Muon> & coll_4, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt) {
    //unsigned int printevery = 100000;
    //first loop over positive tracks
    for (size_t i = 0; i < coll_1.size(); i++) {
        for (size_t j = 0; j < coll_2.size(); j++) {
            //std::cout << "Vertex cand " << (int)i << " " << (int)j << std::endl;
            //reco::TrackRef part_i, part_j;
            //reco::GsfTrackRef part_i, part_j;
            reco::GsfTrackRef part_i, part_j;
            part_i = coll_1[i];
            part_j = coll_2[j];

            TransientVertex tv;
            if (part_i.isNonnull() && part_j.isNonnull() ) { // && i != j) {
                //form 4-vectors for the two electrons
                TLorentzVector ele_i;
                ele_i.SetPtEtaPhiM(part_i->pt(),part_i->eta(),part_i->phi(), ele_mass);
                TLorentzVector ele_j;
                ele_j.SetPtEtaPhiM(part_j->pt(),part_j->eta(),part_j->phi(), ele_mass);
                for(size_t k = 0; k < coll_3.size(); k++) {
                    reco::Track part_k;
                    part_k = *(coll_3[k].bestTrack());
                    //this particle is clearly a muon.
                    TLorentzVector mu_k;
                    mu_k.SetPtEtaPhiM(part_k.pt(),part_k.eta(),part_k.phi(), mu_mass);
                    //finally loop over neg muons
                    for(size_t l = 0; l < coll_4.size(); l++) {
                        reco::Track part_l;
                        part_l = *(coll_4[l].bestTrack());
                        //first make a 4 vector of the 4 particles, see if their invariant mass is in the regime of interest or nah
                        TLorentzVector mu_l;
                        mu_l.SetPtEtaPhiM(part_l.pt(),part_l.eta(),part_l.phi(), mu_mass);

                        //first compute the pT/M of the 4 particle system to make sure it's consistent with eta/eta' meson
                        //first 4-vector assumes tracks are electrons
                        TLorentzVector eta_mmee = ele_i + ele_j + mu_k + mu_l; 
                        float eta_Mmmee = eta_mmee.M();
                        float eta_Ptmmee = eta_mmee.Pt();
                        //if the invar. mass of the 4vector isn't in the range of interest, continue
                        bool passed_mmee = false;
                        bool passed_mmpp = false;
                        if( eta_Mmmee > 0.4 && eta_Mmmee < 1.1 && eta_Ptmmee > 12 )  {
                            passed_mmee = true;
                        }
                        //if( eta_Mmmpp > 0.4 && eta_Mmmpp < 1.1 && eta_Ptmmpp > 12 )  {
                        //    passed_mmpp = true;
                        //}
                        if( ! (passed_mmee || passed_mmpp) ) {
                            continue;
                        }

                        float dr = -9999;
                        TransientVertex tv;
                        vector<reco::TransientTrack> transient_tracks{};

                        // build all 4 transient tracks ( 2 el(or pi) + 2 mu )
                        transient_tracks.push_back(theB->build(part_i));
                        transient_tracks.push_back(theB->build(part_j));
                        transient_tracks.push_back(theB->build(fix_track(&part_k)));
                        transient_tracks.push_back(theB->build(fix_track(&part_l)));

                        tv = kvf.vertex(transient_tracks);
                        float probVtx = calcVertices(transient_tracks, tv, type, nt);
                        //choose arbitrary very loose cutoff of prob .1
                        if ( probVtx > 0.1 ) {
                        //if ( probVtx > 0 ) {
                            //dr of the electron pair
                            dr = reco::deltaR(*part_i, *part_j);
                            //save these temporarily, just for debugging!!
                            //nt.recoVtxM_[type].push_back(eta_Mmmee);    /////
                            //nt.recoVtxPt_[type].push_back(eta_Ptmmee);  /////
                            //nt.mmelelVtxM_.push_back(eta_Mmmee);    /////
                            //nt.mmelelVtxPt_.push_back(eta_Ptmmee);  /////
                            //for this vertex, two different pointers to constituent particles
                            // first for the tracks (electrons or pions), then for muons
                            //  for each of the two, first 4 bits are for the positive particle, last 4 bits for the neg particle
                            //uint8_t tracks = igood[i] + (jgood[j] << 4); 
                            uint8_t muonP = nt.muonsP[k];
                            uint8_t muonN = nt.muonsN[l];

                            if(type == "mmelel") {
                                uint8_t eleP = nt.gsfElsP[i];
                                uint8_t eleN = nt.gsfElsN[j];
                                //nt.recoVtxDr_[type].push_back(dr);
                                nt.mmelelVtxDr_.push_back(dr);
                                nt.mmelelVtxEleP_.push_back(eleP);
                                nt.mmelelVtxEleN_.push_back(eleN);
                                nt.mmelelVtxMuonP_.push_back(muonP);
                                nt.mmelelVtxMuonN_.push_back(muonN);
                                if(nt.mmelelVtxEleP_.size() != nt.mmelelVtxDr_.size()){
                                    std::cout << "Error!!!!!!!!!!" << (int)nt.mmelelVtxEleP_.size() << " vs. " << (int)nt.mmelelVtxDr_.size() << std::endl;
                                    throw std::runtime_error("what the h*ck EleP size diff from Dr size for mmelel???");
                                }
                            } //mmelel type
                            else if(type == "mmlplp") {
                                uint8_t eleP = nt.gsfLowPtElsP[i];
                                uint8_t eleN = nt.gsfLowPtElsN[j];
                                nt.mmlplpVtxDr_.push_back(dr);
                                nt.mmlplpVtxEleP_.push_back(eleP);
                                nt.mmlplpVtxEleN_.push_back(eleN);
                                nt.mmlplpVtxMuonP_.push_back(muonP);
                                nt.mmlplpVtxMuonN_.push_back(muonN);
                                if(nt.mmlplpVtxEleP_.size() != nt.mmlplpVtxDr_.size()){
                                    std::cout << "Error!!!!!!!!!!" << (int)nt.mmlplpVtxEleP_.size() << " vs. " << (int)nt.mmlplpVtxDr_.size() << std::endl;
                                    throw std::runtime_error("what the h*ck EleP size diff from Dr size for mmlplp???");
                                }

                            } //mmlplp type
                            else {
                                std::cout << "Error!!! Unrecognized type " << type << std::endl;
                                throw std::runtime_error("unrecognized type for 4-lepton vertex.");
                            }
                        }
                    } //l loop
                } //k loop
            } //if GsfTracks Nonnull

        } // j loop
    } // i loop
} // computeVertices

//got this code from Sergey Polikarpov
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

