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
    if ( prob > 0.1 ) {
        nt.recoVtxReducedChi2_[type].push_back(vtx_chi2);
        nt.recoVtxVxy_[type].push_back(vxy);
        nt.recoVtxVz_[type].push_back(vz);
        nt.recoVtxSigmaVxy_[type].push_back(sigma_vxy);
        //std::cout << "probability: " << prob << std::endl;
    }
    return prob;
    
} 

//compute vertices for two track vectors coll_1 and coll_2, and add them to the ntuple.
// (overloaded below)
// Return type: struct VertexTracks (defined above)
VertexTracks computeVertices(vector<reco::Track> & coll_1, vector<reco::Track> & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt) {
    //initialize the struct that will be returned
    VertexTracks myVertTracks;
    //initialize the tracks list in the struct that will be returned.
    myVertTracks.tracksP = {};
    myVertTracks.tracksN = {};
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
                dr = reco::deltaR(part_i, part_j);
                nt.recoVtxDr_[type].push_back(dr);
                //there are positive AND negative tracks
                // number pushed back is 4 bits for the first track # (in the good list), followed by 4 bits for the 2nd track #
                //  so max we can handle is 16 good tracks...
                //uint8_t full_val = nt.mmeeTrxP[i] + (nt.mmeeTrxN[j] << 4);
                uint8_t trackP = nt.mmeeTrxP[i];
                uint8_t trackN = nt.mmeeTrxN[j];
                //nt.recoVtxTracks_[type].push_back(full_val);
                nt.recoVtxTrackP_[type].push_back(trackP);
                nt.recoVtxTrackN_[type].push_back(trackN);
                myVertTracks.tracksP.push_back(part_i);
                myVertTracks.tracksN.push_back(part_j);
            }
            
        } // j loop
    } // i loop
    return myVertTracks;
} // computeVertices

//overloaded function
// now compute vertex for two sets of GsfElectron vectors coll_1 and coll_2.
VertexTracks computeVertices(vector<reco::GsfTrackRef> & coll_1, vector<reco::GsfTrackRef> & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf, NtupleContainer & nt) {
    //we're gonna keep all Gsf electrons anyway (there aren't that many of them) so this struct will just be empty.
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
                    dr = reco::deltaR(*part_i, *part_j);
                    nt.recoVtxDr_[type].push_back(dr);
                    //first 4 bits for first electron, last 4 for second electron
                    //uint8_t full_val = nt.gsfElsP[i] + (nt.gsfElsN[j] << 4);
                    uint8_t eleP = nt.gsfElsP[i];
                    uint8_t eleN = nt.gsfElsN[j];
                    nt.recoVtxEleP_[type].push_back(eleP);
                    nt.recoVtxEleN_[type].push_back(eleN);
                }
            } 

        } // j loop
    } // i loop
    return myVertTracks;
} // computeVertices

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
    uint8_t nmuons = 0;
    uint8_t ntracks = 0;
    //unsigned int printevery = 100000;
    //first loop over positive tracks
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

                    //first compute the pT/M of the 4 particle system to make sure it's consistent with eta/eta' meson
                    //first 4-vector assumes tracks are electrons
                    TLorentzVector eta_mmee = ele_i + ele_j + mu_k + mu_l; 
                    //this 4-vector assumes tracks are pions
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

                    // build all 4 transient tracks ( 2 el(or pi) + 2 mu )
                    transient_tracks.push_back(theB->build(fix_track(&part_i)));
                    transient_tracks.push_back(theB->build(fix_track(&part_j)));
                    transient_tracks.push_back(theB->build(fix_track(&part_k)));
                    transient_tracks.push_back(theB->build(fix_track(&part_l)));

                    tv = kvf.vertex(transient_tracks);
                    float probVtx = calcVertices(transient_tracks, tv, type, nt);
                    //if ( probVtx > 0 ) {
                    if ( probVtx > 0.1 ) {
                        //dr of the electron pair
                        dr = reco::deltaR(part_i, part_j);
                        nt.recoVtxDr_[type].push_back(dr);
                        //check if the particle is saved yet or nah before adding it
                            //make sure you also keep track of what particle numbers are involved with each vertex.
                        //if the positive track has not already been added to the list of all good tracks, add it.
                        if( igood[i] < 0 ) {
                            myVertTracks.tracksP.push_back(part_i);
                            nt.mmeeTrxP.push_back(ntracks);
                            igood[i] = ntracks++;
                        }
                        //do the same for negative tracks, etc.
                        if( jgood[j] < 0 ) {
                            myVertTracks.tracksN.push_back(part_j);
                            nt.mmeeTrxN.push_back(ntracks);
                            jgood[j] = ntracks++;
                        }
                        if ( kgood[k] < 0 ) {
                            myVertTracks.muonsP.push_back(coll_3[k]);
                            kgood[k] = nmuons++;
                        }
                        if( lgood[l] < 0 ) {
                            myVertTracks.muonsN.push_back(coll_4[l]);
                            lgood[l] = nmuons++;
                        } 
                        //for this vertex, two different pointers to constituent particles
                        // first for the tracks (electrons or pions), then for muons
                        //  for each of the two, first 4 bits are for the positive particle, last 4 bits for the neg particle
                        //uint8_t tracks = igood[i] + (jgood[j] << 4); 
                        uint8_t trackP = igood[i]; 
                        uint8_t trackN = jgood[j]; 
                        //uint8_t muons = kgood[k] + (lgood[l] << 4);
                        uint8_t muonP = kgood[k];
                        uint8_t muonN = lgood[l];
                        //only 4 bits to represent the particle numbers, so can't handle more than 16 particles!!
                        //if ( igood[i] > 15 || jgood[j] > 15 ) {
                        //    throw std::runtime_error("Too many good tracks!");
                        //}
                        //if ( kgood[k] > 15 || lgood[l] > 15 ) {
                        //    throw std::runtime_error("Too many good muons!");
                        //}
                        //nt.recoVtxTracks_[type].push_back(tracks);
                        nt.recoVtxTrackP_[type].push_back(trackP);
                        nt.recoVtxTrackN_[type].push_back(trackN);
                        //nt.recoVtxMuons_[type].push_back(muons);
                        nt.recoVtxMuonP_[type].push_back(muonP);
                        nt.recoVtxMuonN_[type].push_back(muonN);
                    }
                } //l loop
            } //k loop

        } // j loop
    } // i loop
    return myVertTracks;
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

