#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <chrono>

std::vector<std::string> readlist(std::string fname) {
    //now read the file into the list
    std::ifstream file(fname); // Open the file
    std::vector<std::string> lines; // Vector to store lines

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) { // Read each line
            lines.push_back(line); // Store the line in the vector
        }
        file.close(); // Close the file
    } 
    else {
        std::cout << "Unable to open file " << fname << std::endl;
        throw std::runtime_error("File cannot be opened error");
    }
    return lines;
} //end readlist function

//find the sub-leading pT of muons in the event e
// return the index of it in the pt array (bc that will also allow us to get the corresponding eta, phi, etc) 
//float find_subleadpt(std::vector<float>* Muon_pt, bool printOn=false) {
uint8_t find_subleadpt(std::vector<float>* Muon_pt, bool printOn=false) {
    float leadpt = 0.0;
    float subleadpt = -1.0;
    if(printOn) {
        std::cout << "inside Muon_pt size = " << (int)Muon_pt->size() << std::endl;
    }
    //now find the sub-leading pT to fill the hist
    for( uint8_t j=0; j< Muon_pt->size(); j++) {
        float pt = (*Muon_pt)[j];
        if(printOn) {
            std::cout << "pt: " << pt << std::endl;
        }
        if(pt > leadpt) {
            subleadpt = leadpt;
            leadpt = pt;
            if(printOn) {
                std::cout << "new leadpt: " << leadpt << "; new subleadpt: " << subleadpt << std::endl;
            }
        }
        else if( pt > subleadpt) {
            subleadpt = pt;
            if(printOn) {
                std::cout << "new subleadpt: " << subleadpt << std::endl;
            }
        }
    } //end j loop
    if(printOn) {
        std::cout << "leadpt: " << leadpt << "; subleadpt: " << subleadpt << std::endl;
    }
    return subleadpt;
} //end find_subleadpt

//void measure_trigCorrex(int argc, char* argv[]) {
void measure_trigCorrex(std::string lett, std::string numb, std::string proc) {
    //start
    auto t_start = std::chrono::high_resolution_clock::now();

    float mu_mass = .105658;

    //max number of files to process (for test), -1 for all
    int nmax = -1;

    bool isMC = false;
    bool isSig = false;
    bool isMuMu = false;
    std::string let = "";
    int num = 0;
    int process = 0;
    //C-style arguments
    //if(argc == 1) {
    //    isMC = true;
    //    isSig = true;
    //    let = "sig";
    //}
    //else if(std::strcmp(argv[1], "bkg") == 0) {
    //    isMC = true;
    //    let = "bkg";
    //}
    //else if(std::strcmp(argv[1], "mumu")) {
    //    isMC = true;
    //    isMuMu = true;
    //    let = "mumu";
    //}
    //else {
    //    stringstream sslet;
    //    sslet << argv[1];
    //    let = sslet.str();
    //    num = std::stoi(argv[2]);
    //}
    num = std::stoi(numb);
    process = std::stoi(proc);
    std::stringstream sslet;
    sslet << lett;
    let = sslet.str();

    std::string DEname;
    std::string Muname;
    if(isMC && isSig){
        DEname = "DEfiles0sig.txt";
        Muname = "sigMCList.txt";
    }
    else if(isMC && !isSig && !isMuMu) {
        DEname = "DEfiles0bkg.txt";
        Muname = "bkgMCList.txt";
    }
    else if(isMC && isMuMu) {
        DEname = "DEfiles0mumu.txt";
        Muname = "mumuMCList.txt";
    }
    else {
        //data
        stringstream ssDEname;
        stringstream ssMuname;
        //ssDEname << "trigfilelists/DEfiles" << num << let << ".txt";
        //ssDEname << "trigfilelists/SMfiles" << num << let << ".txt";
        ssDEname << "trigfilelists/EGfiles" << let << ".txt";
        DEname = ssDEname.str();
        //ssMuname << "Mufiles" << let << ".txt";
        ssMuname << "trigfilelists/flist_" << let << "_" << process << ".txt";
        Muname = ssMuname.str();
        //std::cout << "filenames assigned: DEname " << DEname << ", Muname " << Muname << std::endl;
    }
    std::vector<std::string> DEfiles = readlist(DEname);
    std::vector<std::string> Mufiles = readlist(Muname);

    auto t_readlist = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = t_readlist - t_start;
    //std::cout << "readlist: Elapsed time: " << duration.count() << " seconds" << std::endl;

    if(isMC && DEfiles.size() != Mufiles.size()) {
        throw std::runtime_error("Error: for MC DoubleElectron and DoubleMuon files must exactly correspond.");
    }

    //save all the good events of each run num: map from run number to event numbers
    //std::map<ULong64_t, std::vector<std::pair<ULong64_t, float>>> runs;
    std::map<ULong64_t, std::vector<std::tuple<ULong64_t, float, float, float>>> runs;
    //TH1F* hAll = new TH1F("hAll", "subleading p_{T} of DoubleElectron events with 2 muons", 100, 0.0, 100.0);
    //TH1F* hAll = new TH1F("hAll", "subleading p_{T} of SingleMuon events with 2 muons", 100, 0.0, 100.0);
    //TH1F* hAll = new TH1F("hAll", "subleading p_{T} of EGamma events with 2 muons", 100, 0.0, 100.0);
    TH2F* hsubVdRAllCenter = new TH2F("hsubVdRAllCenter", "subleading p_{T} vs. #Delta R separation of muons (all) w/ |#eta|<1.4", 100, 0.0, 100.0, 100, 0.0, 1.0);
    TH2F* hsubVdRAllOuter = new TH2F("hsubVdRAllOuter", "subleading p_{T} vs. #Delta R separation of muons (all) w/ |#eta|>1.4", 100, 0.0, 100.0, 100, 0.0, 1.0);
    TH2F* hsubVdRTrigCenter = new TH2F("hsubVdRTrigCenter", "subleading p_{T} vs. #Delta R separation of muons (passing trigger) w/ |#eta|<1.4", 100, 0.0, 100.0, 100, 0.0, 1.0);
    TH2F* hsubVdRTrigOuter = new TH2F("hsubVdRTrigOuter", "subleading p_{T} vs. #Delta R separation of muons (passing trigger) w/ |#eta|>1.4", 100, 0.0, 100.0, 100, 0.0, 1.0);
    //TH1F* hPassed = new TH1F("hPassed", "subleading p_{T} of passing DoubleMuon events", 100, 0.0, 100.0);

    //for MC, must compare file by file (whereas for data, compare all files to all files)
    if(!isMC) {
        int nf = 0;
        for(std::string deF : DEfiles) {

            auto t_startdeF = std::chrono::high_resolution_clock::now();
            duration = t_startdeF - t_start;
            //std::cout << "startdeF: Elapsed time: " << duration.count() << " seconds" << std::endl;

            if(nmax > 0 && nf >= nmax) break;
            std::cout << "Trying to open file " << nf << " / " << (int)DEfiles.size() << std::endl;
            std::cout << "Trying to open: " << deF << std::endl;
            stringstream fullpath;
            //fullpath << "root://cmseos.fnal.gov//store/user/bgreenbe/BParking_2022/ParkingDoubleElectronLowMass" << num << let << "/" << deF;
            //fullpath << "root://cmseos.fnal.gov//store/user/bgreenbe/BParking_2022/ParkingSingleMuon" << num << let << "/" << deF;
            fullpath << "root://cmseos.fnal.gov//store/user/bgreenbe/BParking_2022/EGamma" << let << "/" << deF;
            //fullpath << "root://cmseos.fnal.gov//store/user/lpcdisptau/eta2mu2e/BParking_2022/EGamma" << let << "/" << deF;
            //try to xrdcp the file instead of opening it remotely.
            stringstream xrdcmd;
            xrdcmd << "xrdcp " << fullpath.str() << " .";
            std::system(xrdcmd.str().c_str());
            //TFile * elF = TFile::Open(fullpath.str().c_str()); 
            TFile * elF = TFile::Open(deF.c_str()); 

            auto t_openF = std::chrono::high_resolution_clock::now();
            duration = t_openF - t_start;
            //std::cout << "openF: Elapsed time: " << duration.count() << " seconds" << std::endl;

            TTree *elT = (TTree*)elF->Get("ntuples/recoT");
            if(!elT) {
                continue;
            }

            // Variables to hold data from the tree
            std::vector<float> *Muon_pt = nullptr;
            std::vector<uint8_t> *Muon_charge = nullptr;
            std::vector<float> *Muon_phi = nullptr;
            std::vector<float> *Muon_eta = nullptr;
            ULong64_t evt = 0;
            ULong64_t run = 0;

            // Set the branch addresses
            elT->SetBranchAddress("Muon_pt", &Muon_pt);
            elT->SetBranchAddress("Muon_charge", &Muon_charge);
            elT->SetBranchAddress("Muon_phi", &Muon_phi);
            elT->SetBranchAddress("Muon_eta", &Muon_eta);
            elT->SetBranchAddress("evt", &evt);
            elT->SetBranchAddress("run", &run);

            //first go thru the dielectron file and find all the events with at least one pair of OS muons
            int nTot = (int)elT->GetEntries();

            auto t_startEvloop = std::chrono::high_resolution_clock::now();
            duration = t_startEvloop - t_start;
            //std::cout << "startEvloop: Elapsed time: " << duration.count() << " seconds" << std::endl;

            for(int i = 0; i < nTot; i++) {
                auto t_starti = std::chrono::high_resolution_clock::now();
                duration = t_starti - t_start;
                //std::cout << "start i loop: Elapsed time: " << duration.count() << " seconds" << std::endl;

                elT->GetEntry(i);

                auto t_getEntry = std::chrono::high_resolution_clock::now();
                duration = t_getEntry - t_start;
                //std::cout << "GetEntry: Elapsed time: " << duration.count() << " seconds" << std::endl;
                //if(i%10000 == 0) { 
                ////if(i%1 == 0) { 
                //    std::cout << "Event " << i << " / " << nTot << " in DoubleElectron file!" << std::endl;
                //}
                //std::cout << "checking nMuons" << std::endl;
                //requiring EXACTLY 2 muons of opposite charge
                //if(Muon_pt->size() != 2) {
                //    continue;
                //}
                //if((*Muon_charge)[0] == (*Muon_charge)[1]) {
                //    continue;
                //}
                //find the highest pT mu+ and mu-, get their invar mass
                uint8_t idxp, idxn;
                float maxp, maxn;
                maxp = -1.0;
                maxn = -1.0;
                idxp = 0;
                idxn = 0;
                for(uint8_t idx = 0; idx < Muon_pt->size(); idx++) {
                    if((*Muon_charge)[idx] == 1 && (*Muon_pt)[idx] > maxp){
                        maxp = (*Muon_pt)[idx];
                        idxp = idx;
                    }
                    else if((*Muon_charge)[idx] != 1 && (*Muon_pt)[idx] > maxn){
                        maxn = (*Muon_pt)[idx];
                        idxn = idx;
                    }
                } //end idx loop
                //make TLorentzVectors 
                TLorentzVector mu0;
                mu0.SetPtEtaPhiM((*Muon_pt)[idxp], (*Muon_eta)[idxp], (*Muon_phi)[idxp], mu_mass);
                TLorentzVector mu1;
                mu1.SetPtEtaPhiM((*Muon_pt)[idxn], (*Muon_eta)[idxn], (*Muon_phi)[idxn], mu_mass);
                float Mmumu = (mu0+mu1).M();
                ////if(!(Mmumu > 0.45 && Mmumu < 0.65)) {
                //if(Mmumu > 2.0) {
                //    continue;
                //}
                ////if there are 2 opposite-sign muons with good invariant mass, this is a good event!
                ////try all the combos of muons to see if any is good
                ////auto t_muloop = std::chrono::high_resolution_clock::now();
                ////duration = t_muloop - t_start;
                ////std::cout << "Muloop: Elapsed time: " << duration.count() << " seconds" << std::endl;
                //bool goodEv = false;
                //for(uint8_t j=0; j<Muon_pt->size(); j++) {
                //    //std::cout << "j = " << (int)j << std::endl;
                //    for(uint8_t k=j+1; k<Muon_pt->size(); k++) {
                //        //std::cout << "k = " << (int)k << std::endl;
                //        if((*Muon_charge)[j] == (*Muon_charge)[k]) {
                //            continue;
                //        }
                //        //make TLorentzVectors 
                //        TLorentzVector mu0;
                //        mu0.SetPtEtaPhiM((*Muon_pt)[j], (*Muon_eta)[j], (*Muon_phi)[j], mu_mass);
                //        TLorentzVector mu1;
                //        mu1.SetPtEtaPhiM((*Muon_pt)[k], (*Muon_eta)[k], (*Muon_phi)[k], mu_mass);
                //        float Mmumu = (mu0+mu1).M();
                //        if(Mmumu > 0.45 && Mmumu < 0.65) {
                //            goodEv = true;
                //            break; 
                //        }
                //    } //end for k loop
                //    if(goodEv) break;
                //} // end for j loop
                ////auto t_mudone = std::chrono::high_resolution_clock::now();
                ////duration = t_mudone - t_start;
                ////std::cout << "Mudone: Elapsed time: " << duration.count() << " seconds" << std::endl;
                ////if it's not a good event, on to the next one.
                //if(!goodEv) continue;

                //std::cout << "goodEv!" << std::endl;
                auto it = runs.find(run);
                if (it == runs.end()) {
                    //add new entry to the runs map
                    //runs[run] = std::vector<std::pair<ULong64_t, float>>();
                    runs[run] = std::vector<std::tuple<ULong64_t, float, float, float>>();
                } 

                //auto t_find = std::chrono::high_resolution_clock::now();
                //if((t_find - t_start).count() - duration.count() > 0.5) {
                //    duration = t_find - t_start;
                //    std::cout << "find: Elapsed time: " << duration.count() << " seconds" << std::endl;
                //}
                //std::cout << "looking for subleadpt..." << std::endl;
                //float subleadpt = find_subleadpt(Muon_pt); //, true);
                float subleadpt = maxp > maxn ? maxn : maxp;
                float subleadeta = maxp > maxn ? (*Muon_eta)[idxn] : (*Muon_eta)[idxp];

                float dR = mu0.DeltaR(mu1);

                //now add the evt to the runs
                runs[run].push_back(std::make_tuple(evt, subleadpt, dR, subleadeta));
                //auto t_findsub = std::chrono::high_resolution_clock::now();
                //if((t_findsub - t_start).count() - duration.count() > 0.5) {
                //    duration = t_findsub - t_start;
                //    std::cout << "findsub: Elapsed time: " << duration.count() << " seconds" << std::endl;
                //}

                //hAll->Fill(subleadpt);
                if( std::fabs(subleadeta) < 1.4 ) {
                    hsubVdRAllCenter->Fill(subleadpt, dR);
                }
                else {
                    hsubVdRAllOuter->Fill(subleadpt, dR);
                }
                if(subleadpt > 30) {
                    std::cout << "DE subleadpt = " << subleadpt << ": run " << run << " evt " << evt << std::endl;
                    std::cout << "     Muons: pT  \t  eta  \t  phi" << std::endl;
                    std::cout << "         0: " << (*Muon_pt)[0] << "\t" << (*Muon_eta)[0] << "\t" << (*Muon_phi)[0] << std::endl;
                    std::cout << "         1: " << (*Muon_pt)[1] << "\t" << (*Muon_eta)[1] << "\t" << (*Muon_phi)[1] << std::endl;
                }
                auto t_filled = std::chrono::high_resolution_clock::now();
                duration = t_filled - t_start;
                //std::cout << "filled: Elapsed time: " << duration.count() << " seconds" << std::endl;
                //std::cout << "filled!" << std::endl;
            } //end event loop i

            auto t_endEvloop = std::chrono::high_resolution_clock::now();
            duration = t_endEvloop - t_start;
            //std::cout << "endEvloop: Elapsed time: " << duration.count() << " seconds" << std::endl;

            //increment the file counter
            nf++;
            //now delet the file to save space
            stringstream rmcmd;
            rmcmd << "rm " << deF;
            std::system(rmcmd.str().c_str());
        } //end file loop

        nf = 0;
        for(std::string muf : Mufiles) {
            if(nmax > 0 && nf >= nmax) break;
            std::cout << "Starting file " << nf << " / " << (int)Mufiles.size() << std::endl;
            std::cout << "Starting file " << muf << std::endl;
            std::stringstream fullpath;
            fullpath << "root://cmseos.fnal.gov/" << muf;
            //fullpath << "root://cmseos.fnal.gov//store/group/lpcdisptau/eta2mu2e/BParking_2022/ParkingDoubleMuonLowMass" << num << let << "/" << muf;

            //try to xrdcp the file instead of opening it remotely.
            stringstream dmF;
            dmF << "dmfile" << nf << ".root";
            stringstream xrdcmd;
            xrdcmd << "xrdcp " << fullpath.str() << " ./" << dmF.str();
            std::system(xrdcmd.str().c_str());

            //TFile* muF = TFile::Open(fullpath.str().c_str());
            TFile* muF = TFile::Open(dmF.str().c_str());
            TTree *muT = (TTree*)muF->Get("ntuples/recoT");

            // Variables to hold data from the tree
            std::vector<float> *Muon_pt = nullptr;
            std::vector<uint8_t> *Muon_charge = nullptr;
            std::vector<float> *Muon_phi = nullptr;
            std::vector<float> *Muon_eta = nullptr;
            ULong_t Triggers_fired = 0;
            ULong64_t evt = 0;
            ULong64_t run = 0;

            // Set the branch addresses
            muT->SetBranchAddress("Muon_pt", &Muon_pt);
            muT->SetBranchAddress("Muon_charge", &Muon_charge);
            muT->SetBranchAddress("Muon_phi", &Muon_phi);
            muT->SetBranchAddress("Muon_eta", &Muon_eta);
            muT->SetBranchAddress("Triggers_fired0", &Triggers_fired);
            muT->SetBranchAddress("evt", &evt);
            muT->SetBranchAddress("run", &run);

            //then go through the dimuon file and see which of those events is in this file
            int nTot = (int)muT->GetEntries();
            for(int i=0; i < nTot; i++) {
                muT->GetEntry(i);
                //if(i%100000 == 0) { 
                //    std::cout << "Event " << i << " / " << nTot << " in DoubleMuon file!" << std::endl;
                //}
                //if(i > 20000) break; //test!!

                auto it = runs.find(run);
                if (it == runs.end()) {
                    continue;
                }
                //auto itr = std::find(runs[run].begin(), runs[run].end(), evt);
                //if(itr == runs[run].end()) {
                //    continue;
                //}
                unsigned int itr = 0;
                for(itr = 0; itr < runs[run].size(); itr++) {
                    //if(runs[run][itr].first == evt) {
                    if(std::get<0>(runs[run][itr]) == evt) {
                        break;
                    }
                }
                if(itr == runs[run].size()) {
                    continue;
                }

//*****            ////for now, for debugging, don't do this!!!!
            //    //now just make sure it passed our particular triggers
            //    if((Triggers_fired & ((1<<11) + (1<<12) + (1<<17) + (1<<25) + (1<<26) + (1<<28))) == 0){
            //        continue;
            //    }

                //std::cout << "run " << run << " event " << evt << std::endl;
                //std::cout << "outside Muon_pt size = " << (int)Muon_pt->size() << std::endl;
                //    std::cout << "     Muons: pT  \t  eta  \t  phi" << std::endl;
                //    std::cout << "         0: " << (*Muon_pt)[0] << "\t" << (*Muon_eta)[0] << "\t" << (*Muon_phi)[0] << std::endl;
                //    std::cout << "         1: " << (*Muon_pt)[1] << "\t" << (*Muon_eta)[1] << "\t" << (*Muon_phi)[1] << std::endl;
                //this event must be in the other file too! Find the subleadpt
                //float subleadpt = find_subleadpt(Muon_pt, true);
                //float subleadpt = runs[run][itr].second;
                float subleadpt = std::get<1>(runs[run][itr]);
                float dR = std::get<2>(runs[run][itr]);
                float subleadeta = std::get<3>(runs[run][itr]); 

                //hPassed->Fill(subleadpt);
                if(std::fabs(subleadeta) < 1.4) {
                    hsubVdRTrigCenter->Fill(subleadpt, dR);
                }
                else {
                    hsubVdRTrigOuter->Fill(subleadpt, dR);
                }

                //std::cout << "back outside Muon_pt size = " << (int)Muon_pt->size() << std::endl;
                //if(subleadpt > 10) {
                    std::cout << "DM subleadpt = " << subleadpt << ": run " << run << " evt " << evt << std::endl;
                    //std::cout << "     Muons: pT  \t  eta  \t  phi" << std::endl;
                    //std::cout << "         0: " << (*Muon_pt)[0] << "\t" << (*Muon_eta)[0] << "\t" << (*Muon_phi)[0] << std::endl;
                    //std::cout << "         1: " << (*Muon_pt)[1] << "\t" << (*Muon_eta)[1] << "\t" << (*Muon_phi)[1] << std::endl;
                //}
            } //end event i loop
            nf++;
            //delet the file to save space.
            stringstream rmcmd;
            rmcmd << "rm " << dmF.str();
            std::system(rmcmd.str().c_str());
        } //end Mufile loop
    } //end isData block
    else {
        throw runtime_error("MC not implemented yet!");
        ////is MC
        //for nf,deF in enumerate(DEfiles):
        //    if nmax > 0 and nf >= nmax: break
        //    runs = {}
        //    print("Trying to open: %s"%deF) 
        //    elF = ROOT.TFile.Open("root://cmseos.fnal.gov/%s"%(deF)) 
        //    elT = elF.Get("ntuples/recoT")

        //    if not elT:
        //        continue
        //    #first go thru the dielectron file and find all the events with at least one pair of OS muons
        //    nTot = elT.GetEntries()
        //    for i,e in enumerate(elT):
        //        if i%10000 == 0: 
        //            print("Event %d / %d in DoubleElectron file!"%(i, nTot))
        //        if e.Triggers_fired0 == 0:
        //            continue
        //        if ord(e.nGoodMuon) < 2:
        //            continue
        //        #if there are 2 muons, this is a good event!
        //        if e.run not in runs:
        //            runs[e.run] = {}
        //            runs[e.run][0] = []
        //        runs[e.run][0].append(e.evt)

        //        subleadpt = find_subleadpt(e)

        //        hAll.Fill(subleadpt)

        //    muf = Mufiles[nf]
        //    print("Starting file %d / %d"%(nf, len(Mufiles))) 
        //    muF = ROOT.TFile.Open("root://cmseos.fnal.gov/%s"%(muf)) 
        //    muT = muF.Get("ntuples/recoT")

        //    #then go through the dimuon file and see which of those events is in this file
        //    nTot = muT.GetEntries()
        //    for i,e in enumerate(muT):
        //        if i%100000 == 0: 
        //            print("Event %d / %d in DoubleMuon file!"%(i, nTot))
        //        #if i > 20000: break #test!!
        //        if e.Triggers_fired0 & ((1<<11) + (1<<12) + (1<<17) + (1<<25) + (1<<26) + (1<<28)) == 0:
        //            continue
        //        if e.run not in runs:
        //            continue
        //        if e.evt not in runs[e.run][0]:
        //            continue

        //        #this event must be in the other file too! Find the subleadpt
        //        subleadpt = find_subleadpt(e)

        //        hPassed.Fill(subleadpt)
    } //end isMC block

    //now open an output file
    std::stringstream foutname;
    foutname << "root://cmseos.fnal.gov//store/user/bgreenbe/BParking2022/trigEff/trigEff" << num << let << "_" << process << "_2022.root";
    std::cout << "Opening output file " << foutname.str() << std::endl;
    TFile *fout = TFile::Open(foutname.str().c_str(), "recreate");
    //hAll->Write();
    hsubVdRAllCenter->Write();
    hsubVdRAllOuter->Write();
    //hPassed->Write();
    hsubVdRTrigCenter->Write();
    hsubVdRTrigOuter->Write();
    fout->Close();
} //end main method
