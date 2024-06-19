#python 2 idk where it's at

#this script will take in the sample names from the file "bpgSamples.csv", and set up a different directory for each of them.
#  (by calling makeCondorsam.py for each sample name/location).


import os
import sys

def main():
    
    #submit jobs within this script (T) or just do the setup (F)
    submit_now = False #True

    #true if running the jobs on lxplus instead of cmslpc
    lxplus = True

    #always delete current directory if it already exists?
    always_del = False

    year = 2023 #2022
    
    #number of root files to run in a single job
    nroot = 1 #0 #5

    #do systematics or nah (will take way longer)
    doSyst = False

    #name for parent directory for all the new directories
    parent = "BParking_{}".format(year)
        
    inname = "bpgSamples.csv"
    infile = open(inname, "r")

    for line in infile:
        words = line.split(',')
        
        sample = words[1].strip() 

        #make the new directory
        #nickname of the sample
        samp_name = words[0].strip()
        
        #True if the files will be stored in my personal eos space, false otherwise.
        my_eos = False
        #if "HToAA" in samp_name and "_M-" in samp_name:
        #    my_eos = True
        
        isMC = False #True
        #if "data" in samp_name: isMC = False

        new_name = "%s/%s"%(parent, samp_name)
        new_sample = True
        if os.path.exists(new_name):
            new_sample = False
        else:
            print("new_name: " + new_name)
            os.system("mkdir -p %s"%(new_name) )
        #on lxplus can't really do this.......just make sure the eos directories already exist!!
        #also make a directory on eos for it.
        if not lxplus:
            eos_path = "/eos/uscms/store/user/bgreenbe/BParking_%d/%s"%(year, samp_name)
            #MUST delete all prior contents in the eos directory if it already exists.
            if os.path.exists(eos_path):
                if not always_del:
                    cont = raw_input("Directory %s already exists. Delete all contents? (Y to delete and continue, N to cancel, O to override and proceed without deleting.) "%(eos_path))
                    if cont in ["O", "o", "override", "Override"]: 
                        #allow to continue without deleting (MAKE SURE TO CHANGE THIS BACK!!)
                        always_del = False
                    elif cont not in ["Y", "y", "Yes", "yes"]:
                        sys.exit()
                    else:
                        always_del = True
                if always_del:
                    os.system("rm %s/*.root"%(eos_path))
            os.system("eos root://cmseos.fnal.gov mkdir /store/user/bgreenbe/BParking_%d/%s"%(year, samp_name))
            #os.system("eos root://cmseos.fnal.gov mkdir /store/user/lpcdisptau/eta2mu2e/BParking_%d/%s"%(year, samp_name))
        lang = "tcsh"
        credname = "x509up_u52949"
        if lxplus:
            lang = "bash"
            credname = "x509up_u104084"
        #run the setup code
        os.system("cd %s ; python ../../makeCondorbpg.py --dataSet %s --nickName %s --csv bpgSamples.csv --mode anaXRD --year %d -c %d -p ~/%s -l %s\n"%(new_name, sample, samp_name, year, nroot, credname, lang))
        
        this_dir = os.getcwd() #"/uscms_data/d3/bgreenbe/CMSSW_12_4_13/src/eta2mu2e/eta2mu2eSkimmer/test"
        #come back to get ready for the next one.
        os.system("cd %s"%(this_dir))
        #cp the submit_jobs file to the directory.
        #os.system("cp %s/submit_jobs.py %s"%(parent, new_name))
        os.system("sed \'s/nfiles = 5/nfiles = %d/\' %s/../submit_jobs.py > %s/submit_jobs.py"%(nroot, parent, new_name))
        os.system("cp %s/mv_files.sh %s"%(parent, new_name))
        print("Directory %s created."%(new_name))
        
        ##try to cp fileList.txt to cern eos (if it's a signal sample).
        #if my_eos:
        #    #get mass of the signal sample
        #    mass = samp_name.split('-')[-1]
        #    #no need to cp the file again if it's already there.
        #    if new_sample:
        #        print("command: scp %s/fileList.txt bgreenbe@lxplus.cern.ch:/eos/user/b/bgreenbe/HAA_ntuples/ggha01a01To4tau_%d_%s/"%(new_name, year, mass))
        #        os.system("scp %s/fileList.txt bgreenbe@lxplus.cern.ch:/eos/user/b/bgreenbe/HAA_ntuples/ggha01a01To4tau_%d_%s/"%(new_name, year, mass))
        #        #create the directory in eos for these files to be moved to
        #        os.system("eos root://cmseos.fnal.gov mkdir /store/user/bgreenbe/haa_4tau_%d/signal_M-%s"%(year, mass))
        #    #add signal sample to nAODv7
        #    #os.system("python /uscms/homes/b/bgreenbe/work/CMSSW_10_2_9/src/ZH_Run2/MC/addPUhisto.py -f /uscms/homes/b/bgreenbe/work/CMSSW_10_2_9/src/ZH_Run2/MC/MC_%d_nAODv7.root -ch ZZZ -n HToAATo4Tau_M-%s"%(year, mass))
        #    os.system("python /uscms/homes/b/bgreenbe/work/CMSSW_10_2_9/src/ZH_Run2/MC/addPUhisto.py -f /uscms/homes/b/bgreenbe/work/CMSSW_10_2_9/src/ZH_Run2/MC/MC_%d_nAODv7.root -n HToAATo4Tau_M-%s"%(year, mass)) #ZHToTauTau default.
        #elif new_sample:
#       #     os.system("python /uscms/homes/b/bgreenbe/work/CMSSW_10_2_9/src/ZH_Run2/MC/addPUhisto.py -f /uscms/homes/b/bgreenbe/work/CMSSW_10_2_9/src/ZH_Run2/MC/MC_%d_nAODv7.root -ch ZZZ -n %s"%(year, samp_name)) 
        #    pass

        #submit the jobs now (if we're supposed to)
        if submit_now:
            os.system("cd %s; python submit_jobs.py"%(new_name))
            
if __name__=="__main__":
    main()

