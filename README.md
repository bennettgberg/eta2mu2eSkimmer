Code for skimming MINIAOD files for use in CMS Run3 2022 BParking SM eta meson -> mu mu e e decay search.

Authors: Bennett Greenberg, Andre Frankenthal

Email: bennett.greenberg@cern.ch

# Set up

```bash
cmsrel CMSSW_12_4_12
cd CMSSW_12_4_12/src
cmsenv
git cms-init
git clone https://github.com/bennettgberg/eta2mu2eSkimmer.git eta2mu2e/eta2mu2eSkimmer
scram b
cmsRun eta2mu2e/eta2mu2eSkimmer/test/run_ntuplizer_cfg.py
```
