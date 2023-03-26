Code for skimming MINIAOD files for use in CMS Run3 2022 BParking SM eta meson -> mu mu e e decay search.

Authors: Bennett Greenberg, Andre Frankenthal

Email: bennett.greenberg@cern.ch

# Set up

```bash
cmsrel CMSSW_12_4_12
cd CMSSW_12_4_12/src
cmsenv
git cms-init
```
Now clone this repository, either via SSH keys:

```bash
git clone git@github.com:bennettgberg/eta2mu2eSkimmer.git eta2mu2e/eta2mu2eSkimmer
```

or via HTTPS if you don't have SSH keys set up:

```bash
git clone https://github.com/bennettgberg/eta2mu2eSkimmer.git eta2mu2e/eta2mu2eSkimmer
```

Then compile everything:

```bash
scram b
```

Finally, run a CMSSW configuration:

```bash
cmsRun eta2mu2e/eta2mu2eSkimmer/test/run_ntuplizer_cfg.py
```
