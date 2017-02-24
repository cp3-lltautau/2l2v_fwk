# Installation for 76X (2015)
```bash 
export SCRAM_ARCH=slc6_amd64_gcc493
scramv1 project CMSSW CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src/
cmsenv
wget -O - --no-check-certificate https://raw.githubusercontent.com/cms2l2v/2l2v_fwk/master/TAGS.txt | sh
```

Now, ONLY after the code has finished to compile insert inside UserCode/llvv_fwk/BuildFile.xml
``` c++
<use name="ZZMatrixElement/MELA"/>
scram b -j 12  
```

#Step to use MELA
```bash 
cd CMSS_X_Y_Z/src
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
cd ZZMatrixElement
sh setup.sh -j 12
```

