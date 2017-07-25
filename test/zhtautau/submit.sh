#!/usr/bin/env bash

#--------------------------------------------------
# Global Code
#--------------------------------------------------

if [[ $# -eq 0 ]]; then
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]"
    printf "\nOPTIONS\n"
    printf "\n\t%-5s  %-40s\n"  "0"  "completely clean up the directory"
    printf "\n\t%-5s  %-40s\n"  "1"  "run 'runZHTauTauAnalysis' on samples.json"
    printf "\n\t%-5s  %-40s\n"  "1.1"  "run 'runZHTauTauAnalysis' on photon_samples.json"
    printf "\n\t%-5s  %-40s\n"  "2"  "compute integrated luminosity from processed samples"
    printf "\n\t%-5s  %-40s\n"  "2.1"  "compute integrated luminosity from processed samples connecting to lxplus (ssh)"
    printf "\n\t%-5s  %-40s\n"  "3"  "make plots and combine root files"
    printf "\n\t%-5s  %-40s\n"  "3.1"  "make plots for photon_smaples"
fi

step=$1   #variable that store the analysis step to run

#Additional arguments to take into account
arguments=''; for var in "${@:2}"; do arguments=$arguments" "$var; done
if [[ $# -ge 4 ]]; then echo "Additional arguments will be considered: "$arguments ;fi

#--------------------------------------------------
# Global Variables
#--------------------------------------------------
SUFFIX=_2017_07_20
#SUFFIX=$(date +"_%Y_%m_%d")
MAINDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau
JSON=$MAINDIR/samples.json
RESULTSDIR=$MAINDIR/results$SUFFIX
PLOTSDIR=$MAINDIR/plots${SUFFIX}
PLOTTER=$MAINDIR/plotter${SUFFIX}

#printf "Result dir is set as: \n\t%s\n" "$RESULTSDIR"


case $step in

    0)  #analysis cleanup
	echo "ALL DATA WILL BE LOST! [N/y]?"
	read answer
	if [[ $answer == "y" ]];
	then
	    echo "CLEANING UP..."
	    rm -rdf $RESULTSDIR $PLOTSDIR LSFJOB_* core.* *.sh.e* *.sh.o*
	fi
	;;

    1)  #submit jobs for 2l2v analysis
	echo "JOB SUBMISSION"
	echo "Input: " $JSON
	echo "Output: " $RESULTSDIR

	queue='8nh'
	#IF CRAB3 is provided in argument, use crab submissiong instead of condor/lsf
	if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi
	runAnalysisOverSamples.py -e runZHTauTauAnalysis -j $JSON -o $RESULTSDIR  -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @data_pileup=datapileup_latest" -s $queue --report True $arguments
	;;

    1.1)  #submit jobs for 2l2v photon jet analysis
	echo "JOB SUBMISSION for Photon + Jet analysis"
	queue='8nh'
	JSON=$MAINDIR/photon_samples.json
	echo "Input: " $JSON
	echo "Output: " $RESULTSDIR
	if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi
	runAnalysisOverSamples.py -e runZHTauTauAnalysis -j $JSON -o $RESULTSDIR -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True $arguments
	;;

    2)  #extract integrated luminosity of the processed lumi blocks
	echo "MISSING LUMI WILL APPEAR AS DIFFERENCE LUMI ONLY IN in.json"
	mergeJSON.py --output=$RESULTSDIR/json_all.json        $RESULTSDIR/Data*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleMu.json   $RESULTSDIR/Data*_DoubleMu*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleEl.json   $RESULTSDIR/Data*_DoubleElectron*.json
	mergeJSON.py --output=$RESULTSDIR/json_muEG.json   $RESULTSDIR/Data*_MuEG*.json
	mergeJSON.py --output=$RESULTSDIR/json_in.json  Cert_*.txt
	echo "MISSING LUMI BLOCKS IN DOUBLE MU DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleMu.json
	echo "MISSING LUMI BLOCKS IN DOUBLE ELECTRON DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleEl.json
	echo "MISSING LUMI BLOCKS IN MUON EGAMMA DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_muEG.json

	echo "COMPUTE INTEGRATED LUMINOSITY"
	export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
	pip install --install-option="--prefix=$HOME/.local" brilws &> /dev/null #will be installed only the first time
	brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i $RESULTSDIR/json_all.json -u /pb -o $RESULTSDIR/LUMI.txt
	tail -n 3 $RESULTSDIR/LUMI.txt
	;;
    2.1)  #extract integrated luminosity of the processed lumi blocks
	CACHE_ALL=$RESULTSDIR/js_all_cache.txt
	CACHE_DMU=$RESULTSDIR/js_doubleMu_cache.txt
	CACHE_SMU=$RESULTSDIR/js_singleMu_cache.txt
	CACHE_DELE=$RESULTSDIR/js_doubleEl_cache.txt
	CACHE_SELE=$RESULTSDIR/js_singleEl_cache.txt
	if [ ! -e $CACHE_ALL ] || [ ! -e $CACHE_DMU ] || [ ! -e $CACHE_DELE ] || [ ! -e $CACHE_SMU ] || [ ! -e $CACHE_SELE ]
		then
			echo "MISSING LUMI WILL APPEAR AS DIFFERENCE LUMI ONLY IN in.json"
			ls -l $RESULTSDIR/Data*.json                 > $CACHE_ALL
			ls -l $RESULTSDIR/Data*_DoubleMu*.json       > $CACHE_DMU
			ls -l $RESULTSDIR/Data*_SingleMu*.json       > $CACHE_SMU
			ls -l $RESULTSDIR/Data*_DoubleElectron*.json > $CACHE_DELE
			ls -l $RESULTSDIR/Data*_SingleElectron*.json > $CACHE_SELE
			mergeJSON.py --output=$RESULTSDIR/json_all.json        $RESULTSDIR/Data*.json
			mergeJSON.py --output=$RESULTSDIR/json_doubleMu.json   $RESULTSDIR/Data*_DoubleMu*.json
			mergeJSON.py --output=$RESULTSDIR/json_singleMu.json   $RESULTSDIR/Data*_SingleMu*.json
			mergeJSON.py --output=$RESULTSDIR/json_doubleEl.json   $RESULTSDIR/Data*_DoubleElectron*.json
			mergeJSON.py --output=$RESULTSDIR/json_singleEl.json   $RESULTSDIR/Data*_SingleElectron*.json
			mergeJSON.py --output=$RESULTSDIR/json_in.json  Cert_*Collisions16*.txt
			echo "MISSING LUMI BLOCKS IN DOUBLE MU DATASET"
			compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleMu.json
			echo "MISSING LUMI BLOCKS IN SINGLE MU DATASET"
			compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_singleMu.json
			echo "MISSING LUMI BLOCKS IN DOUBLE ELECTRON DATASET"
			compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleEl.json
			echo "MISSING LUMI BLOCKS IN SINGLE ELECTRON DATASET"
			compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_singleEl.json
		else
			ls -l $RESULTSDIR/Data*.json                 > ${CACHE_ALL}_tmp
        		ls -l $RESULTSDIR/Data*_DoubleMu*.json       > ${CACHE_DMU}_tmp
			ls -l $RESULTSDIR/Data*_SingleMu*.json       > ${CACHE_SMU}_tmp
			ls -l $RESULTSDIR/Data*_DoubleElectron*.json > ${CACHE_DELE}_tmp
			ls -l $RESULTSDIR/Data*_SingleElectron*.json > ${CACHE_SELE}_tmp
			if [ "`md5sum $CACHE_ALL | cut -d ' ' -f1`" == "`md5sum ${CACHE_ALL}_tmp | cut -d ' ' -f1`" ] && \
                           [ "`md5sum $CACHE_DMU | cut -d ' ' -f1`" == "`md5sum ${CACHE_DMU}_tmp | cut -d ' ' -f1`" ] && \
                           [ "`md5sum $CACHE_SMU | cut -d ' ' -f1`" == "`md5sum ${CACHE_SMU}_tmp | cut -d ' ' -f1`" ] && \
                           [ "`md5sum $CACHE_DELE | cut -d ' ' -f1`" == "`md5sum ${CACHE_DELE}_tmp | cut -d ' ' -f1`" ] && \
                           [ "`md5sum $CACHE_SELE | cut -d ' ' -f1`" == "`md5sum ${CACHE_SELE}_tmp | cut -d ' ' -f1`" ]
				then
					echo "JSON files aren't changed. Please use "$RESULTSDIR"/LUMI.txt"
					exit
				else
					echo "JSON files are changed!!"
					mv ${CACHE_ALL}_tmp $CACHE_ALL
					mv ${CACHE_DMU}_tmp $CACHE_DMU
					mv ${CACHE_SMU}_tmp $CACHE_SMU
					mv ${CACHE_DELE}_tmp $CACHE_DELE
					mv ${CACHE_SELE}_tmp $CACHE_SELE
					mergeJSON.py --output=$RESULTSDIR/json_all.json        $RESULTSDIR/Data*.json
					mergeJSON.py --output=$RESULTSDIR/json_doubleMu.json   $RESULTSDIR/Data*_DoubleMu*.json
					mergeJSON.py --output=$RESULTSDIR/json_singleMu.json   $RESULTSDIR/Data*_SingleMu*.json
					mergeJSON.py --output=$RESULTSDIR/json_doubleEl.json   $RESULTSDIR/Data*_DoubleElectron*.json
					mergeJSON.py --output=$RESULTSDIR/json_singleEl.json   $RESULTSDIR/Data*_SingleElectron*.json
					mergeJSON.py --output=$RESULTSDIR/json_in.json  Cert_*Collisions16*.txt
					echo "MISSING LUMI BLOCKS IN DOUBLE MU DATASET"
					compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleMu.json
					echo "MISSING LUMI BLOCKS IN SINGLE MU DATASET"
					compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_singleMu.json
					echo "MISSING LUMI BLOCKS IN DOUBLE ELECTRON DATASET"
					compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleEl.json
					echo "MISSING LUMI BLOCKS IN SINGLE ELECTRON DATASET"
					compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_singleEl.json
					mergeJSON.py --output=$RESULTSDIR/json_all.json        $RESULTSDIR/Data*.json
			fi
	fi
	echo "COMPUTE INTEGRATED LUMINOSITY"
	echo "Coping json file to lxplus area...."
	LXPLUS_USERNAME=""
	while [ "$LXPLUS_USERNAME" == "" ]; do
		echo -n "Enter your lxplus username --> "
		read LXPLUS_USERNAME
	done
	scp $RESULTSDIR/json_all.json $LXPLUS_USERNAME@lxplus.cern.ch:~/private/
	RESULT=$?
	if [[ $RESULT -eq 0 ]] ; then
       		echo "Copied!"
		ssh $LXPLUS_USERNAME@lxplus.cern.ch 'export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH;
					    cd /afs/cern.ch/user/'${LXPLUS_USERNAME:0:1}'/'$LXPLUS_USERNAME'/private;
					    brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i ./json_all.json -u /pb -o ./LUMI.txt;
					    grep "Summary" -A 2 ./LUMI.txt' > $RESULTSDIR/LUMI.txt
		echo -n "LUMI.txt copied in "$RESULTSDIR
	else
		echo "Impossible to copy json_all.json on your lxplus private area!"
		exit
	fi

	;;
    3)  # make plots and combined root files
        if [ -f $RESULTSDIR/LUMI.txt ]; then
           INTLUMI=`tail -n 1 $RESULTSDIR/LUMI.txt | cut -d ',' -f 6`
        else
           INTLUMI=35866.932 #correspond to the value from DoubleMu OR DoubleEl OR MuEG without jobs failling and golden JSON
           echo "WARNING: $RESULTSDIR/LUMI.txt file is missing so use fixed integrated luminosity value, this might be different than the dataset you ran on"
        fi
	echo "MAKE PLOTS AND SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI"
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/ --outFile $PLOTTER.root  --json $JSON --no2D $arguments
	ln -s -f $PLOTTER.root $MAINDIR/plotter.root
	;;

    4)  # make plots and combined root files
        extractFRWeights --inFile $PLOTTER.root --outFile FR_Weights.root --outDir $PLOTSDIR/FakeRateWeights

esac
