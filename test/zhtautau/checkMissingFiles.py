#!/usr/bin/env python
import os,sys,time
import json
import optparse
import commands
import subprocess, shlex
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def getByLabel(proc,key,defaultVal=None) :
    try :
        return proc[key]
    except KeyError:
        return defaultVal

def getByLabelFromKeyword(proc,keyword,key,defaultVal=None) :
    try:
        print "found keyword option %s" % str(proc[keyword][key])
        return proc[keyword][key]
    except:
        try :
            return proc[key]
        except KeyError:
            return defaultVal

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-f', '--folder'     ,    dest='theFolder'      , help='folder containing root files'  , default=''             )
parser.add_option('-j', '--json'       ,    dest='samplesDB'      , help='samples json file'             , default='./sample.json')
parser.add_option('-t', '--tag'        ,    dest='onlytag'            , help='process only samples matching this tag'    , default='all')
parser.add_option('-k', '--key'        ,    dest='onlykeyword'        , help='process only samples matching this keyword', default='')


(opt, args) = parser.parse_args()

#open the file which describes the sample
jsonFile = open(opt.samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

#lists for plotting
jobsNumber    = []
pendingNumber = []
runningNumber = []
finishNumber  = []
missingNumber = []
names         = []
position      = []
pos = 0

for procBlock in procList :
    #run over processes
    processedDataset = []
    noProcessedDataset = []
    for proc in procBlock[1] :
        data = proc['data']
        keywords = getByLabel(proc,'keys',[])
        #print "{0} and {1}".format(len(keywords), opt.onlykeyword )
        if(opt.onlykeyword!='' and len(keywords)>0 and opt.onlykeyword not in keywords): continue
#        print keywords
        if(getByLabel(proc,'nosample'      , '')!=''):continue
	#run over samples

        for procData in data :
            origdtag = getByLabel(procData,'dtag','')
            if(origdtag=='') : continue
            dtag = origdtag
            if opt.onlytag!='all' and dtag.find(opt.onlytag)<0 : continue
            sample =  dtag+'_[0-9]*[0-9]'

            userName = str(subprocess.check_output("whoami",shell=True))

            squeueCMD     = "squeue -u ccaputo --format=\"%.18i %.5P %.184j %.2t %.10M\" | cut -d\'/\' -f17 | grep -c "+dtag+" "
            squeueCMD_R   = "squeue -u ccaputo --format=\"%.18i %.5P %.184j %.2t %.10M %T\" | cut -d\'/\' -f17 | grep "+dtag+" | grep -c RUNNING"
            squeueCMD_PD  = "squeue -u ccaputo --format=\"%.18i %.5P %.184j %.2t %.10M %T\" | cut -d\'/\' -f17 | grep "+dtag+" | grep -c PENDING"
            command_line  = "ls "+opt.theFolder+"/*.py   | grep -c "+sample+" "
            command_line2 = "ls "+opt.theFolder+"/*.root | grep -c "+sample+" "
            # args = shlex.split(command_line)
            # print command_line
            try:
                njobsCommand     = subprocess.Popen(command_line  , stdout=subprocess.PIPE, shell=True)
                processSqueue    = subprocess.Popen(squeueCMD, stdout=subprocess.PIPE, shell=True)
                processSqueue_R  = subprocess.Popen(squeueCMD_R, stdout=subprocess.PIPE, shell=True)
                processSqueue_PD = subprocess.Popen(squeueCMD_PD, stdout=subprocess.PIPE, shell=True)
                process          = subprocess.Popen(command_line2 , stdout=subprocess.PIPE, shell=True)
                returncode       = process.wait()
                # print('ping returned {0}'.format(returncode))
                jobs    = int( njobsCommand.stdout.read() )
                onFarm  = int( processSqueue.stdout.read() )
                running = int( processSqueue_R.stdout.read() )
                pending = int( processSqueue_PD.stdout.read() )
                files = int( process.stdout.read() )

                #Filling lists
                jobsNumber.append(jobs)
                pendingNumber.append(pending)
                runningNumber.append(running)
                finishNumber.append(files)
                missingNumber.append(jobs-files)
                names.append(dtag)
                position.append(pos)
                pos = pos + 1
            except subprocess.CalledProcessError as e:
                raise RuntimeError("command '{0}' return with error (code {1}): {2}".format(e.cmd, e.returncode, e.output))
                # print "{0} - {1}".format(files,jobs)
                # jobs  = int( os.system('sample='+dtag+'_[0-9]*[0-9]; echo `ls '+opt.theFolder+'/*.py | grep -c "$sample"`') )
                # files = int( os.system('sample='+dtag+'_[0-9]*[0-9]; echo `ls '+opt.theFolder+'/*.root | grep -c "$sample"`') )
            if jobs!= 0 and jobs != files:
                print bcolors.OKGREEN + dtag + bcolors.ENDC
                print 'Files: {0} of {1}'.format(files,jobs)
                print bcolors.BOLD + "  Running(Pending) / FileMissing: "+bcolors.ENDC+" {0}({1}) / {2}".format(running,pending,int(jobs - files) )
            elif jobs!= 0 and files!=0 and jobs == files :
                processedDataset.append(dtag)
            elif jobs == 0:
                noProcessedDataset.append(dtag)

print bcolors.OKBLUE + 'Dataset processed' + bcolors.ENDC
for item in processedDataset:
    print "- "+item

print bcolors.WARNING + 'Dataset NOT processed' + bcolors.ENDC
for item in noProcessedDataset:
    print "- "+item

jobsNumberArray    = np.array(jobsNumber)
pendingNumberArray = np.array(pendingNumber)
runningNumberArray = np.array(runningNumber)
finishNumberArray  = np.array(finishNumber)
missingNumberArray = np.array(missingNumber)
namesArray         = np.array(names)
positionArray      = np.array(position)

# barWidth = 1
#
# # Create brown bars
# plt.bar(position, finishNumber, color='#7f6d5f', edgecolor='white', width=barWidth)
# # Create green bars (middle), on top of the firs ones
# plt.bar(position, missingNumber, bottom=finishNumber, color='#557f2d', edgecolor='white', width=barWidth)
#
# # Custom X axis
# plt.xticks(position, names, fontweight='bold')
# plt.xlabel("group")
#
# # Show graphic
# plt.savefig('job_monitoring.png', dpi=400)
