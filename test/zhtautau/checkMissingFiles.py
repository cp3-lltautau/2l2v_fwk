#!/usr/bin/env python
import os,sys,time
import json
import optparse
import commands
import subprocess, shlex

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

(opt, args) = parser.parse_args()

#open the file which describes the sample
jsonFile = open(opt.samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

for procBlock in procList :
    #run over processes
    for proc in procBlock[1] :
        data = proc['data']
	#run over samples
        for procData in data : 
            origdtag = getByLabel(procData,'dtag','')
            if(origdtag=='') : continue
            dtag = origdtag        
	    sample =  dtag+'_[0-9]*[0-9]'

	    command_line  = "ls "+opt.theFolder+"/*.py   | grep -c "+sample+" "
	    command_line2 = "ls "+opt.theFolder+"/*.root | grep -c "+sample+" "
	    # args = shlex.split(command_line)
	    # print command_line
	    try: 
	    	jobs  = int( subprocess.check_output(command_line  , shell=True) )
	    	process = subprocess.Popen(command_line2 , stdout=subprocess.PIPE, shell=True) 
		returncode = process.wait()
		# print('ping returned {0}'.format(returncode))
		files = int( process.stdout.read() )
	    except subprocess.CalledProcessError as e:
	        raise RuntimeError("command '{0}' return with error (code {1}): {2}".format(e.cmd, e.returncode, e.output))
	    # print "{0} - {1}".format(files,jobs) 
	    # jobs  = int( os.system('sample='+dtag+'_[0-9]*[0-9]; echo `ls '+opt.theFolder+'/*.py | grep -c "$sample"`') )
	    # files = int( os.system('sample='+dtag+'_[0-9]*[0-9]; echo `ls '+opt.theFolder+'/*.root | grep -c "$sample"`') ) 
	    if jobs != files:
	    	print bcolors.OKGREEN + dtag + bcolors.ENDC
		print 'Missing Files: {0} - {1}'.format(files,jobs)
	   # elif jobs == files : 
	   #	print 'Dataset processed!!'
	
