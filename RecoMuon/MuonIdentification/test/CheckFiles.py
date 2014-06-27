import ROOT
import os
import sys

FileList = os.popen("cmsLs  /store/group/upgrade/muon/ME0GlobalReco/PU_ParallelRun/ | grep Matched | gawk '{print $5}'").readlines()
#print FilesList
for line in FileList:
    #print line.replace('\n','')
    #print "edmFileUtil root://eoscms//eos/cms" + line.replace('\n','')
    FileInfo = os.popen("edmFileUtil root://eoscms//eos/cms" + line.replace('\n','')).readlines()
    Unclosed = 0
    Events = -1
    for infoline in FileInfo:
        if "Could not open" in infoline:
            Unclosed = 1
        if "events" in infoline:
            print infoline.split("events")[0].split(",")[-1]+" events in "+line
    if Unclosed==1:
        print line+ " is unclosed"
    

