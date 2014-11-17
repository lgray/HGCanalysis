#!/usr/bin/env python

import os,sys
import optparse
import commands
import time

cmsswBase=os.environ['CMSSW_BASE']
cmsswVersion=os.environ['CMSSW_VERSION']

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                                        , default='local')
parser.add_option('-t', '--tag'        ,    dest='tag'                , help='tag'                                                , default='Single13_%s_v2'%cmsswVersion)
parser.add_option('-s', '--step'       ,    dest='step'               , help='step'                                               , default=1,  type=int)
parser.add_option('-o', '--out'        ,    dest='output'             , help='output directory'                                   , default='/store/cmst3/group/hgcal/CMSSW/ntuples/')
parser.add_option('-c', '--cfg'        ,    dest='cfg'                , help='cfg file'                                           , default='test/runHGCSimHitsAnalyzer_cfg.py')
(opt, args) = parser.parse_args()


#prepare output
os.system('cmsMkdir %s'%opt.output)
jobsDir='%s/src/FARM%s'%(cmsswBase,time.time())
os.system('mkdir -p %s'%jobsDir)

from UserCode.HGCanalysis.storeTools_cff import fillFromStore
allFiles=fillFromStore('/store/cmst3/group/hgcal/CMSSW/'+opt.tag)
if opt.step<=0 : opt.step=len(allFiles)
outputTag=opt.tag.replace('/','_')
print outputTag
for ffile in xrange(0,len(allFiles),opt.step):
    #create a wrapper for standalone cmssw job
    scriptFile = open('%s/runJob_%s_%d.sh'%(jobsDir,outputTag,ffile), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('cd %s/src\n'%cmsswBase)
    scriptFile.write('eval `scram r -sh`\n')
    scriptFile.write('cd %s\n'%jobsDir)
    scriptFile.write('cmsRun %s/src/UserCode/HGCanalysis/%s %s %d %d\n'%(cmsswBase,opt.cfg,opt.tag,ffile,opt.step))
    scriptFile.write('cmsStage -f /tmp/psilva/%s_SimHits_%d.root %s\n'%(outputTag,ffile,opt.output))
    scriptFile.write('rm /tmp/psilva/%s_SimHits_%d.root\n'%(outputTag,ffile))
    scriptFile.write('echo "All done for this job\"\n')
    scriptFile.close()
    os.system('chmod u+rwx %s/runJob_%s_%d.sh'%(jobsDir,outputTag,ffile))

    #submit it to the batch or run it locally
    if opt.queue=='local':
        os.system('sh %s/runJob_%s_%d.sh'%(jobsDir,outputTag,ffile))
    else:
        os.system("bsub -q %s \'%s/runJob_%s_%d.sh\'"%(opt.queue,jobsDir,outputTag,ffile))
