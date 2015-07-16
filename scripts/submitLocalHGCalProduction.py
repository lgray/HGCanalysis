#!/usr/bin/env python

import os,sys
import optparse
import commands
import time


usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                                        , default='')
parser.add_option('-n', '--njobs'      ,    dest='njobs'              , help='number of jobs'                                     , default=1,  type=int)
parser.add_option('-o', '--options'    ,    dest='options'            , help='script options'                                     , default='')
parser.add_option('-s', '--script'     ,    dest='script'             , help='script to run (look under scripts/)'                , default='generateEventsFromCfi.sh')
(opt, args) = parser.parse_args()


#prepare working directory
cmsswBase=os.environ['CMSSW_BASE']
jobsDir=cmsswBase+'/src/FARM%s'%(time.time())
os.system('mkdir -p %s'%jobsDir)

#loop over the required number of jobs
for n in xrange(0,opt.njobs):

    #create a wrapper for standalone cmssw job
    scriptFile = open('%s/runJob_%d.sh'%(jobsDir,n+1), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('cd %s/src\n'%cmsswBase)
    scriptFile.write('eval `scram r -sh`\n')
    scriptFile.write('cd %s\n'%jobsDir)
    scriptFile.write('%s %s -j %d\n'%(opt.script,opt.options,n+1))
    scriptFile.close()
    os.system('chmod u+rwx %s/runJob_%d.sh'%(jobsDir,n+1))

    #submit it to the batch or run it locally
    if opt.queue=='':
        print 'Job #%d will run locally'%(n+1)
        os.system('%s/runJob_%d.sh'%(jobsDir,n+1))
    else:
        print 'Job #%d will run remotely'%(n+1)
        os.system("bsub -o /tmp/%s/std_out -q %s -R \"swp>2000 && pool>30000\" -J HGCSIM%d \'%s/runJob_%d.sh\'"%(os.environ['USER'],opt.queue,n+1,jobsDir,n+1))
    
