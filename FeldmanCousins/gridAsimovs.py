import datetime as dt
import os, time, sys, math
from config.GRIDConfig import gridargs,anaargs
from tools import Utilities

baseDir = os.path.dirname(os.path.abspath(__file__))+"/../"
MacroName = baseDir.split("/")[-4]
MAT = os.environ["PLOTUTILSROOT"].split("/")[-2]
CONFIG=os.environ["PYTHONPATH"].split(":")[0].split("/")[-1]
print(CONFIG)

def createTarball(outDir):
  print("I'm inside createTarball()")
  found = os.path.isfile(outDir)
  if(not found):
    cmd = "tar -czf %s -C %s %s"%(outDir, baseDir+"../", "{} {}".format(MacroName,MAT))# remove /Ana/CCNuE from baseDir bc want to tar the release.
    print(cmd)
    os.system(cmd)

  print("I'm done creating the tarballs")

def unpackTarball( mywrapper):
  # Add lines to wrapper that wil unpack tarball; add additional setup steps here if necessary  
  mywrapper.write("cd $CONDOR_DIR_INPUT\n")
  mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh\n")
  mywrapper.write("spack load root@6.28.12\n")
  mywrapper.write("spack load cmake\n")
  mywrapper.write("spack load gcc\n")
  mywrapper.write("spack load fife-utils@3.7.4\n")
  mywrapper.write("spack load py-numpy\n")
  mywrapper.write("tar -xvzf {}\n".format(outdir_tarball.split("/")[-1]))
  mywrapper.write("export MINERVA_PREFIX=`pwd`/{}\n".format(MAT))
  mywrapper.write("pushd {}/bin\n".format(MAT))
  mywrapper.write("source setup.sh\n")
  mywrapper.write("popd\n")
  mywrapper.write("pushd {}\n".format(MacroName))
  mywrapper.write("source setup_ccnue.sh {}\n".format(CONFIG))
  mywrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}\n")
  mywrapper.write("source FeldmanCousins/py3env/bin/activate\n")
  
  mywrapper.write("popd\n")

def addBashLine( wrapper , command ):

  wrapper.write("echo '---------------'\n")
  wrapper.write("echo '-------%s'\n" % command)
  wrapper.write("%s\n" % command)
  wrapper.write("echo '---------------'\n")


def submitJob(tupleName,tag):

  # Create wrapper
  wrapper_name = "grid_wrappers/%s/%s_wrapper.sh" % ( processingID , tupleName) 
  
  my_wrapper = open(wrapper_name,"w")
  my_wrapper.write("#!/bin/sh\n")
  unpackTarball(my_wrapper)

  # This is the bash line that will be executed on the grid
  my_wrapper.write( "cd $CCNUEROOT/FeldmanCousins\n")
  my_wrapper.write( "export USER=$(whoami)\n")
  my_wrapper.write( "source py3env/bin/activate\n")
  my_wrapper.write( 'py3env/bin/python3 fitAsimovs.py --grid --output $CONDOR_DIR_HISTS 2>> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.err 1>> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log\n %s' % (tupleName,tag,tupleName,tag,argstring))
  my_wrapper.write("exit $?\n")

  my_wrapper.close()
  
  os.system( "chmod 777 %s" % wrapper_name )
  
  cmd = "jobsub_submit --group=minerva -l '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest\\\"' --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --append_condor_requirements='CpuFamily != 6' --role=Analysis --memory %dMB -f %s -d HISTS %s -d LOGS %s -N %d --expected-lifetime=%dh  file://%s/%s" % ( memory , outdir_tarball , outdir_hists , outdir_logs , njobs, 36, os.environ["PWD"] , wrapper_name )
  os.system(cmd)

if __name__ == '__main__':
  PNFS_switch = gridargs.PNFS_switch
  # Automatically generate unique output directory
  processingID = '%s_%s-%s' % ("FHC_RHC", dt.date.today() , dt.datetime.today().strftime("%H%M%S") )
  os.system( "mkdir -p grid_wrappers/%s" % processingID )
  outdir_tarball=gridargs.tarball if gridargs.tarball else "/pnfs/minerva/resilient/tarballs/rhowell-%s.tar.gz" % (processingID)
  createTarball(outdir_tarball)

  if gridargs.memory is None:
    memory = 1000
  else:
    memory = gridargs.memory

  if gridargs.count is None:
    count = 1000
  else:
    count = gridargs.count

  njobs = 400
  if os.path.exists("grid_asimov.sh"):
    os.system( "rm grid_asimov.sh")

  for i in range(4,6):
    argstring=" ".join(anaargs)
    outdir_hists = "/pnfs/minerva/scratch/users/%s/%s_%s_Asimov_dchi2s_texts" % (os.environ["USER"],str(i),processingID)
    os.system( "mkdir -p %s" % outdir_hists )
    outdir_logs = "/pnfs/minerva/scratch/users/%s/%s_%s_Asimov_dchi2s_logs" % (os.environ["USER"],str(i),processingID)
    os.system( "mkdir -p %s" % outdir_logs )
    cmdString = "Asimovs_{}".format(i)
    tag = gridargs.ntuple_tag
    submitJob(cmdString,tag)
