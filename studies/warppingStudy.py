import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process

from tools.PlotLibrary import PLOT_SETTINGS,HistHolder
from config.AnalysisConfig import AnalysisConfig
from config.BackgroundFitConfig import CATEGORY_FACTORS
from config.SignalDef import SIGNAL_DEFINATION
from config import PlotConfig
from tools import Utilities,PlotTools
from tools.unfoldingwrapper import DeOverflowWrapper,WithBackgroundWrapper,IdentityWrapper

ROOT.TH1.AddDirectory(False)
baseDir = os.path.dirname(os.path.abspath(__file__))+"/../"
MacroName = baseDir.split("/")[-4]
MAT = os.environ["PLOTUTILSROOT"].split("/")[-2]
CONFIG=os.environ["PYTHONPATH"].split(":")[0].split("/")[-1]

def createTarball(outDir):
  print("I'm inside createTarball()")
  found = os.path.isfile(outDir)
  if(not found):
    cmd = "tar -czf %s -C %s %s"%(outDir, baseDir+"../", "{} {}".format(MacroName,MAT))# remove /Ana/CCNuE from baseDir bc want to tar the release.
    print(cmd)
    os.system(cmd)

  print("I'm done creating the tarballs")

Iterations = "1,2,3,4,5,6,7,8,9,10,15,20,30,40"
#Iterations= "1"
NUniverses = 100
MaxChi2 =600
stepChi2 = MaxChi2//200
#stepChi2 = 10
remove=""
NUniverses_per_job = NUniverses
K=0
threshold = 5
Sys_on_Data=False
#Additional_uncertainties = [0,1,2,3,4,5,6,7,8,9,10]
unfoldwithFakes = False
LOCAL=True

def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None

def unpackTarball( mywrapper):
  # Add lines to wrapper that wil unpack tarball; add additional setup steps here if necessary  
  mywrapper.write("cd $CONDOR_DIR_INPUT\n")
  mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/products/setup\n")
  mywrapper.write("setup root v6_22_06a -q e19:p383b:prof\n")
  mywrapper.write("tar -xvzf {}\n".format(outdir_tarball.split("/")[-1]))
  mywrapper.write("export MINERVA_PREFIX=`pwd`/{}\n".format(MAT))
  # Set up test release
  #mywrapper.write("pushd Tools/ProductionScriptsLite/cmt/\n")#assuming only one test release
  #mywrapper.write("cmt config\n")
  #mywrapper.write(". setup.sh\n")
  #mywrapper.write("popd\n")
  mywrapper.write("pushd {}/bin\n".format(MAT))
  mywrapper.write("source setup.sh\n")
  #mywrapper.write("cmt make\n")
  #mywrapper.write(". setup.sh\n")
  mywrapper.write("popd\n")
  #mywrapper.write("export LD_LIBRARY_PATH=/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/GSL/1.10/x86_64-slc7-gcc49-opt/lib:$LD_LIBRARY_PATH\n")
  mywrapper.write("pushd {}\n".format(MacroName))
  mywrapper.write("source setup_ccnue.sh {}\n".format(CONFIG))
  #mywrapper.write("make\n")
  
  mywrapper.write("popd\n")

def submitJob( tupleName,tag,njobs):

  # Create wrapper
  wrapper_name = "grid_wrappers/%s/%s_wrapper.sh" % ( processingID , tupleName )

  my_wrapper = open(wrapper_name,"w")
  my_wrapper.write("#!/bin/sh\n")
  unpackTarball(my_wrapper)
    # Write dummy output files to get around dumb copy-back bug
  #makeDummyFile(my_wrapper,'$CONDOR_DIR_LOGS')
  #makeDummyFile(my_wrapper,'$CONDOR_DIR_HISTS')


  # This is the bash line that will be executed on the grid
  my_wrapper.write( "export USER=$(whoami)\n")
  #my_wrapper.write( "export XRD_LOGLEVEL=\"Debug\"\n")

  #my_wrapper.write("readlink -f `which TransWarpExtraction` &> $CONDOR_DIR_OUTPUT/log_{tag}_$PROCESS\n".format(tag=tag)) 
  pot_scale=1
  my_wrapper.write("TransWarpExtraction -o $CONDOR_DIR_OUTPUT/transWrap_{tag}_{PROCESS}.root -d data -D $CONDOR_DIR_INPUT/{iput} -i data_truth -I $CONDOR_DIR_INPUT/{iput} -m migration -M $CONDOR_DIR_INPUT/{iput} -r reco -R $CONDOR_DIR_INPUT/{iput} -t reco_truth -T $CONDOR_DIR_INPUT/{iput} -n {Iter} -z {dimension} -u {NUniverse} -c {max_chi2} -C {step_chi2} -b -f {PROCESS} -L {remove} >> $CONDOR_DIR_OUTPUT/log_{tag}_{PROCESS}\n".format(tag=tag,iput=ifile,Iter=Iterations,dimension=dimension,NUniverse=NUniverses_per_job,max_chi2=MaxChi2,step_chi2=stepChi2,remove="-x {}".format(remove) if remove else "", PROCESS="$PROCESS" if njobs>1 else "-1"))

  my_wrapper.write("exit $?\n")
  #my_wrapper.write( "python eventSelection.py -p %s --grid --%s-only --ntuple_tag %s --count %d %d  --output $CONDOR_DIR_HISTS %s \n" % (playlist, dataSwitch, gridargs.ntuple_tag, start, count, argstring) )

  #addBashLine(my_wrapper,'date')
  my_wrapper.close()

  os.system( "chmod 777 %s" % wrapper_name )
  cmd = "jobsub_submit --group=minerva -l '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\\\"' --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --use-pnfs-dropbox --memory %dMB --tar_file_name dropbox://%s -d OUTPUT %s %s --expected-lifetime=%dh -f dropbox://%s/%s file://%s/%s" % ( memory , outdir_tarball , outdir_hists , "-N {}".format(njobs) if njobs>1 else "", 4,  os.environ["PWD"],ifile , os.environ["PWD"] , wrapper_name )
  print(cmd)
  os.system(cmd)

def findexcludebins(h,threshold):
  result = []
  for i in range(h.GetSize()):
    if 0< h.GetBinContent(i)<threshold:
      result.append(str(i))
  print(result)
  return ",".join(result)

def ConsolidateInputs(Temp_path,mig,reco,reco_truth,data,data_truth):
  def scale(hi): 
    h = hi.Clone("{}_clone".format(name))
    if bignue_file == data_file:
      h.Scale(17379.0/h.Integral(0,-1,0,-1))
    else:
      h.Scale(1e21/mc_pot)
    for i in range(h.GetSize()):
      h.SetBinError(i,math.sqrt(h.GetBinContent(i)))
    return h
  scale_dict = {
    "data","data_truth",
    #"reco","reco_truth","migration" 
  }
  f = ROOT.TFile(Temp_path,"RECREATE")
  for name,h in [("data",data),("data_truth",data_truth),("migration",mig),("reco_truth",reco_truth),("reco",reco)]: 
    if name in scale_dict:
      tmp = scale(h)
      tmp.Write(name)
      if name =="data_truth":
          rm = findexcludebins(tmp,threshold) 
    else:
      h.Write(name)
  f.Close()
  return rm

def RunTransWrapper(ifile,ofile,njobs):
  cmd = "TransWarpExtraction -o {output} -d data -D {iput} -i data_truth -I {iput} -m migration -M {iput} -r reco -R {iput} -t reco_truth -T {iput} -n {Iter} -z {dimension} -u {NUniverse} -c {max_chi2} -C {step_chi2} {remove} -L -b".format(output=ofile,iput=ifile,Iter=Iterations,dimension=dimension,NUniverse=NUniverses,max_chi2=MaxChi2,step_chi2=stepChi2,remove="-x {}".format(remove) if remove else "")
  print(cmd)
  os.system(cmd)


AlternativeModels = {
   "CV":None,
   "LowQ2Pi0": ("LowQ2Pi",0),
   "LowQ2Pi1": ("LowQ2Pi",1),
   "LowQ2Pi2": ("LowQ2Pi",2),
   "LowQ2Pi3": ("LowQ2Pi",3),
   "2p2h0": ("Low_Recoil_2p2h_Tune",0),
   "2p2h1": ("Low_Recoil_2p2h_Tune",1),
   "2p2h2": ("Low_Recoil_2p2h_Tune",2),
   "RPA_highq20" :("RPA_HighQ2",0),
   "RPA_highq21" :("RPA_HighQ2",1),
   "RPA_lowq20" :("RPA_LowQ2",0),
   "RPA_lowq21" :("RPA_LowQ2",1),
   "MK_Model":("MK_model",0),
   "FSI_Weight0":("fsi_weight",0),
   "FSI_Weight1":("fsi_weight",1),
   "FSI_Weight2":("fsi_weight",2),
   "SuSA2p2h":("SuSA_Valencia_Weight",0),
   "GenieMaCCQE_UP":("GENIE_MaCCQE",1),
   "GenieMaCCQE_DOWN":("GENIE_MaCCQE",0),
   "GenieMaRES":("GENIE_MaRES",0),
   "GenieRvn2pi":("GENIE_Rvn2pi",0), 
}

Measurables = [
  "Eavail_q3", 
  "Eavail_Lepton_Pt"
]

if __name__ == '__main__':
  playlist=AnalysisConfig.playlist
  type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
  data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True) 
  bignue_file = ROOT.TFile.Open('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcNueOnly_nx_collab1_fspline.root')
 
  if bignue_file:
    data_file = bignue_file
    mc_file = bignue_file 
  #  #data_file = bignue_file
    print("Using BigNuE file for MC")

  #unfolded_file = ROOT.TFile.Open(AnalysisConfig.UnfoldedHistoPath(playlist,False))
  ifile="fin.root" 
  PNFS_switch = "scratch"
  # Automatically generate unique output directory
  processingID = '%s_%s-%s' % ("CCNUE_Warping", dt.date.today() , dt.datetime.today().strftime("%H%M%S") )
  
  outdir_hists = "/pnfs/minerva/%s/users/%s/%s_hists" % (PNFS_switch,os.environ["USER"],processingID)
  os.system( "mkdir %s" % outdir_hists )
  #outdir_logs = "/pnfs/minerva/%s/users/%s/%s_logs" % (PNFS_switch,os.environ["USER"],processingID)
  #os.system( "mkdir %s" % outdir_logs )
  os.system( "mkdir -p grid_wrappers/%s" % processingID )
  outdir_tarball="/pnfs/minerva/resilient/tarballs/shenry-%s.tar.gz" % (processingID)
  #createTarball(outdir_tarball)
  memory =1500
  for prefix in Measurables:
    MnvMigration = Utilities.GetHistogram(mc_file, prefix+ "_migration")
    MnvMCreco = Utilities.GetHistogram(mc_file, prefix)
    MnvMCtruth = Utilities.GetHistogram(mc_file, prefix+ "_migration_truth")
    if not MnvMCtruth: 
      MnvMCtruth = MnvMigration.ProjectionY("_truth",0,-1) 
    MnvMCSignal = Utilities.GetHistogram(mc_file, prefix+ "_migration_reco") 

    if not MnvMCSignal: 
      MnvMCSignal = MnvMigration.ProjectionX("_reco",0,-1)
    dimension = MnvMCreco.GetDimension()
    MnvDataMigration = Utilities.GetHistogram(data_file, prefix+ "_migration")
    if not MnvDataMigration:  
      MnvDataMigration=MnvMigration
      MnvDatareco = MnvMCreco
      MnvDatatruth = MnvMCtruth
      MnvDataSignal = MnvMCSignal
    else:
      MnvDatareco = Utilities.GetHistogram(data_file,prefix)
      MnvDatatruth = Utilities.GetHistogram(data_file, prefix+ "_migration_truth") 
      if not MnvDatatruth: 
        MnvDatatruth = MnvDataMigration.ProjectionY(prefix + "_truth",0,-1)
      MnvDataSignal = Utilities.GetHistogram(data_file, prefix+ "_migration_reco")
      if not MnvDataSignal: 
        MnvDataSignal = MnvDataMigration.ProjectionX("_reco",0,-1)
   
    if dimension ==1:
      H2D = PlotUtils.MnvH1D
    else:
      H2D = PlotUtils.MnvH2D
      sumh = lambda h:h.Integral(0,-1,0,-1)
      
   
    for model,value in AlternativeModels.items():
      data_truth = H2D(MnvDatatruth.GetCVHistoWithStatError())
      data_reco = H2D(MnvDatareco.GetCVHistoWithStatError())
      data_signal = H2D(MnvDataSignal.GetCVHistoWithStatError())
      model_migration = PlotUtils.MnvH2D(MnvMigration.GetCVHistoWithStatError())
      model_truth=H2D(MnvMCtruth.GetCVHistoWithStatError())
      model_signal=H2D(MnvMCSignal.GetCVHistoWithStatError())
      model_reco=H2D(MnvMCreco.GetCVHistoWithStatError())
      if model!="CV":
        try:
          if Sys_on_Data: 
            data_truth = H2D(MnvDatatruth.GetVertErrorBand(value[0]).GetHist(value[1]))
            data_signal= H2D(MnvDataSignal.GetVertErrorBand(value[0]).GetHist(value[1]))
            data_reco= H2D(MnvDatareco.GetVertErrorBand(value[0]).GetHist(value[1]))
          else:
            model_migration = PlotUtils.MnvH2D(MnvMigration.GetVertErrorBand(value[0]).GetHist(value[1]))
            model_truth = H2D(MnvMCtruth.GetVertErrorBand(value[0]).GetHist(value[1]))
            model_signal= H2D(MnvMCSignal.GetVertErrorBand(value[0]).GetHist(value[1]))
            model_reco= H2D(MnvMCreco.GetVertErrorBand(value[0]).GetHist(value[1]))
        except TypeError as e:
          print(e)
          continue
  
      if unfoldwithFakes:
        remove=ConsolidateInputs(ifile,model_migration,model_reco,model_truth,data_reco,data_truth)
      else:
        remove=ConsolidateInputs(ifile,model_migration,model_signal,model_truth,data_signal,data_truth)

      #ConsolidateInputs(ifile,model_migration,model_signal,model_truth,data_signal,data_truth)
      njobs = 11#int(math.ceil(1.0*NUniverses/NUniverses_per_job))
      cmdString = "CCNuE-%s-%s" % (playlist,"mc")
      #wrapper.setup(migration,cv_reco,cv_bkg)
      #MyUnfoldWrapper(ifile)
      #MnvUnfoldWrapper(ifile)
      if not LOCAL:
        submitJob(cmdString,model+prefix,1)
        submitJob(cmdString,model+prefix,3)
      else:
        RunTransWrapper(ifile,"transWarp_{}.root".format(model+prefix),-1)

      #submitJob(cmdString,model+prefix)
      #RunTransWrapper(ifile,"transWarp_{}.root".format(model+prefix))
      #break
      # mnvunfold=UnfoldUtils.MnvUnfold()
      # h_test = MnvMCtruth.Clone()
      # cov=ROOT.TMatrixD(1,1)
      # mnvunfold.UnfoldHistoWithFakes(h_test,cov,MnvMigration,data,cv_signal,cv_truth,cv_bkg,5,False,False)

 
