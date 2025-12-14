# Macro Files Structure:
Here is a list of files currently used in the macro:
- file_option/: store playlist files
    - playlist.txt: txt files that consist of paths of ntuple files of given playlist
    
- tools/: store codes that do specific tasks, intended to be share between different nu_e analyses.
    - SystematicsUniverse.py : define a variety of systematics universes using centralized C++ CVUniverses, supersede CVUniverse.py.
    - CVUniverse.py : obsolete CVUniverse.
    - mcweight.py : obsolete mc weights.
    - CutLibary.py : define a variety of selection cuts.
    - EventClassification.py : classify reco events according to selection cuts and classify truth events according to signal defination.
    - MyHistograms.py : define PlotProcessor to handle filling different histograms.
    - PlotLibary.py : define a variety of plots, containing settings to set up a PlotProcessor to fill requested plot.
    - PlotTools.py: functions that used to make plots.
    - Utilities.py: utility functions that might be used in more than one place.
    - KinematicsCalculator.py: Calculate reco/truth event kinematics variables.
    
- config/: configure the macro to fit different analysis, intended to be different for each analysis
    - AnalysisConfig.py: define default and command line options for event selection and probably background fits, unfolding, etc
    - CutConfig.py: specify cuts to be used for signal region and each side band. cuts should be defined in tools/CutLibrary
    - GRIDConfig.py : define default and command line options for running event seletions in GRID.
    - POT.json : specify location of playlist file for different playlist, version of ntuple, data or mc. store POT infomation of given playlist file for future usage.
    - PlotConfig.py : configure binnings of plots in tools/PlotLibary.py, histograms to be make in eventSelection, and color/name of different truth categories. 
    - SignalDef.py: define the categories to be classify in tools/EventClassification.py. 
    - SystematicsConfig.py : control variable shifts in tools/SystematicsUniverse, and grouping of different Errorbands when plotting with systematic errors
    
- selection/: entry points for event selection, hopefully to be the same for all nu_e analyses.
    - eventSelection.py: entry point to perform event selection (i.e. `python eventSelection.py --options`)
    - plotSelection.py: entry points to make plots for selected events (i.e. `python plotSelection.py --options`)
    - gridSelection.py: enntry points to perform event seletion on grid (i.e. `python gridSelection.py --options`)
    
-background_fit/: entry points for background fittings, hopefully to be the same for all nu_e analysis.
    - backgroundFit.py: entry point to perform backgroundFit (i.e. python backgroundFit.py --options)
    

# Changes between v1.2 to v1.1
## Changes in Ana stage
- Added variables for birk's constant shifted universes, birk's constant shifted dEdX and  electron score.
- Relaxed electron score cut to 0.6 because birk's constant shifted electron score could be 0.7, don't want to cut these events.
- Relaxed primayprong cut == 1 cut to accommodate events that have more than one candidate due to relaxed cut (e.g. two prongs with score  0.7 and 0.6)
- Sorted the electron candidate prongs in desending order of electron score for selecting best candidate for macro stage.
- Put back number of startpoint >=1 cut, otherwise startpoint match cut will crash

## Changes in macro stage
### Underlying mechinery
- Changed systematics universe implementation from my customized version to centralized c++ version. Changed interface in eventSelection as well
- Changed delivery of electron momentum in cv universes from TVector4D to Energy and Theta, eliminated unused kinematics calculation and changed TVector4D to modern GenVector implementation, for improving efficiency of event selection.
- Merged PlotProcessor1D and PlotProcessor2D. Added cut options for PlotProcessors. Changed PLOT_SETTING formats accordingly. Name of plots are modified to avoid empty space.
- Added tags to PLOT_SETTING, for creating a bunch of similar plots.
- Splitted PlotConfig to PlotConfig and PlotLibary, and splitted CutCong_inc to CutConfig and CutLibrary, for sharing Libaries among different analyses.
- Introduced config/SignalDef.py for easily adding/removing categories, modified EventClassification accodingly.


### User interface
- Control knobs in eventSelection were moved to AnalysisConfig(I/O paths, playlist nicknames, etc), PlotConfig(HISTS_TO_MAKE), or other places in /config. Added commandline options for controling these knobs. use `--help` to check available options.
- Removed the ability of running over multiple playlists of eventSelection. Running over multiple playlists should be done by gridSelection.
- CCNuE/macros/Low_Recoil/ would be added to python path by setting up CCNuE package, and was hard-coded in gridSelection. As a result, users are recommanded to work on Low_Recoil directory rather than create thier own directory. Add back `sys.path.insert(0,os.environ["CCNUEROOT"]+"path/to/your/low_recoil")` to eventSelection and change path in Line 74 of gridSelection.py if user want to use a different directory
- Modified the playlist files to xrootd urls for GRID access.


# Notes on using gridSelection
- By default, gridSelection assume the mc files are merged by official merge tool to one file per run, and not merging data files. As a result, it submit a job per 1000 data files and a job per 1 mc file (500 mc files if unmerged). Users can change number of files per job by `--count` option. It is recommanded to merge mc files to imporve efficiency.
- The output file/log live in users' dCache persisent area, named by CCNuE_selection_YYYY_MM_DD_HH00_hists and CCNuE_selection_20XX_MM_DD_HH00_logs. Users can merge output root files by `madd.exe`. It is shipped with PlotUtils and the usage is the same as `hadd`
- Be sure to remove previous output directories if you changed the code and resubmit in an hour, the dCache files can't be overwritten.
- most options for eventSelection will work with gridSelection, even though it won't show up in `--help`. The help message shows options that overwritten in gridSelection, unrecognized option will be forwarded to eventSeletion
- The grid processing runs a playlist in separated jobs, hence counting POT in each job doesn't make sense. User should used `--cal_POT` option of gridSelection to count POT.

# Changes from AL9 Upgrade
- This framework now supports AL9! Specifically, it supports more modern versions of ROOT that are built with Python3 and change how multi-cross inheritance works.

# How to Install in Alma9 OS
## Install MAT-MINERvA
```
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load root@6.28.12
spack load cmake
spack load gcc
spack load fife-utils
export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}

mkdir /exp/minerva/app/users/$USER/MAT_AL9 && cd /exp/minerva/app/users/$USER/MAT_AL9
git clone https://github.com/MinervaExpt/MAT.git
git clone https://github.com/MinervaExpt/MAT-MINERvA.git
git clone https://github.com/MinervaExpt/UnfoldUtils.git

mkdir opt && cd opt #install prefixed for OPTimized version of MAT/tutorial libraries and executables

#Build MAT-MINERvA which builds MAT & UnfoldUtils, and also finds the centralized flux files
mkdir buildMAT-MINERvA && cd buildMAT-MINERvA
cmake ../../MAT-MINERvA/bootstrap -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Release #installs libraries in ".."
make install
source /exp/minerva/app/users/$USER/MAT_AL9/opt/bin/setup.sh

cd /exp/minerva/app/users/$USER/MAT_AL9/
git clone https://github.com/rhowell42/CC-NuE-XSec.git
cd CC-NuE-XSec/CVUniversePythonBinding/
rm -r build && mkdir build/ && cd build/
cmake ../ -DCMAKE_INSTALL_PREFIX=/exp/minerva/app/users/$USER/MAT_AL9/opt -DCMAKE_BUILD_TYPE=Release
make install
```
The following shell script can be run on login to the gpvms to set everything up:
```
source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh

spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3
spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3
spack load gcc
spack load python@3.9.15
spack load ifdhc-config@2.6.20%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3
spack load py-numpy@1.24.3%gcc@12.2.0
export JOBSUB_GROUP=minerva

htgettoken -a htvaultprod.fnal.gov -i minerva
export BEARER_TOKEN_FILE=/run/user/`id -u`/bt_u`id -u`

source /exp/minerva/app/users/$USER/MAT_AL9/opt/bin/setup.sh

export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
cd /exp/minerva/app/users/$USER/MAT_AL9/CC-NuE-XSec/
```

# Making an Event Selection
## Running Interactively
- `python eventSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband dEdX --truth --cal_POT --selection_tag thesis --mc_only`
    - `--playlist me5A_p4`  tells you to run over the anatuples in the me5A playlist, which is a MINERvA RHC run. This specifically tells you to use the tuples whose full paths are stored in a dictionary in configs/POT.json which points to the files in file_option/
    - change POT.json to point to text files containing paths to any additional anatuples you want to run over, and add the text files in file_option/ to actually point to those LE files. See files in those areas for examples.
    - `--ntuple_tag MAD`  further identifies which files you want to use in configs/POT.json, it specifies the tool used to create the anatuples, in this case (and for yours) it stands for MasterAnaDev
    - `--use-sideband dEdX`  uses a sideband named dEdX which is defined with a series of kinematic cuts in $CONFIGPATH/config/CutConfig.py.
    - `--selection_tag thesis` puts a string at the end the output root files to further help identify them, this case "thesis"
    - `--mc_onl`y tells you to only run over the monte carlo production. Just remove this to run with both data and MC, or put `--data_only` if you just want to run over data.
    - `--test` put this in your command to just run over 1000 events to test that everything is working before you devote a lot of time to running over the full thing
## Running on the Grid
- Enter an SL7 container and run through your setup same as before
- `python gridSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband dEdX --truth --cal_POT --selection_tag thesis --mc_only` , where you can just use the same command line options as you used for the interactive run (without the `--test` option)
- This will, among other things, create a tarball of your environment that will be used to run the event selection on a grid node, and create the jobsub_submit commands to actually submit those jobs
- This will spit out a bash script in jobsub_commands/ called jobsub_wrapper.sh. Because of recent OS changes from SL7 to Alma9, you can't actually submit jobs inside an SL7 container. What you can do is run jobsub_submit commands in a native Alma9 environment (see next step)
- Exit your SL7 container and submit your jobs by running source jobsub_wrapper.sh. If everything has worked this will submit your jobs and you can check they're running with jobsub_q -G minerva --user $USER some time later. These take up to an hour to finish running all the way, usually.
## Combine Files
- `python combine_file.py --playlist NewPlaylistName --i  /pnfs/minerva/persistent/users/$USER/{selection_directory}_hists/ --cal_POT --ntuple_tag MAD --selection_tag thesis`
- You pass it a new playlist name that you create yourself, to represent this fully combined event selection. For me I named them to denote them as either RHC or FHC samples. Make sure you pass the same playlist name for MC and data selections of the same sample.
- This will prompt you for additional paths if you want to add more samples. If you're running this over your base MC selection sample, this is where you would add the path to your special samples selection.
- If you're running over data you don't need to do anything here, just press enter again.
- This merges all of your samples together with madd, which is the MINERvA wrapper for hadd that takes care of all the universe sub-histograms in your event selection. It spits out one root file for each run and saves it to output_dir, which is set in your analysis configuration. I think you can pass this on the command line as well if you want to change this in situ.
## Background Scaling
- Run background tuning and scaling: this tunes your signal region and sideband region to data in a histogram type provided by you in . This takes the output root file from the event selection as the input, and uses the histograms in the list HISTOGRAMS_TO_UNFOLD in $CONFIGPATH/config/UnfoldingConfig.py.
- `python background_fit/backgroundFit.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband dEdX --truth --selection_tag thesis --bkgTune_tag N4_tune`
- `--bkgTune_tag N4_tune` refers to which a dictionary object in $CONFIGPATH/config/BackgroundConfig.py that describes the categories with which to tune to data, as well as a few other things
- This outputs root files containing your background subtracted data, the scaled-to-data histograms, and the base MC prediction with which to compare to the background subtracted data.
- This will also make plots as given in $CONFIGPATH/config/DrawingConfig.py on the scaled and background subtracted histograms

# Sterile Neutrino Oscillations
The framework to perform sterile neutrion oscillation analyses is kept in oscillations/. This was built with Medium Energy in mind, but it could be modifed to work on Low Energy samples as well. This is heavily reliant on python libraries that are not available on the gpvms, so a virtual environment is used.
```
cd oscillations/
python3 -m venv py3env
source py3env/bin/activate
pip install -r requirements.txt
deactivate
```
## Special Samples
- If accounting for numu -> nue oscillations, one needs to have a flavor swapped sample available of nue events with the parent decay positions and energies of numu events (and polarizations and other applied conservations)
- There exists a sample for the Medium Energy sample, of flavor swapped ME1A, ME5A, and ME6A, at roughly 1/20 total POT. One may be able to reweight this to form a modified flavor swapped sample to use for Low Energy studies
- The same selection should be applied to the flavor swapped sample and care should be made to properly POT scale it when added back in for an oscillation hypothesis
## Configuration
- You need two histograms to perform an oscillation analysis: "Biased Neutrino Energy", and "Reco Energy vs L/E"
    - "Biased Neutrino Energy" is the 1D reconstructed energy estimator, binned in GeV in unequal bin widths from 0 to 20 GeV
    - "Reco Energy vs L/E" is the 2D template used to assigned oscillation probabilities, it gives the energy estimator and the corresponding true L/E
- If running an electron neutrino selection, you also need to run with the dE/dX sideband for the background scaling and subtraction procedure. This is not necessary for the flavor swapped sample, since we take the base MC prediction of both samples
- The cuts to select the low-nu-like selections are defined in each of the corresponding config/CutConfig.py files. This is true for the FHC, RHC, nue, and numu selections
- Numu selections do not need to be background scaled/subtracted because there's neglibible background events for CC numu in the NuMI beam
## Samples
- Many samples used for oscillation measurements in MINERvA have strong correlations between each other, which must be taken into account when computing the chi2. A StitchedHistogram object is used to combine every sample together into one histogram, with correctly handled systematics, to get an accurate full covariance matrix.
- `py3env/bin/python3 stitchHistograms.py`
    - This will combine all the samples you want for an oscillation analysis. Some options can be given in the command line for different samples; `--ratio` to use flavor ratios instead of base selections, `--exclude "{Your Samples}"` to exclude any samples from an analysis, `--fit_muons` to keep in muon neutrino selections when using flavor ratios, etc...
- `py3env/bin/python3 fitData.py` fit an oscillation analysis to data, returns best fit parameters and test statistic
- `py3env/bin/python3 gridSurface.py` send jobs to the grid to compute sensitivity and exclusion regions in oscillation space for your sample
    - `merge_files.py` to merge these files in the proper order to make plotting easier (see https://github.com/rhowell42/FeldmanCousinsMacros)  
- `py3env/bin/python3 gridAsimovs.py` send jobs to the grid to compute 100,000 pseudo-experiment delta chi2s to determine Feldman Cousins contours for null hypothesis
