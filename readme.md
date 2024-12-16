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

# How to make a selection
- `python eventSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband dEdX --truth --cal_POT --selection_tag thesis --mc_only`
  -- `--playlist me5A_p4`  tells you to run over the anatuples in the me5A playlist, which is a MINERvA RHC run. This specifically tells you to use the tuples whose full paths are stored in a dictionary in configs/POT.json which points to the files in file_option/
    --- change POT.json to point to text files containing paths to any additional anatuples you want to run over, and add the text files in file_option/ to actually point to those LE files. See files in those areas for examples.
  -- `--ntuple_tag MAD`  further identifies which files you want to use in configs/POT.json, it specifies the tool used to create the anatuples, in this case (and for yours) it stands for MasterAnaDev
 -- `--use-sideband dEdX`  uses a sideband named dEdX which is defined with a series of kinematic cuts in $CONFIGPATH/config/CutConfig.py.
-- `--selection_tag thesis` puts a string at the end the output root files to further help identify them, this case "thesis"
-- `--mc_onl`y tells you to only run over the monte carlo production. Just remove this to run with both data and MC, or put `--data_only` if you just want to run over data.
-- `--test` put this in your command to just run over 1000 events to test that everything is working before you devote a lot of time to running over the full thing
