#!/bin/sh
jobsub_submit --group=minerva -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --memory 5000MB -f /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_1000Universes.tar.gz -d HISTS /pnfs/minerva/persistent/users/rhowell/CCNUE_selection_2024-07-10-135129_hists -d LOGS /pnfs/minerva/persistent/users/rhowell/CCNUE_selection_2024-07-10-135129_logs -N 20 --expected-lifetime=36h  file:///exp/minerva/app/users/rhowell/cmtuser/CCNue/selection/grid_wrappers/CCNUE_selection_2024-07-10-135129/CCNuE-me5A_swap-mc_wrapper.sh
