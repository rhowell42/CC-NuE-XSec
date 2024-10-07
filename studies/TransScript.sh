#!/bin/bash
#echo "You provided the arguments:" "$@"
#Scrip for run the TransWarpExtraction tool
#=========================================
RootFileDir="/exp/minerva/data/users/${USER}/nu_e"
#=========================================
#data is the distribution to be unfolded
#data_truth is the true unfolded distribution
data="Eavail_q3_bkgSubed"
data_truth="Eavail_q3_migration_truth"
#=========================================
#reco is the reco distribution that will be used to unfold the “data” dist.
#reco_truth is the truth distribution that will be used to unfold the “data” dist.
#migration is the migration matrix between the reco distand truth dist.
reco="Eavail_q3_migration_reco"
reco_truth="Eavail_q3_migration_truth"
migration="Eavail_q3_migration"
#=========================================
#DirData="/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcmastercombine_nx_new_fspline.root"
#DirReco="/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcmastercombine_nx_new_fspline.root"
#DirMigrationMatrix="/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcmastercombine_nx_new_fspline.root"

#DirData="/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcoptimization_nx_new_fspline.root"
#DirReco="/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcoptimization_nx_new_fspline.root"
#DirMigrationMatrix="/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcoptimization_nx_new_fspline.root"

#DirData="/exp/minerva/data/users/${USER}/nu_e/unfolded_full_col3_Sarah.root"
DirData="/exp/minerva/data/users/${USER}/nu_e/unfolded_fullfakedata_col3_Sarah.root"
DirDataTruth="/exp/minerva/data/users/${USER}/nu_e/kin_dist_datafullfakedata_col3_Sarah.root"
DirReco="/exp/minerva/data/users/${USER}/nu_e/kin_dist_mcfullfakedata_col3_Sarah.root"
DirMigrationMatrix="/exp/minerva/data/users/${USER}/nu_e/kin_dist_mcfull_col3_Sarah.root"

#DirData="/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcmastercombine_Noq3_nx_new_fspline.root"
#DirReco="/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcmastercombine_Noq3_nx_new_fspline.root"
#DirMigrationMatrix="/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcmastercombine_Noq3_nx_new_fspline.root"


iterations="1,2,3,4,5,6,7,8,10,15,20,25"
#iterations="1,2,3,5,8,9,10,20,50,100,120"
#iterations="1,2,3,4,5,6,7,8,9,10,11,12,13,14,20,26,32,38,44,50,90,130,144"
#iterations="1,2,3,4,5,6,7,8,9,10,11,12,13,14,20,26,32,38,44,50"
#iterations="1,2,3,4,5,6,7,8,9,10,20,40,80"

univers=100
#univers=200
#univers=200
#univers=25
#univers=25
#univers=100
#univers=25
dimension=2
#MaxChi2=100000
#MaxChi2=60000000
#stepChi2=500000

#MaxChi2=600
MaxChi2=1000
stepChi2=10
#stepChi2=500
remove=0
#remove="4,5,6,15,16"
#remove="20,24,6,10,13"
#pot_norm="0.38429"
#pot_norm="2.145365"
pot_norm="1"
#pot_norm="1"
#pot_norm="2.6"
pathExecutableOutput="q0_q3_2.0ExcludeBins2.root"

echo "Your argument fake_data is :" "$1"
echo "Your argument migration is :" "$2"

#./TransWarpExtraction -o $pathExecutableOutput -d $data -D $DirData -i $data_truth -I $DirData -m $migration -M $DirMigrationMatrix -r $reco -R $DirReco -t $reco_truth -T $DirReco -n $iterations -z $dimension -u $univers -c $MaxChi2 -C $stepChi2 -x $remove
#./TransWarpExtraction -o $pathExecutableOutput -d $data -D $DirData -i $data_truth -I $DirData -m $migration -M $DirMigrationMatrix -r $reco -R $DirReco -t $reco_truth -T $DirReco -n $iterations -z $dimension -u $univers -c $MaxChi2 -C $stepChi2 --data_pot_norm $pot_norm -x $remove
#./TransWarpExtraction -o $pathExecutableOutput -d $data -D $DirData -i $data_truth -I $DirData -m $migration -M $DirMigrationMatrix -r $reco -R $DirReco -t $reco_truth -T $DirReco -n $iterations -z $dimension -u $univers -c $MaxChi2 -C $stepChi2 --data_pot_norm $pot_norm --exclude_bins 17,25,34,35,44,53

./TransWarpExtraction -o $pathExecutableOutput -d $data -D $DirData -i $data_truth -I $DirDataTruth -m $migration -M $DirMigrationMatrix -r $reco -R $DirReco -t $reco_truth -T $DirReco -n $iterations -z $dimension -u $univers -c $MaxChi2 -C $stepChi2 -P $pot_norm --exclude_bins $remove

#echo -o $pathExecutableOutput -d $data -D $DirData -i $data_truth -I $DirData -m $migration -M $DirMigrationMatrix -r $reco -R $DirReco -t $reco_truth -T $DirReco -n $iterations -z $dimension -u $univers -c $MaxChi2 -C $stepChi2 --data_pot_norm $pot_norm

#./TransWarpExtraction -o $pathExecutableOutput -d $data -D $DirData -i $data_truth -I $DirData -m $migration -M $DirMigrationMatrix -r $reco -R $DirReco -t $reco_truth -T $DirReco -n $iterations -z $dimension -u $univers -c $MaxChi2 -C $stepChi2 --data_pot_norm $pot_norm --exclude_bins  
#./TransWarpExtraction -o $pathExecutableOutput -d $data -D $DirData -i $data_truth -I $DirData -m $migration -M $DirMigrationMatrix -r $reco -R $DirReco -t $reco_truth -T $DirReco -n $iterations -z $dimension -u $univers -c $MaxChi2 -C $stepChi2 --data_pot_norm $pot_norm 
#./TransWarpExtraction -o $pathExecutableOutput -d $data -D $DirData -i $data_truth -I $DirData -m $migration -M $DirMigrationMatrix -r $reco -R $DirReco -t $reco_truth -T $DirReco -n $iterations -z $dimension -u $univers -c $MaxChi2 -C $stepChi2 --data_pot_norm $pot_norm
#./TransWarpExtraction -o $pathExecutableOutput -d $data -D $DirData -i $data_truth -I $DirData -m $migration -M $DirMigrationMatrix -r $reco -R $DirMigrationMatrixTest -t $reco_truth -T $DirMigrationMatrixTest -n $iterations -z $dimension -u $univers -c $MaxChi2 -C $stepChi2 --data_pot_norm $pot_norm
