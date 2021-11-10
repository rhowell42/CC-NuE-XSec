#col11 : rhcpt tune cv
#col11.1 : rhc weights alone.
#col11.2 : rhc post tune
#col11.3 : additional systemtatics for warping
#col11.4 : pre-tune
#col11.5 : post-tune 
#col11.6 : proper migration (hopefully)
stag="col11.6"
ntag="MAD"
#tarball="/pnfs/minerva/resilient/tarballs/hsu-gridselection-macro$(date +%H%M%S).tar.gz"
tarball="/minerva/data/users/hsu/tarballs/hsu-gridselection-macro$(date +%H%M%S).tar.gz"
if [[ -f $tarball ]]; then
    rm  $tarball
fi

common_flag="--selection_tag ${stag} --ntuple_tag ${ntag}  --use_sideband Excess_High_Inline Excess_Low_Inline Pi0 --tarball ${tarball} --extra_weighter RHCEelPt_tune"
common_flag+=" --exclude_universes SuSA_Valencia_Weight fsi_weight MK_model"
#common_flag+=" --exclude_universe all --scratch"
echo python gridSelection.py $common_flag
python gridSelection.py $common_flag --truth
#echo python gridSelection.py $common_flag --ext_2p2h --count 10
#python gridSelection.py $common_flag --ext_2p2h --count 10
echo python gridSelection.py $common_flag --NCDIF --count 2
python gridSelection.py $common_flag --NCDIF --count 1
