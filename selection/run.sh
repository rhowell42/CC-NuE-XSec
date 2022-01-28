stag="col12.4"
ntag="MAD"
#tarball="/pnfs/minerva/resilient/tarballs/hsu-gridselection-macro$(date +%H%M%S).tar.gz"
tarball="/minerva/data/users/hsu/tarballs/hsu-gridselection-macro$(date +%H%M%S).tar.gz"
if [[ -f $tarball ]]; then
    rm  $tarball
fi

common_flag="--selection_tag ${stag} --ntuple_tag ${ntag}  --use_sideband Excess_High_Inline Excess_Low_Inline Pi0 --extra_weighter FHCPt_tune2 --tarball ${tarball}"
#common_flag+=" --exclude_universes SuSA_Valencia_Weight fsi_weight MK_model"
#common_flag+=" --exclude_universe all --scratch"
echo python gridSelection.py $common_flag
python gridSelection.py $common_flag
#diffractive pi0 sample
echo python gridSelection.py $common_flag --NCDIF --count 2
python gridSelection.py $common_flag --NCDIF --count 1
#signal rich sample
#python gridSelection.py $common_flag --bignue --truth --count 1
