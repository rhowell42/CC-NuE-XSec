#12.5 include 50 Mev bin , 12.6 without
stag="col13.0"
ntag="MAD"
#tarball="/pnfs/minerva/resilient/tarballs/hsu-gridselection-macro$(date +%H%M%S).tar.gz"
tarball="/minerva/data/users/hsu/tarballs/hsu-gridselection-macro$(date +%H%M%S).tar.gz"
if [[ -f $tarball ]]; then
    rm  $tarball
fi

common_flag="--selection_tag ${stag} --ntuple_tag ${ntag} --tarball ${tarball}"
common_flag+=" --extra_weighter rhc_weight --exclude_universes bkg_tune"
#common_flag+=" --extra_weighter CV_tune"
#common_flag+=" --exclude_universes SuSA_Valencia_Weight fsi_weight MK_model"
#common_flag+=" --exclude_universe all --scratch"
echo python gridSelection.py $common_flag
python gridSelection.py $common_flag --truth --skip_2p2h 
#diffractive pi0 sample
python gridSelection.py $common_flag --NCDIF --count 5
#extended 2p2h sample
python gridSelection.py $common_flag --ext_2p2h --count 5 --truth
#signal rich sample
#python gridSelection.py $common_flag --bignue --truth --count 5
