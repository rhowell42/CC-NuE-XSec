#col11 : rhcpt tune cv
#col11.1 : rhc weights alone.
stag="col11.1"
ntag="MAD"
tarball="/pnfs/minerva/resilient/tarballs/hsu-gridselection-macro$(date +%H%M%S).tar.gz"
if [[ -f $tarball ]]; then
    rm  $tarball
fi

common_flag="--selection_tag ${stag} --ntuple_tag ${ntag}  --use_sideband Excess_High_Inline Excess_Low_Inline --tarball ${tarball}"
#common_flag="--selection_tag ${stag} --ntuple_tag ${ntag}  --use_sideband Excess_High_Inline Excess_Low_Inline --tarball ${tarball}"
#common_flag+=" --exclude_universe all --scratch"
echo python gridSelection.py $common_flag
python gridSelection.py $common_flag --truth
#echo python gridSelection.py $common_flag --ext_2p2h --count 10
#python gridSelection.py $common_flag --ext_2p2h --count 10
echo python gridSelection.py $common_flag --NCDIF --count 1
python gridSelection.py $common_flag --NCDIF --count 1
