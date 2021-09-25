stag="col9"
ntag="Had"
tarball="/pnfs/minerva/resilient/tarballs/hsu-gridselection-macro.tar.gz"
if [[ -f $tarball ]]; then
    rm  $tarball
fi

common_flag="--selection_tag ${stag} --ntuple_tag ${ntag}  --use_sideband Excess_High_Inline Excess_Low_Inline --tarball ${tarball}"
echo python gridSelection.py $common_flag
python gridSelection.py $common_flag --truth
#echo python gridSelection.py $common_flag --ext_2p2h --count 10
#python gridSelection.py $common_flag --ext_2p2h --count 10
echo python gridSelection.py $common_flag --NCDIF --count 10
python gridSelection.py $common_flag --NCDIF --count 10
