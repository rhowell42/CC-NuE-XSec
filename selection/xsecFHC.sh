declare -a playlists=("B" "C" "D" "E" "F" "G")

echo "RUNNING AUDIT JOBS"

for i in "${playlists[@]}"
do
   python gridSelection.py --playlist me1"$i"_p4 --ntuple_tag MAD --use-sideband XSec --cal_POT --mc_only --truth --no-reco --selection_tag xsec --tarball  /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_2024-03-27-144104.tar.gz --memory 10000
done

echo "DONE RUNNING AUDIT JOBS"
