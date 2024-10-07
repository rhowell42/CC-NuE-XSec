declare -a ncdiffplaylists=("A" "B" "C" "D" "E" "F" "G")
declare -a playlists=("A" "B" "C" "D" "E" "F" "G" "H" "I" "J")

echo "RUNNING NCDIFFSIGNAL JOBS"

python gridSelection.py --playlist me5A_p4_NCDiff --ntuple_tag MAD --use-sideband XSec --truth --no-reco --cal_POT --mc_only --selection_tag xsec --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_2024-03-13-144656.tar.gz --memory 15000

for i in "${ncdiffplaylists[@]}"
do
   python gridSelection.py --playlist me6"$i"_p4_NCDiff --ntuple_tag MAD --use-sideband XSec --truth --no-reco --cal_POT --mc_only --selection_tag xsec --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_2024-03-13-144656.tar.gz --memory 15000
done

echo "RUNNING AUDIT JOBS"

#python gridSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband XSec --truth --no-reco --cal_POT --mc_only --selection_tag xsec --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_2024-03-13-144656.tar.gz --memory 15000

for i in "${playlists[@]}"
do
   python gridSelection.py --playlist me6"$i"_p4 --ntuple_tag MAD --use-sideband XSec --truth --no-reco --cal_POT --mc_only --selection_tag xsec --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_2024-03-13-144656.tar.gz --memory 15000
done

echo "DONE RUNNING AUDIT JOBS"
