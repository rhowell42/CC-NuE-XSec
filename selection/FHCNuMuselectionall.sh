declare -a playlists=("A" "B" "C" "D" "E" "F" "G" "L" "M" "N" "O" "P")

echo "RUNNING AUDIT JOBS"

for i in "${playlists[@]}"
do
   python gridSelection.py --playlist me1"$i"_p4 --ntuple_tag MAD --use-sideband  --cal_POT --mc_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20251005.tar.gz
   python gridSelection.py --playlist me1"$i"_p4 --ntuple_tag MAD --use-sideband  --cal_POT --data_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20251005.tar.gz
done

echo "DONE RUNNING AUDIT JOBS"
