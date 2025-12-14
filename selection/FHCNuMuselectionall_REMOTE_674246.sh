declare -a ncdiffplaylists=("A" "B" "C" "D" "E" "F" "G" "L" "M" "N" "O" "P")
declare -a playlists=("A" "B" "C" "D" "E" "F" "G" "L" "M" "N" "O" "P")

echo "RUNNING NCDIFFSIGNAL JOBS"

for i in "${ncdiffplaylists[@]}"
do
   python gridSelection.py --playlist me1"$i"_p4_NCDiff --ntuple_tag MAD --use-sideband   --truth --cal_POT --mc_only --selection_tag thesis_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_100Universes.tar.gz
done

echo "DONE RUNNING NCDIFFSIGNAL JOBS"
echo "RUNNING AUDIT JOBS"

for i in "${playlists[@]}"
do
   python gridSelection.py --playlist me1"$i"_p4 --ntuple_tag MAD --use-sideband  --cal_POT --mc_only --truth --selection_tag thesis_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_100Universes.tar.gz
   python gridSelection.py --playlist me1"$i"_p4 --ntuple_tag MAD --use-sideband  --cal_POT --data_only --truth --selection_tag thesis_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_100Universes.tar.gz
done

echo "DONE RUNNING AUDIT JOBS"
