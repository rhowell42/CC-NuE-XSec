declare -a ncdiffplaylists=("A" "B" "C" "D" "E" "F" "G" "L" "M" "N" "O" "P")
declare -a playlists=("A" "B" "C" "D" "E" "F" "G" "L" "M" "N" "O" "P")

echo "RUNNING NCDIFFSIGNAL JOBS"

for i in "${ncdiffplaylists[@]}"
do
   python gridSelection.py --playlist me1"$i"_p4_NCDiff --ntuple_tag MAD --use-sideband  dEdX --cal_POT --mc_only --selection_tag paper --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20250605.tar.gz 
done

echo "DONE RUNNING NCDIFFSIGNAL JOBS"
echo "RUNNING AUDIT JOBS"

for i in "${playlists[@]}"
do
   python gridSelection.py --playlist me1"$i"_p4 --ntuple_tag MAD --use-sideband dEdX --cal_POT --mc_only --selection_tag paper --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20250605.tar.gz
   python gridSelection.py --playlist me1"$i"_p4_BigGenie --ntuple_tag MAD --use-sideband dEdX --cal_POT --mc_only --selection_tag paper --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20250605.tar.gz 
   python gridSelection.py --playlist me1"$i"_p4 --ntuple_tag MAD --use-sideband dEdX --cal_POT --data_only --selection_tag paper --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20250605.tar.gz 
done

echo "DONE RUNNING AUDIT JOBS"
