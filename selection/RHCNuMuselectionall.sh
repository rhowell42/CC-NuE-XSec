declare -a playlists=("A" "B" "C" "D" "E" "F" "G" "H" "I" "J")

echo "RUNNING AUDIT JOBS"

python gridSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --mc_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20250218.tar.gz 
python gridSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --data_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20250218.tar.gz 

for i in "${playlists[@]}"
do
   python gridSelection.py --playlist me6"$i"_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --mc_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20250218.tar.gz 
   python gridSelection.py --playlist me6"$i"_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --data_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20250218.tar.gz 
done

echo "DONE RUNNING AUDIT JOBS"
