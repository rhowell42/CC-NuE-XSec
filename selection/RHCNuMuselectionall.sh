declare -a playlists=("A" "B" "C" "D" "E" "F" "G" "H" "I" "J")

echo "RUNNING AUDIT JOBS"

python gridSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --mc_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20251012.tar.gz --scratch
python gridSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --data_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20251012.tar.gz --scratch 

for i in "${playlists[@]}"
do
   python gridSelection.py --playlist me6"$i"_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --mc_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20251012.tar.gz --scratch 
   python gridSelection.py --playlist me6"$i"_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --data_only --selection_tag paper_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-EventSelection_1000Univ_20251012.tar.gz --scratch 
=======
declare -a ncdiffplaylists=("A" "B" "C" "D" "E" "F" "G")
declare -a playlists=("A" "B" "C" "D" "E" "F" "G" "H" "I" "J")

echo "RUNNING NCDIFFSIGNAL JOBS"

python gridSelection.py --playlist me5A_p4_NCDiff --ntuple_tag MAD --use-sideband  --truth --cal_POT --mc_only --selection_tag thesis_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_100Universes.tar.gz 

for i in "${ncdiffplaylists[@]}"
do
   python gridSelection.py --playlist me6"$i"_p4_NCDiff --ntuple_tag MAD --use-sideband  --truth --cal_POT --mc_only --selection_tag thesis_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_100Universes.tar.gz 
done

echo "DONE RUNNING NCDIFFSIGNAL JOBS"
echo "RUNNING AUDIT JOBS"

python gridSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --mc_only --selection_tag thesis_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_100Universes.tar.gz 
python gridSelection.py --playlist me5A_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --data_only --selection_tag thesis_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_100Universes.tar.gz 

for i in "${playlists[@]}"
do
   python gridSelection.py --playlist me6"$i"_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --mc_only --selection_tag thesis_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_100Universes.tar.gz 
   python gridSelection.py --playlist me6"$i"_p4 --ntuple_tag MAD --use-sideband  --truth --cal_POT --data_only --selection_tag thesis_muon --tarball /pnfs/minerva/resilient/tarballs/rhowell-CCNUE_selection_100Universes.tar.gz 
>>>>>>> feature/sterile_neutrino
done

echo "DONE RUNNING AUDIT JOBS"
