#!/bin/bash
tag=Had
for i in me1A me1B me1C me1D me1E me1F me1G me1L me1M me1N me1O me1P 
do
#data 
  find /pnfs/minerva/persistent/users/hsu/NuECCQE-Had-data/${i}/grid/minerva/ana/ -name *.root | pnfs2xrootd.sh > playlist_ccqenue_data${i}_nx_${tag}.txt
#mc
  find /pnfs/minerva/persistent/users/hsu/NuECCQE-Had-mc/${i}_merged/ -name *.root | pnfs2xrootd.sh > playlist_ccqenue_mc${i}_nx_${tag}.txt
#ncdif
  find /pnfs/minerva/persistent/users/hsu/NuECCQE-Had-mc/${i}-NCDIF/grid/NCDiffSignal/minerva/ana/ -name *.root | pnfs2xrootd.sh > playlist_ccqenue_mc${i}_NCDIF_${tag}.txt
done
