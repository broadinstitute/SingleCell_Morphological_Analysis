#!/usr/bin/env bash

plates="A B C D E F H J K L M N O"

channels="DNA Mito Phalloidin WGA ER"

for w in 1 2 3 4 5 6;
do

for p in $plates;
do

for c in $channels;
do

   echo python3 ~/imaging-backup-scripts/restore_intelligent.py pooled-cell-painting "projects/2018_11_20_Periscope_Calico/20210422_6W_CP257/images_corrected_cropped/CP257${p}_Well$w/Corr$c"
   python3 ~/imaging-backup-scripts/restore_intelligent.py pooled-cell-painting "projects/2018_11_20_Periscope_Calico/20210422_6W_CP257/images_corrected_cropped/CP257${p}_Well$w/Corr$c"
   
done;
done;
done;