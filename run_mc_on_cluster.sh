#! /bin/tcsh

cnt=0
gb_folder="all_Al_gbs"
full_gb_path="~/mc_simulations"

for i in `ls -d $full_gb_path/$gb_folder/*/`;

do
   let cnt=cnt+1
   kk=$gb_folder"_"$(basename $i)

   bsub -n 8 -R span[ptile=8] -W 5700 -o o.$kk -e e.$kk -q single_chassis -J $kk 'matlab -nodisplay -r "GBMC_RUN('\'$i\''); quit"'

   echo $i
done

