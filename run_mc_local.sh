#! /bin/tcsh

cnt=0
gb_folder="mc_code"
full_gb_path="~"

for i in `ls -d $full_gb_path/$gb_folder/*/`;

do
   let cnt=cnt+1
   kk=$gb_folder"_"$(basename $i)
   echo $kk

   matlab -nodisplay -r "GBMC_RUN('$i'); quit"

   echo $i
done

