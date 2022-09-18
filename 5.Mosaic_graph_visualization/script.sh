#!/bib/bash

source ../parameters.sh

while read CHR;do
  plotH=$(gawk 'NR==1{print (NF-3)/5}' ../04-Bayesian-smoothing/${CHR}.smoothed.grp)
  Rscript ${script_dir}/5.Mosaic_graph_visualization/draw_AHG.R -b ${bin_size} -H ${plotH} -i ../04-Bayesian-smoothing/${CHR}.smoothed.grp -o ${CHR}.pdf -s ../03-Ancestry-inference/${CHR}.fi -n 20 -c ${script_dir}/5.Mosaic_graph_visualization/20_distinct_colors.txt
done < ../chr_list.txt
