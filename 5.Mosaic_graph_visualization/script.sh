#!/bib/bash
set -euxo pipefail

WD=/data2/rawdata2/tetraintro/

# plot the mosai graph for each chromosomes
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
  plotH=$(gawk 'NR==1{print (NF-3)/5}' ${CHR}.grp)
  Rscript ${WD}/bin/draw_haplo_block.R -H ${plotH} -m /data2/rawdata2/sample_metadata/tetra_intro/metadata_201221.txt -i ${CHR}.contri.grp -o ${CHR}.pdf -s ${CHR}.contri.fi -n 20 -c ${WD}/bin/20_distinct_colors2.txt
done
