#!/bin/bash
set -euxo pipefail

if false;then
# gawk -F"\t" -vOFS="\t" 'ARGIND==1{A[$3]=$1} ARGIND==2{if($1 in A){print A[$1]}else{print $1}}' /data2/rawdata2/sample_metadata/tetra_intro/tetra-introgress_200905.txt /data2/rawdata2/tetraintro/201224/repr_dist_AABB/name.txt > WGS_samplelist_l2.txt

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
(while read line;do
  # change from 0.5-1.3 to 0.5-1.5 on 201209
  gawk 'NR>1{if(((NR-1)%5)==0){if(sum/5<0.5 || sum/5>1.5){print 0}else{print 1};sum=$1;i=0}else{sum+=$1;i+=1}} END{if(((NR-1)%5)!=0){if(sum/i<0.5 || sum/i>1.3){print 0}else{print 1}}}' /data2/rawdata2/readDepth/${line}/${CHR}.1M.norm | datamash transpose
done < WGS_samplelist_l2.txt ) | datamash transpose > ${CHR}.CNVfilter.5M.txt
done

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
(while read line;do
  gawk 'NR>1{if($1<1.5&&$1>0.5){print 1}else{print 0}}' /data2/rawdata2/readDepth/${line}/${CHR}.1M.norm | datamash transpose
done < WGS_samplelist_l2.txt ) | datamash transpose > ${CHR}.CNVfilter.1M.txt
done
fi

# gawk -F"\t" -vOFS="\t" 'ARGIND==1{a[$3]=$1} ARGIND==2{print a[$1]}' /data2/rawdata2/sample_metadata/tetra_intro/tetra-introgress_200905.txt /data2/rawdata2/tetraintro/201224/repr_dist_DD/name.txt > WGS_DD_samplelist_l2.txt
for CHR in chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
  ((while read line;do
    gawk 'NR>1{if(((NR-1)%5)==0){if(sum/5<0.5 || sum/5>1.5){print 0}else{print 1};sum=$1;i=0}else{sum+=$1;i+=1}} END{if(((NR-1)%5)!=0){if(sum/i<0.5 || sum/i>1.3){print 0}else{print 1}}}' /data2/rawdata2/readDepth/${line}/${CHR}.1M.norm | datamash transpose
  done < WGS_DD_samplelist_l2.txt ) | datamash transpose > ${CHR}.CNVfilter.5M.txt) &
done
wait
