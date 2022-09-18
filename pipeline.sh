#!/bin/bash
set -euxo pipefail

source ./parameters.sh

# prepare basic files
cd ${working_dir}
bash ${script_dir}/1.Binwise_genetic_distance/prep_basic_file.sh

# calculate pairwise dist
cd ${working_dir}
mkdir -p 01-Binwise-genetic-distance
cd 01-Binwise-genetic-distance
bash ${script_dir}/1.Binwise_genetic_distance/calculate_dist.sh

# detect CNV bins
cd ${working_dir}
mkdir -p cnv_masker
cd cnv_masker
bash ${script_dir}/1.Binwise_genetic_distance/cnv_masker.sh

cd ${working_dir}
mkdir -p 02-Initial-grouping
cd 02-Initial-grouping
bash ${script_dir}/2.Initial_grouping/script.sh

cd ${working_dir}
mkdir -p 03-Ancestry-inference
cd 03-Ancestry-inference
bash ${script_dir}/3.Ancestry_inference/script.sh

cd ${working_dir}
mkdir -p 04-Bayesian-smoothing
cd 04-Bayesian-smoothing
bash ${script_dir}/4.Bayesian_smoothing/script.sh

cd ${working_dir}
mkdir -p 05-Mosaic-graph-visualization
cd 05-Mosaic-graph-visualization
bash ${script_dir}/5.Mosaic_graph_visualization/script.sh
