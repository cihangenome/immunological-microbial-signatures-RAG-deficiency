#!/bin/bash

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate base
conda activate /data/BCBB_microbiome_db/software/st2/

biomTable=${1%%.*}
metaTable=${2%%.*}
outFolder=${3}


sourcetracker2 -i ${biomTable}.json.biom -m ${metaTable}.txt -o ${outFolder} --jobs 12 \
--source_sink_column SRC_SNK --source_column_value SRC --sink_column_value SNK \
--source_category_column SRC_Type # --source_rarefaction_depth 4000 --sink_rarefaction_depth 4000