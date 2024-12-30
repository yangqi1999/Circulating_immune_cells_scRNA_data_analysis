eggnog_annotations="/dellfsqd2/ST_OCEAN/USER/liyanan2/03_job/11_anno/05_qsm/eggnog2/eggnog-result.emapper.annotations"
kegg_info="/dellfsqd2/ST_OCEAN/USER/liyanan2/03_job/single_cell_seq_analysis/make.OrgDb/kegg_info.RData"
tax_id="682880"
genus="Lampetra"
species="morii"
sed '/^##/d' $eggnog_annotations |sed 's/#//' >eggnog.anno
/dellfsqd2/ST_OCEAN/USER/liyanan2/01_Software/conda/envs/R-4.0/bin/Rscript /dellfsqd2/ST_OCEAN/USER/liyanan2/03_job/single_cell_seq_analysis/make.OrgDb/make.OrgDb.r eggnog.anno  $kegg_info $tax_id $genus $species
