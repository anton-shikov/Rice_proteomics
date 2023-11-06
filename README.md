# Rice proteomics
The repository with working codes and supplementary materials for the study of rice proteomes in shoots and roots under anoxia and post-anoxic re-aeration using two-dimensional difference gel electrophoresis (2D-DIGE). 

## Contents 

This repository contains scripts used for statistical processing, MS-based protein identification, visualization of the results, and raw data. Please consult the methods section of the article for extra details:

## Figures
Figures from the main text are available in the `pics/` directory. For the description, please consult the Results section of the article.

## Supplementary
All supplementary material is located in the `supplementary/` directory. 
`Supplementary_figures.docx` is a Microsoft Word document with a description of all supplementary figures.
`Supplementary_tables.xlsx` is an Exel table with all supplementary tables.

* The `supplementary/pics/` folder contains individual supplementary figures.
* For convenience, in the `supplementary/tables/` folder individual supporting tables in CSV format are provided. For a detailed description of the tables, please consult the `Supplementary_tables.xlsx` file.

## Data
Analyzed data are included in the `data/` directory. The directory contains raw data as well as tables obtained when using scripts required for data analysis and vizualization. The directory has the following structure:

- The `/Spectra` directory 
  * `trypsin_mass.txt` and `trypsin_mass.tsv`;
  * `*/MS1` ;
  * `*/MS_MS/MGF` ;
  * `*/MS_MS/MZID` ;
- The `/IPG_db` directory
  
- The `/PQquest_data` directory
  * `_exported.xlsx` and `new_spots`;
  * `PD_quest` and `_annot.csv`;
  * `PD_quest_shoot_annot_no_rubisco.csv` ;
  * `signif_spots_annot.csv` ;
  * `with_new_log.csv` ;
  * `*err_mean_ODS_genes.csv` ;
- The `/Fasta` directory
  * `heatmaps` and `_hclust_` and `_k-means_`;
  * `Protein` ;
  * `Upstream` ;
- The `/Clusterization` directory
  * `Nucleotide` and `new_spots`;
  * `All_clusterizations_mearged.csv` ;
  * `Inter_clusters_similarity.tsv` ;
  * `Clusters_intensity.tsv` ;
- The `/Functional_annotation` directory
  * `emapper.annotations.tsv`;
  * `all_annots.tsv` ;
- The `/TF_predictions` directory
  * `/Rice/Shoot` and `/Rice/Root`;
  * `/Rice/Coords` and `_all.tsv` and `signif.tsv`;
  * `/Rice/numbers_signif.tsv` and `numbers_all.tsv`; 
  * `/AT/All` `all_sites.csv` `common_sites.csv` `merged`;
  * `/AT/signif` `All_`;
  * `/AT/Conditions` ;
  * `/AT/AT_to_rice_accessions.csv` ;
- The `/For_processing_scripts` directory 
  * `AA_percent_all_dist.csv` ;
  * `Lysin_coords` ;
  * `extr_amino_acid_content.csv` ;
  * `proteins_mass_pI.csv` ;
  * `AT_content_sum.csv` ;
  * `AT_content_per_site` ;
  * `Coordinates.csv`


## Scripts
The `scripts/` directory includes all code used for proteome analysis.

- `IO_lib.py` 
- `filter_trypsin.py` - To run the script use the following command:

  ``` python3 filter_trypsin.py ../data/Spectra/trypsin_mass.txt ../data/Spectra/roots/MS_MS/MGF/1_1.mgf``` 

- `amino_acid_distr.py` - To run the script use the following command:

  ``` python3 amino_acid_distr.py -f ../data/Fasta/Protein/root_protein_sequences_IPG.fasta``` 
  
- `calculate_clusters_similarity.py` - To run the script use the following command:
    ``` python3 calculate_clusters_similarity.py -c ../data/Clusterization/All_clusterizations_mearged.csv```
  
- `get_251_nt_promotors.py` - To run the script use the following command:

    ``` python3 get_251_nt_promotors.py -f ../data/Fasta/Upstream/root_nucleotide_sequences_all_500_upstream.fasta```
    
- `get_AT_content.py` - To run the script use the following command:

    ``` python3  get_AT_content.py -f ../data/Fasta/Nucleotide/root_nucleotide_sequences_IPG.fasta```

- `build_promotor_graph.py` - To run the script use the following command:

``` python3 build_promotor_graph.py -r ../data/PQquest_data/roots_with_new_log.csv -s ../data/PQquest_data/shoots_with_new_log.csv -sr ../data/PQquest_data/rice_roots_signif_spots_annot.csv  -ss ../PQquest_data/rice_roots_signif_spots_annot.csv -rd ../data/TF_predictions/Rice/Root -sd ../data/TF_predictions/Rice/Shoot ```

- `MS_processing.R` 
- `pdquest_process.R` 
- `Annotation_and_clusterization_plot_genes.r` 

