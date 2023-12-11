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
* For convenience, in the `supplementary/tables/` folder individual supporting tables in CSV format are provided. For a detailed description of the tables, please consult the tables’ headers provided in the `Supplementary_tables.xlsx` file.

## Data
Analyzed data are included in the `data/` directory. The directory contains raw data as well as tables obtained when using scripts required for data analysis and visualization. The directory has the following structure:

- The `Spectra` directory 
The directory contains both MS1 and MS/MS spectra with m/z values for protein spots excised from gels with rice shoot and root proteins. The names of the files correspond to the spots’ numbers as presented in <b>Figure 1a-b</b>. 
  * `trypsin_mass.txt` and `trypsin_mass.tsv` represent mass-to-charge ratios of the trypsin sequence. The respective parental MS1 peaks were then excluded from the analysis; 
  * `*MS1` directory contains lists of m/z peaks for proteins in different protein spots. In case several samples were obtained from one spot, files are marked with underscore symbols;
  * `*MS_MS/MGF` - the MS/MS spectra for each protein spot in the MGF (mascot generic format) format;
  * `*MS_MS/MZID` folder harbors XML files with protein annotations based on MS/MS spectra using the MSGFplus software on the basis of rice protein sequences obtained from the IPG (identical protein groups) database;

- The `IPG_db` directory includes all rice protein sequences downloaded from IPG used as a database for MS1- and MS/MS-based protein annotation utilizing the Mascot server and MSGFplus package.

- The `PQquest_data` directory contains the raw data from the PDQuest software with relative intensities of protein spots and the results of further statistical processing of these data. Hereinafter, the files are split into different organs marked with `roots` and `shoots` prefixes/postfixes. 
  *The files with `_exported.xlsx` and `new_spots` prefixes/postfixes represent raw data with optical density calculated with PDQuest;
  * The tables with `PD_quest` prefix and `_annot.csv` postfix are spot-wise annotations;
  * `PD_quest_shoot_annot_no_rubisco.csv` is the same as the aforementioned table devoid of protein spots containing the large subunit of Rubisco;
* `with_new_log.csv` postfix mark files with logarithmized and normalized optical densities for each protein spot in a gel-wise manner coupled with the respective protein annotations;
  *Tables with the `signif_spots_annot.csv` postfix represent the same information for protein spots showing significant differences between various conditions;
  * `*err_mean_ODS_genes.csv` - merged optical densities for all spots corresponding to a certain protein annotation from all gels grouped according to experimental conditions.

- The `Fasta` directory contains sequences of the identified proteins based on the IPG database.
  * `Protein` - identified protein sequences;
  * `Nucleotide` -  sequences of genes from the reference rice genomes encoding proteins from annotated spots;
  * `Upstream` directory includes 500 b.p. long upstream sequences in the reference rice genome for genes encoding identified proteins.

- The `Clusterization` directory represents the results of the k-means procedure and hierarchical clustering of individual protein spots based on logarithmized and normalized optical densities. 
  * The `heatmaps` directory includes the pair-wise similarity between clusterizations within four spot groups (all, significant, annotated, both annotated and significant, see the Methods section). The hierarchical clustering and k-means procedures are marked with `_hclust_` and `_k-means_ infixes, respectively`;
  * `All_clusterizations_mearged.csv` - the table representing condition-wise optical densities per protein spot assigned to a certain cluster using two clustering methods within four spot groups;
  * `Inter_clusters_similarity.tsv` - similarity between clustering approaches within a certain spot group;
  * `Clusters_intensity.tsv` - mean optical density of all spots from an individual cluster cluster within distinct clustering method and spot group.

- The `Functional_annotation` directory contains the results of the over-representation test of functional annotations belonging to identified proteins.
  *The tables with the `emapper.annotations.tsv` postfix represent functional annotation terms within different systems and ontologies obtained using the eggNOG mapper utility;
  * `all_annots.tsv` - the table with COG (cluster of orthologous genes) codes for proteins grouped according to plant organs.

- The `TF_predictions` directory includes the results of TF binding sites (TFBS) predictions in the upstream sequences of genes encoding identified proteins in the reference rice genomes and the closest homologs in the reference <i>Arabidopsis thaliana</i> reference genome. 
  * `Rice/Shoot` and `/Rice/Root` directories contain putative TF binding sites found with the PlantPAN server. The tables are named according to the accession numbers of the identified proteins; 
  * Tables in the `Rice/Coords` directory represent the coordinate-wise number of sites and genes within a certain TF family for all and significantly different proteins only denoted with `_all.tsv` and `signif.tsv` postfixes, respectively;
  * `Rice/numbers_signif.tsv` and `numbers_all.tsv` - summaries of the number of sites and genes for a certain TF family with genes grouped according to the condition in which the respective protein spots were of the highest optical density; 
  * `AT/All` - predicted number of TFBS for <i>A. thaliana</i> genes encoding homologs of identified rice proteins using the AthaMap resource;
  * `AT/signif` - the same as the aforementioned data for homologs of significantly different rice proteins;
  * `AT/Conditions` - TFBS for <i>A. thaliana</i> genes grouped according to conditions to which the respective rice homologs were attributed;
  * `AT/AT_to_rice_accessions.csv` - correspondence between accession numbers of rice and <i>A. thaliana</i> genes.

- The `For_processing_scripts` directory includes the remaining data used for generating figures in the article.
  * `AA_percent_all_dist.csv` - the cumulative percentage of certain amino acid residues for all identified proteins in roots and soots;
  * `Lysin_coords` prefix encodes files with the coordinates of lysine residues in the sequences of identified proteins enriched with lysines;
  * `proteins_mass_pI.csv` - predicted and experimentally obtained protein mass and isoelectric points for identified proteins;
  * `AT_content_sum.csv` - gene-wise AT content for upstream regions and genes encoding identified proteins;
  * The tables with the `AT_content_per_site` postfix represent the coordinate-wise percentage of AT nucleotides in promotor regions within genes encoding identified proteins;
  * `Coordinates.csv` - the coordinates of genes on the reference rice genome coupled with the condition to which the respective encoded proteins were attributed.

## Scripts
The `scripts/` directory includes all code used for proteome analysis. All required input files are presented in the ` data /` folder. For convenience, the paths to the input files are given relative to the directory with scripts. 
- `IO_lib.py` - the ancillary script with functions for parsing tables and writing files.
- The `filter_trypsin.py` script was used to clean MS/MS spectra from peaks corresponding to trypsin. To run the script use the following command:

  ``` python3 filter_trypsin.py ../data/Spectra/trypsin_mass.txt ../data/Spectra/roots/MS_MS/MGF/1_1.mgf``` 

- `amino_acid_distr.py` - the script for summarizing the distribution of amino acid residues in sequences of identified proteins. The script could be launched with the following command:

  ``` python3 amino_acid_distr.py -f ../data/Fasta/Protein/root_protein_sequences_IPG.fasta``` 
  
- `calculate_clusters_similarity.py` - the script performs a comparison of clustering patterns using different methods and sets of protein spots. The script is run as follows:
    ``` python3 calculate_clusters_similarity.py -c ../data/Clusterization/All_clusterizations_mearged.csv```
  
- `get_251_nt_promotors.py` - the script for splitting 500 b.p. long sequences into two parts to make the resulting files compatible with the CrProm software. The script is launched with the following command:
    ``` python3 get_251_nt_promotors.py -f ../data/Fasta/Upstream/root_nucleotide_sequences_all_500_upstream.fasta```
    
- `get_AT_content.py` - the script for summarizing AT content in the upstream region and genes encoding identified proteins. The script could be launched as follows:

    ``` python3  get_AT_content.py -f ../data/Fasta/Nucleotide/root_nucleotide_sequences_IPG.fasta```

- `build_promotor_graph.py` - the script for summarizing TF binding sites obtained with the PlantPan software. The script is run with the following command:
``` python3 build_promotor_graph.py -r ../data/PQquest_data/roots_with_new_log.csv -s ../data/PQquest_data/shoots_with_new_log.csv -sr ../data/PQquest_data/rice_roots_signif_spots_annot.csv  -ss ../PQquest_data/rice_roots_signif_spots_annot.csv -rd ../data/TF_predictions/Rice/Root -sd ../data/TF_predictions/Rice/Shoot ```

- The `MS_processing.R` script was applied to perform spots’ annotation based on MS/MS spectra with the MSGFplus package;
- `pdquest_process.R` - the code for processing raw data from PDQuest, including logarithmization, normalization, and clusterization of spots and whole proteomes; 
- `Annotation_and_clusterization_plot_genes.r` - the script for visualizing the properties of proteins, TF prediction results, functional enrichments, and similarity between clusters. 


