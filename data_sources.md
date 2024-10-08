# Data Sources

1. CPTAC data- Proteomic, transcriptomic, and copy number datasets accessed from https://pdc.cancer.gov/pdc/cptac-pancancer (Proteome_BCM_GENCODE_v34_harmonized_v1.zip- all tumour datasets concatenated and Ensembl IDs mapped to gene symbols using HUGO mapping (https://www.genenames.org/download/custom/. Gene-level RNAseq files from tumour datasets and gene-level GISTIC copy number calls for tumor datasets).

2. Paralog pairs from Ensembl 93 (https://ftp.ensembl.org/pub/release-93/) annotated with sequence identity as described in De Kegel et al 2019 (https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008466). Chromosome annotations from HUGO (https://www.genenames.org/download/custom/).

3. DepMap synthetic lethal set and reprocessed CERES scores obtained from De Kegel et al 2021 (https://ars.els-cdn.com/content/image/1-s2.0-S240547122100329X-mmc5.csv and https://ars.els-cdn.com/content/image/1-s2.0-S240547122100329X-mmc3.csv)

4. STRING physical subnetwork accessed from https://stringdb-downloads.org/download/protein.physical.links.v12.0.txt.gz. Mappings from Protein stable IDs (ENSP...) accessed from Ensembl Biomart on 20th June 2024).

5. Annotations of gene loss in DepMap cancer cell lines (A2_lost_essentiality.csv) made as outlined in De Kegel 2021 (https://www.sciencedirect.com/science/article/pii/S240547122100329X)

6. BioGRID database accessed from https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.236/ (Version 4.4.236).

7. CORUM version 4.1 accessed from https://mips.helmholtz-muenchen.de/corum/download

8. EBI Complex Portal accessed from https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv

9. HuMap complexes accessed from http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt

10. Combinatorial synthetic lethal dataset (all_screened_paralog_pairs_25_04_22.csv) generated by processing synthetic lethal data from Parrish et al 2020 (https://linkinghub.elsevier.com/retrieve/pii/S2211-1247(21)01035-4), Gonatopolous-Pournatzis et al 2020 (https://www.nature.com/articles/s41587-020-0437-z), Thompson et al 2021 (https://www.nature.com/articles/s41467-021-21478-9), Dede et al 2020 (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02173-2), and Ito et al 2021 (https://www.nature.com/articles/s41588-021-00967-z) as outlined in methods. To ensure only synthetic lethal pairs, rather than negative genetic interaction pairs, were retained from Ito et al we filtered based on log fold change as performed for the Gonatopoulos-Pournatzis dataset in De Kegel et al 2021. An LFC threshold of -0.8 was used for the Ito et al screen.

11. Ubiquitination site data (Ubiquitination_site_dataset) accessed from Phosphosite (https://www.phosphosite.org/staticDownloads#) on 7th Aug 2024. 

12. Gene ID mapping file (geneidmap_sep24.txt) accessed from HUGO (https://www.genenames.org/download/custom/) on 22nd Sep 2024.
