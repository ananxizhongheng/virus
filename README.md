This directory contains scripts related to the manuscript "Diversity and ecological functions of viruses inhabiting the oil reservoirs".

The scripts (ORIGC analysis pipeline.sh) of metagenomic analysis are placed in "Pipeline" directory. There are two main modules in the pipeline, the construction of the viral catalog and metagenome-assembled genomes. The processes before assembly are same.

######Assembly of reads into contigs
------------------------------------
##1. Trim
The metagenomic raw reads were examined using FastQC v0.11.9 (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), low-quality sequences, primers, and adaptors were trimmed using the Trimmomatic v0.39.
java -jar trimmomatic-0.39.jar PE -threads 8 ./sample_1.fastq.gz ./sample_2.fastq.gz ./sample_1.qc.fq.gz ./sample_unpair_1.qc.fq.gz ./sample_2.qc.fq.gz ./sample_unpair_2.qc.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:20 MINLEN:50

##2. Assembly
Left and right reads (1, and 2) were used as input for assembly with spades and MEGAHIT as follows:

#2.1 For large data (>20G), the trimmed reads were independently assembled using MEGAHIT v1.2.9.
megahit -1 sample _1.qc.fq.gz -2 sample _2.qc.fq.gz -o sample _megahit_output -t 32

#2.2 For small data (<20G), the trimmed reads were independently assembled using metaspades.py.
metaspades.py -k 21,33,55,77,99,127 -1 sample_1.qc.fq.gz -2 sample_2.qc.fq.gz -o sample_spades_output -t 32 -m 2000

####Generation of prokaryotic metagenome-assembled genomes (MAGs)
------------------------------------------------------------------
##1. bin
For each assembly, contigs were binned using the binning module and consolidated into a final bin set using the Bin_refinement module within metaWRAP v1.2.1.
##1.1 metawrap binning
metawrap binning -o sample_INITIAL_BINNING -t 48 -a samplecontigs.fasta --metabat2 --maxbin2 --concoct sample_1.fastq sample_2.fastq

##1.2 metawrap bin_refinement
metawrap bin_refinement -o sample_BIN_REFINEMENT_50_10 -t 48 -A ./sample_INITIAL_BINNING/metabat2_bins/ -B ./sample_INITIAL_BINNING/maxbin2_bins/ -C ./ sample_INITIAL_BINNING/concoct_bins/ -c 50 -x 10

##2. dRep
All produced bin sets were aggregated and de-replicated at 95% average nucleotide identity (ANI) using dRep v3.2.2
dRep dereplicate ./drep95/ -g ./*.fa -p 15 -d -comp 50 -con 10 -nc 0.30 -pa 0.9 -sa 0.95

##3. GTDB
The taxonomy of each MAG was assigned using GTDB-Tk v1.5.0
gtdbtk classify_wf --genome_dir ./00MAG_50_10/ --out_dir ./ --extension fa --prefix bin --cpus 32

####The maximum-likelihood phylogenetic trees of MAGs
-----------------------------------------------------
The maximum-likelihood phylogenetic trees of MAGs were constructed based on a concatenated dataset of 400 universally conserved marker proteins using PhyloPhlAn v3.0.64
phylophlan -i ./ 01drep95_prodigal -d /user/db/phylophlan --diversity high -f my_genome_cell.cfg --accurate -o ./drep95-tree --nproc 15 --min_num_markers 80

###The RPKM values of the MAGs were calculated using CoverM v0.6.1
------------------------------------------------------------------
coverm genome --coupled sample_1.fastq.gz sample_2.fastq.gz -d ./dereplicated_genomes -x fa -t 10 --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 --contig-end-exclusion 0 -m rpkm -o sample_rpkm.txt

### Functional annotation of MAGs
---------------------------------
Open reading frames (ORFs) of these MAGs were predicted with Prodigal v2.6.3
prodigal -a ./MAG_protein_seq.faa -d ./MAG_nucleotide_seq.fasta -o ./MAG_genes.gff -s ./MAG_poteintial.stat -i ./MAG.fna -p meta

The predicted ORFs were annotated using eggNOG-mapper v2.0.1
emapper.py -i ./MAG_protein_seq.faa -o ./MAG_out -m diamond --cpu 32

####Phylogenetic analysis of DsrAB sequences
--------------------------------------------
##1. muscle
The DsrAB sequences were aligned using MUSCLE v3.8 with default parameters.
muscle -align dsrAB.fasta -output ./muscle_dsrAB.aln

##2. trimal
The alignments were then filtered using TrimAL v1.4
trimal -in muscle_dsrAB.aln -out trimal_dsrAB.phy -phylip -cons 50

##3. raxmlHPC
The concatenated DsrAB tree was constructed using RAxML
raxmlHPC-AVX -s trimal_dsrAB.phy -N 100 -n dsrABraxml -f a -p 12345 -x 12345 -m PROTGAMMAIJTT -T 50

###Viral contig identification
------------------------------
Viral contigs were recovered from assembled contigs using VirSorter v2.1 and DeepVirFinder v1.0, Only viral contigs ≥ 10 kb were retained.

##1. DeepVirFinder v1.0
python /home/user/App/DeepVirFinder/dvf.py -i ./samplecontigs.fasta -l 10000 -c 10 -o ./dvf/
The identified viral contigs by DeepVirFinder v1.0 are filtered using R code (00posdvfR.R).

##2. VirSorter v2.1
virsorter run --exclude-lt2gene -i ./samplecontigs.fasta -w ./vs2/sample_virsorter --min-length 10000 --min-score 0.5 -j 10 all

The identified viral contigs from each sample by DeepVirFinder v1.0 and VirSorter v2.1 were aggregated using python script (viral_catalog_02dvf-vs2-for-spades.py, viral_catalog_03dvf-vs2bingji-for-spades.py). The identified viral contigs from all samples were first pooled into a single FASTA file, namely ORIGC_viralcontigs_10K.fasta 

###Viral contig dereplication, virus operational taxonomic unit (vOTU) clustering, and calculation of abundances
----------------------------------------------------------------------------------------------------------------
The identified viral contigs from each sample were clustered into virus operational taxonomic units (vOTUs) using the parameters 95% average nucleotide identity (ANI) and 85% alignment fraction of the smallest scaffolds based on the scripts (https://bitbucket.org/berkeleylab/checkv/src/master/) provided in CheckV v0.8.1
##1. makeblastdb
makeblastdb -in ./ORIGC_viralcontigs_10K.fasta -dbtype nucl -out ./oil_vdb

##2. blastn
blastn -query ./ORIGC_viralcontigs_10K.fasta -db ./oil_vdb -outfmt '6 std qlen slen' -max_target_seqs 10000 -perc_identity 90 -out oil_vblast.tsv -num_threads 4

##3. pairwise ANI
anicalc.py -i oil_vblast.tsv -o oil_vani.tsv
aniclust.py --fna ./ORIGC_viralcontigs_10K.fasta --ani oil_vani.tsv --out oil_vclusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0

##4. Exart vOTUs
cat oil_vclusters.tsv | awk '{print $1}' > vOTU_name.txt
seqkit grep --pattern-file ./vOTU_name.txt ./ORIGC_viralcontigs_10K.fasta > ./vOTU.fasta

##5. Calculation of abundances
coverm contig --coupled ./sample_1.fastq.gz ./sample_2.fastq.gz -r ./vOTU.fasta -t 10 --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 --contig-end-exclusion 0 -m rpkm -o sample_rpkm.txt

####Viral lifestyle
-------------------
Viral lifestyle was predicted by both VIBRANT v1.2.1 and CheckV v0.8.1, while the remaining vOTUs with at least 90% completeness that display no prophage signals or lysogeny-specific genes were considered as potential virulent viruses.
##1. VIBRANT v1.2.1
python3 VIBRANT_run.py -i ./ORIGC_viralcontigs_10K.fasta -t 5

##2. CheckV v0.8.1
checkv end_to_end ORIGC_viralcontigs_10K.fasta ./checkv_output_directory -t 10 -d /media/user/DataBank/database/checkv-db-v1.0

#### Viral taxonomic assignments, viral function annotation, and identification of auxiliary metabolic genes (vAMGs)
--------------------------------------------------------------------------------------------------------------------
##1. Viral taxonomic assignments
To understand the taxonomy of vOTUs, as suggested previously72, we used PhaGCN2.0 and geNomad v1.9 with the ICTV classification to explore the taxonomic affiliation of vOTUs at the family level. The results from these two tools were considered; for a given genome, (1) it was assigned as ‘unclassified’ if both tools failed to assign it, or it was assigned to different taxa, and (2) it was assigned to the taxonomic level determined by one of the tools if the other failed to assign. 
##1. PhaGCN2.0
python run_Speed_up.py --contigs ORIGC_viralcontigs_10K.fasta --len 1700

##2. geNomad v1.9
genomad annotate --cleanup --splits 10 ORIGC_viralcontigs_10K.fasta genomad_output genomad_db

##2. Viral function annotation
Open reading frames (ORFs) of vOTUs were predicted with Prodigal v2.6.3, to understand the function of vOTUs, the predicted viral proteins were first merged and dereplicated using CD-HIT v4.7. The dereplicated viral proteins were assigned to the eggNOG Orthologous Groups database (version 5.0) using eggNOG-mapper v2.0.1.
prodigal -a ./vprotein.fasta -d ./vnucleotide.fasta -i ./ORIGC_viralcontigs_10K.fasta -p meta -g 11 -f gff -q -m
cd-hit -i ./vprotein.fasta -o ./ORIGC_vprotein_90.faa -c 0.90 -s 0.8 -n 5 -M 80000 -g 1 -d 0 -T 32
emapper.py -i ./ ORIGC_vprotein_90.faa -o ./ORIGC_vprotein_emapper -m diamond --cpu 10

##3. Identification of auxiliary metabolic genes (vAMGs)
DRAM-v v1.3.5 were used to recover putative AMGs from vOTUs.
DRAM-v.py annotate -i final-viral-combined-for-dramv.fa -v viral-affi-contigs-for-dramv.tab --threads 20 -o ./annotation
DRAM-v.py distill -i ./annotation/annotations.tsv -o ./annotation/distilled

#### Network analysis
---------------------
Protein-sharing network analysis of vOTUs was performed by vConTACT v.2.0.
vcontact2 --raw-proteins ./refoil_vprotein_10K_vcontact2_newID.fasta --rel-mode Diamond --proteins-fp ./refoil_vprotein_10K_vcontact2_newID_name.csv --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /user/miniconda3/envs/vContact2/bin/cluster_one-1.0.jar --output-dir vcontact2v --threads 72 --db 'ProkaryoticViralRefSeq201-Merged'

#### Virus-host prediction
--------------------------
CRISPR spacer and tRNA were used to predict virus-host interactions.
##1. CRISPR predict
java -cp metaCRT.jar crt ./samplecontigs.fasta ./samplecrispr.txt
Then use R code to filter CRISPR spacer (crisper_spacer_found_cycle_use_viralmeta.R)
##2. tRNA predict
aragorn -t -w -o ./samplevtRNAw.txt ./sample_viralcontigs_10K.fasta
aragorn -t -fons -o ./samplevtRNAfons.txt ./sample_viralcontigs_10K.fasta
Then use python script to merge (tRNAmerge.py).
##3. Fuzznuc
fuzznuc -sequence sample_viralcontigs_10K.fasta -pattern @samplecontigs-crt-clean.txt -outfile ./samplecrisprinv -complement T
fuzznuc -sequence sample_contigs.fasta -pattern @samplevtRNAclean.txt -outfile ./sampletrnaim -complement T
Then use python script to merge (fuzznuc_crisper_trna.py).

