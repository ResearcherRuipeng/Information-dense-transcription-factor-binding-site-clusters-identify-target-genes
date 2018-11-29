# Information-dense-transcription-factor-binding-site-clusters-identify-target-genes-with-similar-tiss
This software suite, which was written with Java and object-oriented programming, identifies other genes with similar tissue-wide expression profiles to a given gene using Bray-Curtis Similarity, and transcription factor target genes using knockdown data. It also scans and clusters transcription factor binding sites in gene promoters using information theory-based position weight matrices, then derives machine learning features from the clusters, which are used to predict gene expression profiles and target genes.

Consistent with Figure 1 illustrating the machine learning framework in the main text of the manuscript, the following documentation describes the usage of the programs involved in each of Panels A,B,C,D. The program for scanning TFBSs were described in the publication "Discovery and validation of information theory-based transcription factor and cofactor binding site motifs", the program for the information density-based clustering (IDBC) algorithm was described in the publication "Tandem machine learning for the identification of genes regulated by transcription factors", thus they are included in this repository.

Programs involved in Panel B of Figure 1, which are stored in the 'IdentifyGenesWithSimilarExpressionProfiles' folder:
*************************************************************************************************************
Program: ComputeBrayCurtisSimilarityValuesBasedOnExpressionProfiles.java
Function: This program computes the Bray-Curtis Similarity value between the expression profile of a given protein-coding (PC) gene and every other PC gene.
Parameters:
1. Path of the median RPKM file (i.e. GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct) containing the expression profiles of 56,238 genes that was obtained from GTEx. This file is included in the 'Data' folder.
2. Path of the output file (Format: .csv, format of one line: Ensembl No.,Gene symbol,Similarity value)
3. The name of the given gene

Example command: 
java ComputeBrayCurtisSimilarityValuesBasedOnExpressionProfiles GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct BrayCurtis_NR3C1.csv NR3C1
*************************************************************************************************************



Programs involved in Panel C of Figure 1, which are stored in the 'IdentifyTFTargetsBasedOnCRISPRKnockdownData' folder:
*************************************************************************************************************
Program: GetDEGenesAndNegativesFromCRISPRdata.java
Function: Using the criteria described in the main text, this program obtains protein-coding differentially expressed (DE) genes and negative non-target genes for a given TF from the CRISPR-generated knockdown data from Dixit et al.
Parameters:
1. Path of the .csv file containing the regulatory matrix from the supplementary files of Dixit et al. This file is included in the 'Data' folder.
2. The name of the given TF
3. Path of the file containing TSS coordinates of all transcripts of all known genes in the human genome that was obtained from Ensembl Biomart. This file is included in the 'Data' folder.
4. Path of the output file for differentially expressed (DE) targets (Format: .csv)
Format of one line: Ensembl No.,Gene name,The sum of regulatory coefficients under all guide RNAs of the TF
5. Path of the output file for negative (i.e. non-DE) genes (Format: .csv)
Format of one line: Gene name
6. Value of the threshold on the average fold change in gene expression level (e.g. 1.01, 1.05, 1.1)

Example command:
java GetDEGenesAndNegativesFromCRISPRdata RegulatoryMatrix.csv EGR1 AllGenesFromBiomart.txt DETargets.csv AllNegatives.csv 1.05
*************************************************************************************************************

*************************************************************************************************************
Program: GeneratePromoterIntervalFile.java
Function: This program generates 10kb promoter intervals and outputs them into one single file for all DE targets.
Parameters:
1. Path of the file containing TSS coordinates of all transcripts of all known genes in the human genome that was obtained from Ensembl Biomart.
This file is the 1st parameter of GetDEGenesAndNegativesFromCRISPRdata.java and included in the 'Data' folder.
2. Path of the file containing the list of DE genes, which is one output file (i.e. 4th parameter) of the GetDEGenesAndNegativesFromCRISPRdata.java program.
3. Path of the output file (Format: .bed)
An example line: chr2	84895624	84905624	TMSB10	ENSG00000034510	0.279712119	protein_coding	1
This file will be intersected with ChIP-seq peaks of the TF from the K562 cell line, which results in an output file containing a subset of promoter intervals.

Example command:
java GeneratePromoterIntervalFile AllGenesFromBiomart.txt DETargets.csv PromoterIntervalFile.bed
*************************************************************************************************************

*************************************************************************************************************
Program: GeneratePositives.java
Function: This program generates the list of all positives (i.e. direct DE targets) for a TF.
Parameters:
1. Path of the input file, which is the output file of intersecting the output file (i.e. 3rd parameter) of GeneratePromoterIntervalFile.java with K562 ChIP-seq peaks of the TF. (Format: .bed)
2. Path of the output file (Format:.csv)
Format of one line: The sum of regulatory coefficients under all guide RNAs of the TF,Gene name,Ensembl No.

Example command:
java GeneratePositives PromoterIntervalIntersectedWithPeaksFile.bed Positives.csv
*************************************************************************************************************

*************************************************************************************************************
Program: SelectNegatives.java
Function: This program selects a subset of negatives that have the smallest similarity values to the positive with the greatest average coefficient in terms of expression profiles.
Parameters:
1. Path of the median RPKM file (i.e. GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct) containing the expression profiles of 56,238 genes that was obtained from GTEx
This file is the 1st parameter of ComputeBrayCurtisSimilarityValuesBasedOnExpressionProfiles.java, and included in the 'Data' subfolder of the current main folder.
2. Name of the positive that has the greatest average coefficient (i.e. greatest average fold change), which is obtained from the output file of GeneratePositives.java
3. Path of the file containing all negative genes.
This file is an output file (i.e. 5th parameter) of GetDEGenesAndNegativesFromCRISPRdata.java.
4. Path of the output file (Format: .csv)
Format of one line: Ensembl No.,Gene name,Bray-Curtis similarity value
5. Path of the file containing TSS coordinates of all transcripts of all known genes in the human genome that was obtained from Ensembl Biomart.
This file is the 1st parameter of GetDEGenesAndNegativesFromCRISPRdata.java and included in the 'Data' folder.
6. A number of Integer type: the number of positives, which is obtained from the output file of GenerateTruePositives.java

Example command:
java SelectNegatives GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct TMSB10 AllNegatives.csv Negatives.csv AllGenesFromBiomart.txt 16
*************************************************************************************************************



Programs involved in Panel D of Figure 1, which are stored in the 'IdentifyTFTargetsBasedOnsiRNAKnockdownData' folder:
*************************************************************************************************************
Program: GetPositivesAndNegativesFromsiRNAData.java
Function: Using the criteria described in the main text, this program obtains positive genes (i.e. targets) and negative genes (i.e. non-targets) for a given TF from the siRNA-generated knockdown data from Cusanovich et al.
Parameters:
1. Path of the file containing P-values of 8,872 genes and the overlapping information (i.e. numbers of peaks) of these genes' promoters and ChIP-seq peaks from the GM12878 cell line, which is a supplementary file from Cusanovich et al.
This file is the 1st parameter of GeneratePromoterIntervalFile.java and included in the 'Data' folder.
2. Name of the given TF
3. Path of the output file for positive genes. Each line contains a positive gene. (Format: .csv)
Format of one line: Ensembl No.,Gene name,P value,number of GM12878 ChIP-seq peaks of the TF overlapping the 10kb promoters of this gene
4. Path of the output file for negative genes. Each line contains a negative gene. (Format: .csv)
Format of one line: Ensembl No.,Gene name,P value
5. Path of the file containing transcription start site (TSS) coordinates of all transcripts of all known genes in the human genome that was obtained from Ensembl Biomart
This file is the 1st parameter of GetDEGenesAndNegativesFromCRISPRdata.java and included in the 'Data' folder.

Example command:
java GetPositivesAndNegativesFromsiRNAData PValueMatrix.txt BATF Positives.csv Negatives.csv AllGenesFromBiomart.txt
*************************************************************************************************************



Programs involved in Panel A of Figure 1
*************************************************************************************************************
Program: GeneratePromoterIntervals.java
Function: This program generates the merged intevals of 10kb promoters upstream of all TSSs for each gene in the list of positives or negatives.
Parameters:
1. Path of the file containing transcription start site (TSS) coordinates of all transcripts of all known genes in the human genome that was obtained from Ensembl Biomart
This file is the 1st parameter of GetDEGenesAndNegativesFromCRISPRdata.java and included in the 'Data' folder.
(Format of one line: Ensembl Gene ID,Ensembl Transcript ID,Chromosome Name,Transcription Start Site (TSS),Associated Gene Name,Associated Transcript Name,Strand,Transcript Start (bp),Transcript End (bp),Transcript length (including UTRs and CDS),Gene Start (bp),Gene End (bp),Gene type,Transcript type)
2. Path of the file containing the list of positives/negatives
(This file is obtained from, and of the same format as the output file of ComputeBrayCurtisSimilarityValuesBasedOnExpressionProfiles.java.)
3. Path of the output folder
(Each merged promoter interval of one gene is stored in a separate file of the BED format. 
If Gene A has two merged promoter intervals, then its output files are named as A_0.txt and A_1.txt
Example format in an output file: chrY\t22428939\t22438939\tAC008175.1\tENSG00000235059\tlincRNA\t1)

Example command:
java GeneratePromoterIntervals AllGenesFromBiomart.txt Positives.csv PromoterIntervalPositiveFolder
*************************************************************************************************************

*************************************************************************************************************
Program: fastaFromBed
Function: This program converts a promoter interval of the BED format to a DNA sequence.
Parameters:
1. Path of the hg38 genome assembly file of the FASTA format
2. Path of a promoter interval file that is the output of GeneratePromoterIntervals.java
2. Path of the output file containing the DNA sequence of the promoter interval

Example command:
./fastaFromBed -fi hg38.fasta -bed PromoterIntervalFile.txt -fo PromoterSequence.txt
*************************************************************************************************************

*************************************************************************************************************
Program: IntersectPromotersWithDHSs
Function: This program intersects each promoter interval in a folder with one or multiple DNase I hypersensitive files by filtering for binding sites within DHSs.
Paramenter:
1. Path of the folder containing binding site files, which is the output folder of the ScanPromoterSequencesWithiPWMs program
As stated above in the description of the ScanPromoterSequencesWithiPWMs program, a binding site file contains all sites in a promoter interval.
2. Path of the folder containing promoter interval files, which is the output folder of GeneratePromoterIntervals.java
3. Path of the output folder in which each file contains all binding sites within DHSs in a promoter interval.
4. Path of one or multiple DNase I hypersensitive files. Each file's path is a separate parameter. 
Each file is of BED format (format: Chromosome\tStarting coordinate\tEnding coordinate), and each line is a DHS (example: chr1\t9980\t10410).
All the files and their sources used in this machine learning framework are described in the main text.

Example command:
java Main BindingSitePositiveFolder PromoterIntervalPositiveFolder BindingSitePositiveDNaseFolder wgEncodeAwgDnaseDukeGm19238UniPk.narrowPeak
*************************************************************************************************************

*************************************************************************************************************
Program: ExtractFeaturesFromClusters
Function: In identifying genes with similar tissue-wide expression profiles, this program extracts the 7 features from the cluster files for each TSS of each positive/negative; in predicting direct DE (differentially expressed) TF targets, it extracts the 6 features.
Parameters:
1. Path of the folder containing cluster files, which is the output folder of the IDBC program
2. Path of the file containing TSS coordinates of all transcripts of all known genes in the human genome that was obtained from Ensembl Biomart
This file is the 1st parameter of GetDEGenesAndNegativesFromCRISPRdata.java and included in the 'Data' folder.
3. Path of the folder containing promoter interval files, which is the output folder of GeneratePromoterIntervals.java
4. Path of the folder storing iPWMs used to detect binding sites
This folder is the 2nd parameter of the ScanPromoterSequencesWithiPWMs program.
5. A number of Double type
(This number * Rsequence of each iPWM) is the minimum threshold on Ri values of binding sites for determining sites taken into account in the following two features.
The number of binding sites of each TF within this cluster;
The sum of Ri values of binding sites of each TF within this cluster
6. A number of Double type
(This number * Rsequence of each iPWM) is the minimum threshold on Ri values of binding sites for determining strong sites taken into account in the following two features.
The number of strong binding sites (Ri > Rsequence) of each TF within this cluster;
The sum of Ri values of strong binding sites (Ri > Rsequence) of each TF within this cluster
7. Path of the output folder in which each positive/negative has a feature set file. This feature set file contains the feature sets consisting of 7/6 features for all TSSs of this gene.
8. Boolean type: true / false
This boolean parameter indicates that the binding sites in cluster files use absolute coordinates in the chromosome or relative coordinates in the promoter interval.
true: Absolute coordinates are used. In this machine learning framework, after intersecting promoter intervals with DNase I Hypersensitive Sites (DHSs), the binding sites in cluster files used absolute coordinates.
false: Relative coordinates are used. In this machine learning framework, if promoter intervals were not intersected with (DHSs), the binding sites in cluster files used absolute coordinates.
9. Boolean type: true / false
This boolean parameter indicates that the binding sites in cluster files were detected by iPWMs of multiple TFs or one iPWM of one single TF.
false: iPWMs of multiple TFs were used. In predicting genes with similar expression profiles to a given gene, 94 iPWMs of 82 TFs were used to detect binding sites, and 7 features were extracted from each cluster.
true: One iPWM of one single TF was used. In predicting direct TF targets, one single iPWM of the TF was used to detect binding sites, and 6 features were extracted from each cluster.

Example command:
java Main ClusterFolder/BindingSitePositiveFolder_d25_i939.0_x0.1 AllGenesFromBiomart.txt PromoterIntervalPositiveFolder iPWMs 0.1 1 FeatureSetPositiveFolder false false
*************************************************************************************************************

*************************************************************************************************************
Program: GenerateInputFilesForClassifiers
Function: This program generates an input file for machine learning classifiers in MATLAB, containing feature vectors of all positives and negatives for a TF
Parameters:
1. Path of the folder containing the feature set files of all positives and negatives
This folder is the union of the two output folders generated by two runs of the ExtractFeaturesFromClusters program on positives and negatives.
2. Path of the file containing the list of positives
This parameter is the 2nd parameter of GeneratePromoterIntervals.java.
3. Path of the file containing the list of negatives
This parameter is the 2nd parameter of GeneratePromoterIntervals.java.
4. Path of the output file, which is the input file for  machine learning classifiers in MATLAB (.csv)
5, oolean type: true / false
This boolean parameter indicates that the binding sites in cluster files were detected by iPWMs of multiple TFs or one iPWM of one single TF.
false: iPWMs of multiple TFs were used. In predicting genes with similar expression profiles to a given gene, 94 iPWMs of 82 TFs were used to detect binding sites, and 7 features were extracted from each cluster.
true: One iPWM of one single TF was used. In predicting direct TF targets, one single iPWM of the TF was used to detect binding sites, and 6 features were extracted from each cluster.

Example command:
java Main FeatureSetFolder Positives.csv Negatives.csv OutputFile.csv true
*************************************************************************************************************
