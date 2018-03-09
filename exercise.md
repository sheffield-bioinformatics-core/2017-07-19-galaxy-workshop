---
output:
  pdf_document: default
  html_document: default
---
<style> @import url('style.css'); </style>


<img src="media/TUOS_PRIMARY_LOGO_FULL COLOUR.png" height=100px align=right>

# RNA-Seq Extra Practice

## Preamble

> In the study of Brooks et al. 2011, the Pasilla (PS) gene, Drosophila homologue of the Human splicing regulators Nova-1 and Nova-2 Proteins, was depleted in Drosophila melanogaster by RNAi. The authors wanted to identify exons that are regulated by Pasilla gene using RNA sequencing data.
Total RNA was isolated and used for preparing either single-end or paired-end RNA-seq libraries for treated (PS depleted) samples and untreated samples. These libraries were sequenced to obtain a collection of RNA sequencing reads for each sample. The effects of Pasilla gene depletion on splicing events can then be analyzed by comparison of RNA sequencing data of the treated (PS depleted) and the untreated samples.
The genome of Drosophila melanogaster is known and assembled. It can be used as reference genome to ease this analysis. In a reference based RNA-seq data analysis, the reads are aligned (or mapped) against a reference genome, Drosophila melanogaster here, to significantly improve the ability to reconstruct transcripts and then identify differences of expression between several conditions.

## Exercise

To save time and memory, pre-counted files are available online ([https://zenodo.org/record/1185122#.WqGfbnW0PCI](https://zenodo.org/record/1185122#.WqGfbnW0PCI))

- GSM461176_untreat_single.counts
- GSM461177_untreat_paired.counts
- GSM461178_untreat_paired.counts
- GSM461179_treat_single.counts
- GSM461180_treat_paired.counts
- GSM461181_treat_paired.counts
- GSM461182_untreat_single.counts


Your exercise is to identify differentially-expressed genes in the dataset, and which pathways are enriched among the results. Please send a short (1 page) report of your analysis of the dataset including

- a brief description of the methods you used
- the number of differentially-expressed genes with a p-value less than 0.05 and absolute fold-change greater than 1
- how many enriched pathways you found in total, and the names of the five most-enriched pathways

## Optional (in your own time)

The online archive for this dataset also includes the fastq files for each sample.

[https://zenodo.org/record/1185122#.WqGXzHW0PCI](https://zenodo.org/record/1185122#.WqGXzHW0PCI)

You can use these files to practice the QC, alignment and quantification steps used in the workshop. However, be aware that these operations may take a while to run on the public galaxy server.

