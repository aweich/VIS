#!/usr/bin/env Rscript

cran_packages <- c("stringi", "BiocManager", "argparse")
bio_packages <- c("GenomicFeatures", "Gviz", "Rsamtools","biomaRt", "GenomicRanges", "rtracklayer")
for (i in cran_packages) {
  if (!require(i, character.only = TRUE))
    install.packages(i, repos="https://cloud.r-project.org/", quiet = TRUE) #point to the CRAN mirror
}
for (i in bio_packages) {
  if (!require(i, quietly = TRUE, character.only = TRUE))
    invisible(BiocManager::install(i)) 
}

library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(biomaRt)
library(rtracklayer)
library(argparse)

#take input from the command line
parser <- ArgumentParser(description= 'For every BED entry, a chromosomal overview of its location is plotted.')

parser$add_argument('--input_bam', '-ibam', help= 'input BAM file')
parser$add_argument('--input_bed', '-ibed', help= 'input BED file')
parser$add_argument('--outputpath', '-o', help= 'output file', type= 'character')
parser$add_argument('--buffer', '-buffer', help= 'Number of bases that are added to each end of the read window. Good for the format.', type= 'integer')

xargs<- parser$parse_args()


#or use manual data and parameters
bedpath <- xargs$input_bed #"~/Data/VIS_data/GenmoicLocation_100_full_ads.bed"
bampath <- xargs$input_bam #"/home/weichan/Data/VIS_data/BasicMapping_full_ads.bam"
outputpath <- xargs$outputpath #"~/Projects/VIS/PLOTS/"
buffer <- xargs$buffer #50000

bm <- useEnsembl(host = "https://grch37.ensembl.org", #should change hg37 to hg38
                 biomart = "ENSEMBL_MART_ENSEMBL", 
                 dataset = "hsapiens_gene_ensembl")

Cov <- coverage(bampath)
Cov <- GRanges(Cov)
Cov

#data from bed file for the coordinates
bed_data <- read.table(bedpath, header = FALSE, sep = "\t")
bed_data <- bed_data[grep("^chr", bed_data$V1), ] #drop rows that do not start with chr

#create second bed as granges
bed <- makeGRangesFromDataFrame(bed_data, seqnames.field=c("V1"), start.field="V2", end.field=c("V3"), keep.extra.columns = TRUE)

# Extract and print the chromosome, start, and stop coordinates for each line
for (i in 1:nrow(bed_data)) {
  chromosome <- bed_data[i, 1]
  start <- bed_data[i, 2]
  stop <- bed_data[i, 3]


  start_coord <- bed_data[i, 2] - buffer
  end_coord <- bed_data[i, 3] + buffer
  chromosome <- bed_data[i, 1]
  ref_genome <- "hg38"
    
  
  
  #0: Get Gene Region data from biomart
  biomTrack <- BiomartGeneRegionTrack(genome = ref_genome, chromosome = chromosome, 
                                      start = start_coord, end = end_coord,
                                      name = "ENSEMBL", biomart = bm, stacking = "squish")
  
  bed
  #1: Load BED of Insertions
  coord <- AnnotationTrack(bed, name = "Insertions")
  
  #2: Load BAM Coverage
  altrack <- AlignmentsTrack(bampath, isPaired = FALSE) #bampath, isPaired = FALSE
	
  cov <- DataTrack(Cov, type=("heatmap"), name="Coverage")
  #3: Use respective Ideogram and GenomeAxis (automatic)
  itrack <- IdeogramTrack(genome = ref_genome, chromosome = chromosome)
  gtrack <- GenomeAxisTrack()
  
  #4: Plot all
  filename <- sprintf("%s_start%s_end%s.jpeg",chromosome, start, stop)
  jpeg(file=paste(outputpath,filename, sep="/"))
  plotTracks(c(itrack,gtrack,altrack,biomTrack, cov, coord), #altrack, 
             chromosome =chromosome, from = start_coord,
             to = end_coord, transcriptAnnotation="symbol", shape="box", col="black")
  dev.off()
}
