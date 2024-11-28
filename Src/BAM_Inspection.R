#!/usr/bin/env Rscript
gc()
getwd()
# Create personal library folder if it doesn't exist
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)

# Add the personal library folder to the library path
.libPaths(Sys.getenv("R_LIBS_USER"))

# Function to install packages if missing
install_if_missing <- function(packages, bio = FALSE) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (bio) {
        invisible(BiocManager::install(pkg, quiet = TRUE))
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
      }
      library(pkg, character.only = TRUE)
    }
  }
}

# Define package lists
cran_packages <- c("stringi", "BiocManager", "argparse")
bio_packages <- c("GenomicFeatures", "Gviz", "Rsamtools", "biomaRt", "GenomicRanges", "rtracklayer", "Homo.sapiens")

# Install and load packages
install_if_missing(cran_packages)
install_if_missing(bio_packages, bio = TRUE)
library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(biomaRt)
library(rtracklayer)
library(argparse)
library(Homo.sapiens)

#take input from the command line
parser <- ArgumentParser(description= 'For every BED entry, a chromosomal overview of its location is plotted.')

parser$add_argument('--input_bam', '-ibam', help= 'input BAM file')
parser$add_argument('--input_bed', '-ibed', help= 'input BED file')
parser$add_argument('--input_H3K4Me1', '-iH3K4Me1', help= 'input H3K4Me1 file')
parser$add_argument('--input_H3K4Me3', '-iH3K4Me3', help= 'input H3K4Me3 file')
parser$add_argument('--input_H3K27Ac', '-iH3K27Ac', help= 'input H3K27Ac file')
parser$add_argument('--input_gtf', '-igtf', help= 'input gtf file')
parser$add_argument('--input_TF', '-iTF', help= 'input TF')
parser$add_argument('--outputpath', '-o', help= 'output file', type= 'character')
parser$add_argument('--buffer', '-buffer', help= 'Number of bases that are added to each end of the read window. Good for the format.', type= 'integer')

xargs<- parser$parse_args()

options(ucscChromosomeNames=FALSE)

#or use manual data and parameters
bedpath <- xargs$input_bed #"~/Data/VIS_data/GenmoicLocation_100_full_ads.bed"
bampath <- xargs$input_bam #"/home/weichan/Data/VIS_data/BasicMapping_full_ads.bam"
H3K4Me1path <- xargs$input_H3K4Me1 #"~/Data/VIS_data/GenmoicLocation_100_full_ads.bed"
H3K4Me3path <- xargs$input_H3K4Me3 #"~/Data/VIS_data/GenmoicLocation_100_full_ads.bed"
H3K27Acpath <- xargs$input_H3K27Ac #"~/Data/VIS_data/GenmoicLocation_100_full_ads.bed"
gtfpath <- xargs$input_gtf
TFpath <- xargs$input_TF
outputpath <- xargs$outputpath #"~/Projects/VIS/PLOTS/"
buffer <- xargs$buffer #50000

bm <- useEnsembl(host = "https://grch37.ensembl.org", #should change hg37 to hg38
                 biomart = "ENSEMBL_MART_ENSEMBL", 
                 dataset = "hsapiens_gene_ensembl")

# Create a TxDb object from the GTF file
txdb <- makeTxDbFromGFF(gtfpath)
columns(txdb)
head(keys(txdb))
#functional genomics
H3K4Me1 <- read.table(H3K4Me1path, header = FALSE, sep = "\t",skip = 1)
H3K4Me1 <- makeGRangesFromDataFrame(H3K4Me1, seqnames.field=c("V1"), start.field="V2", end.field=c("V3"), keep.extra.columns = T)

H3K27Ac <- read.table(H3K27Acpath, header = FALSE, sep = "\t",skip = 1)
H3K27Ac <- makeGRangesFromDataFrame(H3K27Ac, seqnames.field=c("V1"), start.field="V2", end.field=c("V3"), keep.extra.columns = TRUE)

H3K4Me3 <- read.table(H3K4Me3path, header = FALSE, sep = "\t",skip = 1)
H3K4Me3 <- makeGRangesFromDataFrame(H3K4Me3, seqnames.field=c("V1"), start.field="V2", end.field=c("V3"), keep.extra.columns = TRUE)

#DNaseH <- read.table(DNaseHpath, header = F, skip=1)
#DNaseH <- makeGRangesFromDataFrame(DNaseH, seqnames.field="V1", start.field="V2", end.field="V3", keep.extra.columns = TRUE)

#combine histone mods
#H3K4Me1$ID <- "H3K4Me1"
#H3K4Me3$ID <- "H3K4Me3"
#H3K27Ac$ID <- "H3K27Ac"
#hm <- rbind(H3K4Me1, H3K27Ac,H3K4Me3)
#names(hm) <- c("chromosome","start","end", "signal", "ID")
#hm <- makeGRangesFromDataFrame(hm, seqnames.field="chromosome", start.field="start", end.field="end", keep.extra.columns = TRUE)
#head(hm)
#colnames(mcols(hm))

#TF data
TF <- read.table(TFpath, header = FALSE)
colnames(TF) <- c("chromosome", "start","end", "id")
print(head(TF))

#data from bed file for the coordinates: ChatGPT create the following beauty to prevent my snakemake script from failing in case the BED is
tryCatch({
  bed_data <- read.table(bedpath, header = FALSE, sep = "\t")
}, error = function(e) {
  if (inherits(e, "error") && grepl("no lines available", e$message)) {
    # file is empty
    cat("The BED file is empty. No data was read.\n")
    q(save = "no", status = 0)
  } else {
    # different error
    stop(e)
  }
})
#bed_data <- read.table(bedpath, header = FALSE, sep = "\t")
bed_data <- bed_data[grep("^chr", bed_data$V1), ] #drop rows that do not start with chr
bed_data <- bed_data[!grepl("_", bed_data$V1), ] #drop rows that contain _

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
  #biomTrack <- BiomartGeneRegionTrack(genome = ref_genome, chromosome = chromosome, 
 #                                     start = start_coord, end = end_coord, #filter = list(with_hgnc = TRUE),
  #                                    name = "ENSEMBL", biomart = bm,stacking = "squish",col="black",collapseTranscripts=TRUE,
  #                                    transcriptAnnotation="gene", fill="#fba247")
  #biomTrack <- BiomartGeneRegionTrack(genome = ref_genome, chromosome = chromosome, start = start_coord, end = end_coord, biomart = bm, stacking ="squish",stackHeight=0.5,transcriptAnnotation="symbol",fill = "salmon")
  
  gtfTrack <- GeneRegionTrack(txdb, chromosome = chromosome, start=start_coord, end=end_coord, collapseTranscripts=TRUE, geneSymbol = TRUE, transcriptAnnotation = "symbol", fill="#fba247",stacking = "squish",just.group="above")
   z <- ranges(gtfTrack)
   z$gene <- sub("\\.\\d+", "", z$gene)
   print(z)
   
   tryCatch({
   	z$symbol <- mapIds(Homo.sapiens, z$gene, column = "SYMBOL", keytype = "ENSEMBL")
   	keep <- !is.na(mcols(z)$symbol)
   	z <- z[keep,]
   	}, error = function(e) {
	   # Handle the error (optional)
	   cat("An error occurred:", conditionMessage(e), "\n")
	   # Continue with the script or provide a default behavior
	   }, finally = {
	   # Code to execute regardless of success or failure
	   cat("Continuing with the script...\n")
	   })
   ranges(gtfTrack) <- z

  #bed
  #1: Load BED of Insertions
  coord <- AnnotationTrack(bed, name = "Insertions")

  #functional genomics
  #H3K4Me1_track <- DataTrack(H3K4Me1, name="H3K4Me1", type=("heatmap"), ylim = c(0,30))
  #H3K4Me3_track <- DataTrack(H3K4Me3, name="H3K4Me3", type=("heatmap"), ylim = c(0,30))
  #H3K27Ac_track <- DataTrack(H3K27Ac, name="H3K27Ac", type=("heatmap"), ylim = c(0,30))
  
  
  #TF
  #TF_track <- AnnotationTrack(TF, name="TF Clusters", stacking="squish", featureAnnotation = "id",fontcolor.item="red",
  #                          cex=0.5, fill="white", lty=0)
  #print(TF_track)
  #2: Load BAM Coverage
  altrack <- AlignmentsTrack(bampath, isPaired = FALSE, type = c("coverage")) #bampath, isPaired = FALSE
	
  #3: Use respective Ideogram and GenomeAxis (automatic)
  itrack <- IdeogramTrack(genome = ref_genome, chromosome = chromosome)
  gtrack <- GenomeAxisTrack()
  
  #4: Plot all
  filename <- sprintf("%s_start%s_end%s.pdf",chromosome, start, stop) #itrack, altrack,, DNaseH
  pdf(file=paste(outputpath,filename, sep="/"), width=10, height=10)
  plotTracks(c(itrack,gtrack,gtfTrack,altrack,coord), #H3K4Me1_track,H3K4Me3_track,H3K27Ac_track,TF_track #,0.25,0.25,0.25,5
             chromosome =chromosome, from = start_coord,
             to = end_coord, transcriptAnnotation="symbol", col="black",background.title = "red4",col.main="red", outerMargin=0, innerMargin=5, sizes = c(0.25,0.25,0.25,0.5,0.25)) 
  #plotTracks(list(itrack,gtrack,biomTrack,altrack,coord,H3K4Me1_track,H3K4Me3_track,H3K27Ac_track,DNaseH), chromosome =chromosome, from = start_coord, to = end_coord, background.title = "lightblue", sizes = c(0.25,0.25,1,1,0.25,0.25,0.25,0.25,0.25), #red2,red4,royalblue2
  #col.main="red", innerMargin=10) 
  dev.off()
}
