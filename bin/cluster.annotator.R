#!/usr/bin/env Rscript
# SShekarriz Jun03,2021
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "out.gff"
  args[2] = ".cluster-names"
  args[3] = "_contig.names"
  args[4] = "output"
}

#################
library(tidyverse)
#################

gff_file=args[1]
cluster.name=args[2]
contig.name=args[3]

#cluster names
cluster.name= read.csv(cluster.name, header = F)
cluster.name %>%
  mutate(V1=gsub(">", "", V1),
         V2=gsub(">", "", V2)) -> cluster.name
colnames(cluster.name) <- c("Cluster", "seqname")
#contig and genome name
contig.name=read.csv(contig.name, header = F)
contig.name %>%
  mutate(V1= gsub(".*genomes/", "", V1)) %>%
  separate(V1, c("genome_name", "contig_name"), sep = ":>") %>%
  mutate(genome_name=gsub(".fasta", "", genome_name))-> contig.name
#reading gff file
gff_cols <- c("seqname","source","feature","start","end","score","strand",
              "frame","attributes")
read.delim(gff_file, header=F, comment.char="#") -> gff
gff %>%
  filter(grepl("^cluster_.*", V1)) -> gff
colnames(gff) <- gff_cols

gff %>%
  left_join(cluster.name, by = "seqname") %>%
  separate(Cluster, c("contig_name", "start-stop"), sep = ":") %>%
  right_join(contig.name) %>%
  write.table(args[4], sep="\t", quote=F, row.names=F)


