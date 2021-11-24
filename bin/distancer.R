#!/usr/bin/env Rscript
# SShekarriz Jun03,2021
# find distance between folp and glm
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "input.txt"
  args[2] = "output.txt"
}

#################
library(tidyverse)
#################
cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
          "gapopen", "qstart", "qend", "sstart", "send", 
          "evalue", "bitscore")
input= read.csv(args[1], sep = "\t",
                header = F)
colnames(input) <- cols
input %>%
  filter(pident >= 90) %>%
  select(qseqid, sseqid, sstart, send) %>%
  gather(order, position, -qseqid, -sseqid) %>%
  mutate(gene= paste(qseqid, order, sep = "_")) %>%
  select(sseqid, gene, position) %>%
  spread(gene, position) %>%
  mutate(distance= folP_send - glmM_sstart) %>%
  write.table(args[2], sep="\t", quote=F, row.names=F)	
