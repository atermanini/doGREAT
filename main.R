# Copyright 2016 Alberto Termanini

# USAGE: R --vanilla --slave --args --help < main.R;



#----------------------- GLOBALS:
rm(list=ls());
Sys.setlocale("LC_TIME", "en_US.UTF-8");
options(stringsAsFactors = F); 

#----------------------- LIBRARIES:
library(getopt);
library(rGREAT);

#----------------------- CONSTANTS:
ontology = c("GO Molecular Function", 
             "GO Biological Process", 
             "GO Cellular Component", 
             "Mouse Phenotype", 
             "Human Phenotype", 
             "Disease Ontology", 
             "MSigDB Oncogenic Signatures",
             "MSigDB Immunologic Signatures",
             "PANTHER Pathway",
             "BioCyc Pathway", 
             "MSigDB Pathway",
             "MGI Expression: Detected",
             "MSigDB Perturbation",
             "MSigDB Predicted Promoter Motifs",
             "MSigDB miRNA Motifs",
             "InterPro",
             "TreeFam");

#----------------------- FUNCTIONS:
verbose = function(txt) {
  dt = format(Sys.time(), "%Y-%m-%d %H:%M:%S");
  write(paste0("[", dt, "] ", txt), stderr());
}

#----------------------- PARAMETERS:
#0: no argument
#1: required argument
#2: optional argument
#args types: logical, integer, double, character
m = matrix(c(
  "help"   , "h", "0", "logical", "this help",
  "verbose", "v", "0", "logical", "verbose mode on",
  "workdir",  "w", "2", "character", "working directory (default = current directory)",
  "save-image", "s", "2", "character", "save workspace to: workdir/data/image.Rda",
  "infile",  "i", "1", "character", "input file (BED file with genomic regions)",
  "release",  "r", "1", "character", "genome release (eg. mm9)",
  "outdir", "o", "1", "character", "output directory (default = workdir/results)",
  "FDR", "F", "2", "integer", "FDR threshold (default = 0.05)",
  "fold", "f", "2", "integer", "binomial fold enrichment threshold (default = 2)"
  
), byrow=TRUE, ncol=5);
opt = getopt(spec = m, opt = commandArgs(TRUE));

# help:
if ( !is.null(opt$help) ) {
  cat(getopt(m, usage=TRUE));
  quit(status=1);
}

# defaults:
if ( is.null(opt$"verbose") )      { opt$"verbose" = FALSE; }
if ( is.null(opt$"workdir") )      { opt$"workdir" = getwd(); }
if ( is.null(opt$"outdir") )       { opt$"outdir" = paste0(opt$"workdir","/results"); }
if ( is.null(opt$"save-image") )   { opt$"save-image" = FALSE; } else if ( opt$"save-image"==TRUE ) { opt$"save-image" = paste0(opt$"workdir", "/cummerbund/data/image.Rda"); }
if ( is.null(opt$"FDR") )          { opt$"FDR" = 0.05; } else { opt$"FDR" = as.numeric(opt$"FDR"); }
if ( is.null(opt$"fold") )         { opt$"fold" = 2; } else { opt$"fold" = as.numeric(opt$"fold"); }

# print values:
if (opt$"verbose" == TRUE)  { print(opt); }

# requirements:
if ( is.null(opt$"infile") )		{ quit(status = 1); }
if ( is.null(opt$"release") )		{ quit(status = 1); }


#----------------------- WORKING DIRECTORY:
if (opt$"verbose" == TRUE) { verbose(paste0("Setting working dir to: ", opt$workdir)); }
setwd(opt$workdir);


#----------------------- READING DATA:
if (opt$"verbose" == TRUE) { verbose("Reading data"); };

if(file.exists(opt$"infile")) {
  
  data = read.table(opt$"infile", sep = "\t", header = FALSE);
  colnames(data) = c("chr","start","end");

} else {
  
  verbose(paste0("ERROR: file does not exist, exit. File: ", opt$infile));
  quit(status = 1);
}


#----------------------- DATA PROCESSING:
if (opt$"verbose" == TRUE) { verbose("Processing data"); }

job = submitGreatJob(gr = data, species = opt$"release"); # submit job to GREAT server
if (opt$"verbose" == TRUE) print(job); # print job details
l = getEnrichmentTables(job, ontology = ontology); # download job results
for (n in names(l)) {
  df = l[[n]]; 
  df$FDR_Binom = p.adjust(df$"Binom_Raw_PValue", method="fdr");
  df$FDR_Hyper = p.adjust(df$"Hyper_Raw_PValue", method="fdr");
  l[[n]] = subset.data.frame(df, FDR_Binom <= opt$"FDR" & FDR_Hyper <= opt$"FDR" & Binom_Fold_Enrichment >= opt$"fold"); # filtering results (as in GREAT web page)
}


#----------------------- WRITING DATA:
if (opt$"verbose" == TRUE) { verbose("Writing data"); }
dir.create(opt$"outdir", showWarnings = T);
for (n in names(l)) {
  write.table(l[[n]], file = paste0(opt$"outdir", "/", n, ".xls"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t");
}


#----------------------- SAVE IMAGE:
if(is.character(opt$"save-image")) {
  if(opt$"verbose") { verbose(paste0("Saving image to: ", opt$"save-image")); }
  save.image(file = opt$"save-image", compress = TRUE, safe = TRUE);
}


#----------------------- SESSION INFO:
if(opt$"verbose") {

  sink(file = stderr());
  print(sessionInfo());    
}

#----------------------- EXIT:
quit(status=0);
