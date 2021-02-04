
rm(list = ls())

library(reshape2)
library(sleuth)
library(dplyr)
# library(car)
library(ggplot2)
# library(broom)

## -----------------------------------------------------

## settings
if(TRUE) {
  
  use.ir.geneninfo.grch37.flag = 0
  use.ir.geneninfo.grch38.flag = 1
  if(use.ir.geneninfo.grch37.flag == 1) {use.ir.geneninfo.grch38.flag = 0}
  if(use.ir.geneninfo.grch38.flag == 1) {use.ir.geneninfo.grch37.flag = 0}
  
  ir.string = ''
  if(use.ir.geneninfo.grch37.flag == 1)  { ir.string = '.geneinfo_grch37'  }
  if(use.ir.geneninfo.grch38.flag == 1)  { ir.string = '.geneinfo_grch38'  }
  
  # ir.ratio.min = 0
  
  cov.min = 0 ## recommended by the IRFinder paper; the paer also uses IRratio>0.1; reset after I used lowDP6
  flt.string = ''
  if(cov.min >0) {
    flt.string = paste0('.IntronDepth',cov.min)
  }
  
}

## -----------------------------------------------------

path = '/gpfs/data/pitroda-lab/RZ73/irfinder'
setwd(path)

## -----------------------------------------------------

## data files 
if(TRUE) {
  transxinfo.file = 'gencode.v28.primary_assembly.annotation.maskPAR.gtf.transcriptinfo'
  # sleuth.file = 'CRC.kallisto_sleuth.files' ## grch38
  # metadata.file = '../sampleinfo/Gide.rnaseq.metadata.txt'
  # clinical.file = 'C:/Users/rbaoc/Desktop/current_proj/00_databases/pancancer/clinical/TCGA-CDR.txt'
  ir.path = path
  ir.file = 'lungs.irfinder.clean_LowCovDP6.tsv' ## grch37
  
  if(use.ir.geneninfo.grch37.flag == 1) {
    ir.geneinfo.file = 'gencode.v19.chr_patch_hapl_scaff.annotation.gtf.geneinfo'
  } else if(use.ir.geneninfo.grch38.flag == 1) {
    ir.geneinfo.file = 'gencode.v28.primary_assembly.annotation.maskPAR.gtf.geneinfo'
  }
  
  # kallisto.file = 'TCGASKCM.rnaseq.kallisto.raw.txi.coding.CPM.sm93.csv' ## grch38
  # gene.length.file = 'TCGASKCM.rnaseq.kallisto.length.txt.max_by_gene.csv' ## grch38

}

## -----------------------------------------------------

## import data
if(TRUE) {
  
  ir = read.delim(file.path(ir.path, ir.file),stringsAsFactors = F,row.names=NULL)
  newnames = colnames(ir)
  newnames[1] = "Sample"
  colnames(ir)<-newnames
  transxinfo = read.delim(file.path("/gpfs/data/pitroda-lab/ReferenceData/kallisto_indexes/GRCh38.primary_Gencode28_slim_maskPAR", transxinfo.file), stringsAsFactors = F)
  # subtype = read.delim(subtype.file, stringsAsFactors = F)
  ir.geneinfo = read.delim(file.path("/gpfs/data/pitroda-lab/ReferenceData/kallisto_indexes/GRCh38.primary_Gencode28_slim_maskPAR", ir.geneinfo.file), stringsAsFactors = F)
  # kallisto = read.delim(kallisto.file, stringsAsFactors = F)
  # gene.length = read.csv(gene.length.file, stringsAsFactors = F)
  
  # load(paste0(ir.file,'.ir.RData'))
  # load(paste0(ir.file, '.ir.wPt.RData'))
  
  ## -----------------------------------------------------
  
  # load(paste0(sleuth.file,'.so.RData'))
  # so.norm.matrix = sleuth_to_matrix(so, "obs_norm", "tpm")$data
  # save(so.norm.matrix, file = paste0(sleuth.file, '.so.norm.matrix.RData'))
  
  # load(paste0(sleuth.file, '.so.norm.matrix.RData'))
  # sleuth = so.norm.matrix
  # rm(so, so.norm.matrix)
  
}

## -----------------------------------------------------

## preprocessing 
if(TRUE) {
  
  # head(subtype) ## 28 patients
  
  sample.ex = c('SRR3184292', ## on-treatment
                'SRR3184299' ## dup rnaseq sample for sample patient pt27
  )
  
  # subtype = subtype[!subtype$Sample %in% sample.ex,]
  # dim(subtype) ## 26
  
  ## -----------------------------------------------------

  dim(ir) ## 33666666       25
  # dim(ir.wPt) ## 2356703      26
  # ir.wPt.stats = data.frame(table(ir.wPt$PatientID)) ## 95 patients ...
  
  ## -----------------------------------------------------
  
  ir[1:3,]
  table(ir$Warnings)
  # -              LowCover NonUniformIntronCover 
  # 19000115               5558545               9108006 
  
  table(ir$StaticWarning)
  # clean 
  # 33666666 
  
  ## -----------------------------------------------------
  
  ## check distribution of coverage....
  if(FALSE) {
    
    summary(ir$Coverage)
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 0.0000  0.3731  0.6509  0.5926  0.8476  1.0000 
    summary(ir$IntronDepth)
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    # 0.0000   0.1959   1.0000   1.5546   1.8946 414.1250 
    
    p0 = ggplot(ir, aes(log10(IntronDepth))) +
      geom_density(aes(group = Sample), linetype = 'dotted') +
      geom_vline(xintercept = log10(3), color = '#CC0000')
    
    pdf(file = paste0(ir.file, '.IntronDepth.density.pdf'),width = 8, heigh =6)
    print(p0)
    dev.off()
    
  }
  
  ## -----------------------------------------------------
  
  # ## 10/07/2018: filter by intron depth (cov >=3 in the paper)
  # if(cov.min > 0) {
  #   dim(ir) ## 3410421
  #   
  #   ir = ir[!is.na(ir$IntronDepth) & ir$IntronDepth >= cov.min,]
  #   
  #   dim(ir) ## 406335
  #   ## ~82% of the IRs were removed????
  #   
  #   write.table(ir,
  #               file = paste0(ir.file,flt.string,'.tsv'),
  #               sep = '\t', row.names = F, col.names = T, quote = F)
  #   
  # }
  # 
  
  ## -----------------------------------------------------
  
  head(ir.geneinfo)
  sum(duplicated(ir.geneinfo$gene_id))
  
  ir.geneinfo$gene_id = gsub('[.]\\d+$','',ir.geneinfo$gene_id)
  sum(duplicated(ir.geneinfo$gene_id)) ## 0
  
  table(ir.geneinfo$gene_type)
  table(ir.geneinfo$transcript_type) ## no retained_intron type in grch37 gtf....
  
  ir.geneinfo.coding = ir.geneinfo[ir.geneinfo$gene_type=='protein_coding',]
  # table(ir.geneinfo.coding$transcript_type)
  dim(ir.geneinfo.coding)
  # protein_coding 
  # 22705
  
  ## -----------------------------------------------------
  
  ## select only IRs identified inprotein coding genes 
  ir.coding = ir[ir$EnsemblID %in% ir.geneinfo.coding$gene_id,]
  dim(ir.coding) ## 3401291      25; after filtering IntronDepth>=3: 403272 25
  
  # ## select sleuth samples in subtype list...
  # ir.coding = ir.coding[ir.coding$Sample %in% subtype$Sample,]
  # dim(ir.coding) ## 2307703      25; after filtering IntronDepth>=3: 403272     25
  
}

## -----------------------------------------------------

## how many IRs per patient??
if(TRUE) {
  
  ir.coding.stats = cbind(
    data.frame(table(ir.coding$Sample[ir.coding$IRratio  >= 0.00 & ir.coding$IRratio <= 0.80])),
    data.frame(table(ir.coding$Sample[ir.coding$IRratio  >= 0.05 & ir.coding$IRratio <= 0.80])),
    data.frame(table(ir.coding$Sample[ir.coding$IRratio  >= 0.10 & ir.coding$IRratio <= 0.80])),
    data.frame(table(ir.coding$Sample[ir.coding$IRratio  >= 0.15 & ir.coding$IRratio <= 0.80])),
    data.frame(table(ir.coding$Sample[ir.coding$IRratio  >= 0.20 & ir.coding$IRratio <= 0.80])))
  
  dim(ir.coding.stats) ## 26
  
  print(all.equal(ir.coding.stats[,1], ir.coding.stats[,3]))
  print(all.equal(ir.coding.stats[,1], ir.coding.stats[,5]))
  print(all.equal(ir.coding.stats[,1], ir.coding.stats[,7]))
  print(all.equal(ir.coding.stats[,1], ir.coding.stats[,9]))
  
  ir.coding.stats = ir.coding.stats[,c(1,2,4,6,8,10)]
  
  colnames( ir.coding.stats) = c('Sample','IRratio_0.00','IRratio_0.05',
                                 'IRratio_0.10','IRratio_0.15','IRratio_0.20')
  
  # ir.coding.stats = merge(subtype, ir.coding.stats, by = 'Sample')
  dim(ir.coding.stats)
  head(ir.coding.stats)
  
  write.csv(ir.coding.stats,
              file = paste0(ir.file,ir.string,flt.string,'.ir.coding.stats.csv'),
              row.names = F)
}

## -----------------------------------------------------

## does lengthScaledTPM correlate with gene length???
## it seems so................... :( 
if(FALSE) {
  
  dim(kallisto)
  dim(gene.length)
  
  gene.length$gene_id = gsub('[.]\\d+$','',gene.length$gene_id)
  sum(duplicated(gene.length$gene_id))
  
  row.names(kallisto) = kallisto$X
  
  kallisto.stats = merge(data.frame(mean = 
                                      apply(kallisto[,-1]+1, 1, mean, na.rm=T)),
                         data.frame(median =
                                      apply(kallisto[,-1]+1, 1, median, na.rm=T)),
                         by = 'row.names')
  dim(kallisto.stats) ##  19883     3 
  
  row.names(kallisto.stats) = kallisto.stats $Row.names
  colnames(kallisto.stats)[1] = 'Gene'
  
  kallisto.stats$gene_id = gsub('\\S+[!]','',kallisto.stats$Gene)
  kallisto.stats = merge(kallisto.stats, gene.length, by = 'gene_id')
  dim(kallisto.stats) ## 19883 genes 
  
  cor.test(kallisto.stats$mean, kallisto.stats$length, method = 'spearman')
  ## 0.535731
  
  cor.test(kallisto.stats$median, kallisto.stats$length, method = 'spearman')
  ## 0.5236445
  
  p1 = ggplot(kallisto.stats, aes(length, mean)) +
    geom_point() +
    # geom_smooth(method = lm,formula = y ~ poly(x, 3), se = FALSE) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks()  +
    theme_classic()
  
  print(p1)
}

## -----------------------------------------------------

## how many genes are expressed per patient??
if(FALSE) {
  
  dim(kallisto) ## 19883    94
  row.names(kallisto) = kallisto$Gene
  
  kallisto.gene.expressed = 
    cbind(CPM_0 = apply(kallisto[,-1], 2, function(x) sum(x>=0)),
          CPM_1 = apply(kallisto[,-1], 2, function(x) sum(x>=1)),
          CPM_2 = apply(kallisto[,-1], 2, function(x) sum(x>=2)),
          CPM_3 = apply(kallisto[,-1], 2, function(x) sum(x>=3)),
          CPM_4 = apply(kallisto[,-1], 2, function(x) sum(x>=4)),
          CPM_5 = apply(kallisto[,-1], 2, function(x) sum(x>=5)))
  
  kallisto.gene.expressed = data.frame(Sample = row.names(kallisto.gene.expressed),
                                       kallisto.gene.expressed,
                                       stringsAsFactors = F)
  
  
  
}

## -----------------------------------------------------

# ir.coding.stats = merge(ir.coding.stats, kallisto.gene.expressed, by = 'Sample')
# dim(ir.coding.stats)
# 
# write.csv(ir.coding.stats,
#             file = paste0(ir.file,ir.string,flt.string,'.ir.coding.stats.csv'),
#             row.names = F)


## -----------------------------------------------------

sessionInfo()



