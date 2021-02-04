
rm(list = ls())

library(reshape2)
library(sleuth)
library(dplyr)
library(car)
library(ggplot2)
library(broom)

## -----------------------------------------------------

sample.ex = c()
sample.ex.string = ''
if(length(sample.ex)>0) { sample.ex.string = paste0('.ex',paste(sample.ex,collapse = '_') )}
print(sample.ex.string)

PR.ex.flag = 0
out.string = ''
if(PR.ex.flag==1) { out.string = '.exPR' }

## -----------------------------------------------------

path = '/gpfs/data/pitroda-lab/ncbi/GSE91061/RI'
setwd(path)

## -----------------------------------------------------

transx.file = '/gpfs/data/pitroda-lab/ReferenceData/kallisto_indexes/GRCh38.primary_Gencode28_slim_maskPAR/gencode.v28.primary_assembly.annotation.maskPAR.gtf.transcriptinfo'
rnaseq.file = 'GSE91061.kallisto_sleuth.files'
# subtype.file = 'C:/Users/rbaoc/Desktop/current_proj/spitrodaProj/CRC_neoantigen/data/sean/Subtypes for CRC.txt'
sampleinfo.file = 'Riaz.rnaseq.metadata.txt.wClinical.txt'

## -----------------------------------------------------

transx = read.delim(transx.file, stringsAsFactors = F)
sampleinfo = read.delim(sampleinfo.file, stringsAsFactors = F)

## -----------------------------------------------------

## only need to run once!
if(FALSE) {
  load(paste0(rnaseq.file,'.so.RData'))
  so.norm.matrix = sleuth_to_matrix(so, "obs_norm", "tpm")
  # so.norm.matrix = so.norm.matrix$data
  so.norm.matrix = data.frame(so.norm.matrix, stringsAsFactors = F)
  save(so.norm.matrix, file = paste0(rnaseq.file, '.so.norm.matrix.RData'))

}

load(paste0(rnaseq.file, '.so.norm.matrix.RData'))
rnaseq = so.norm.matrix
rm(so, so.norm.matrix)

## -----------------------------------------------------

## join sampleinfo and group
if(TRUE) {

  sampleinfo = sampleinfo[! sampleinfo$Sample %in% sample.ex,]
  dim(sampleinfo) ## 109

  table(sampleinfo$Response)
  # CR NE PD PR SD
  # 6  3 44 14 34

  table(sampleinfo$Timepoint)
  # On Pre
  # 58  51

  table(sampleinfo$Subtype)
  # ACRAL    CUTANEOUS      MUCOSAL OCULAR/UVEAL        OTHER
  # 2           66           12            8           13

  sampleinfo$BORR = NA
  sampleinfo$BORR[which(sampleinfo$Response %in%
                          c('CR','PR'))] = 'R'
  sampleinfo$BORR[which(sampleinfo$Response %in%
                          c('PD','SD'))] = 'NR'
  table(sampleinfo$BORR)
  # NR  R
  # 78 20

  write.csv(sampleinfo,
            file = paste0(sampleinfo.file,'.wGroup',
                          sample.ex.string,
                          '.csv'),
            row.names = F)
}


if(PR.ex.flag == 1) {

  sampleinfo = sampleinfo[sampleinfo$anti_pd_1_response !='Partial Response',]
}
dim(sampleinfo) ## 109

## -----------------------------------------------------

dim(transx) ## 203742
dim(rnaseq) ## 203675

## transcripts dropped by sleuth??? why???
dim(transx[!transx$transcript_id %in% row.names(rnaseq),]) ## 67

## -----------------------------------------------------

## all transcripts from protein_coding genes
transx.coding = transx[transx$gene_type %in% c('protein_coding'),]
## ri transcripts (regardless of gene_type)
transx.ri = transx[transx$gene_name %in% c('retained_intron'),]
## note gene_type=protein_coding may have transcript_type=retained_intron
## coz some transcripts from protein_coding genes may have retained_intron

dim(transx.coding) ## 149520
dim(transx.ri) ## 27819

## -----------------------------------------------------

## all transcripts from protein_coding genes
rnaseq.coding = rnaseq[row.names(rnaseq) %in% transx.coding$transcript_id,]
## ri transcripts (regardless of gene_type)
rnaseq.ri= rnaseq[row.names(rnaseq) %in% transx.ri$transcript_id,]
## ri transcripts from protein_coding genes
rnaseq.coding.ri= rnaseq[row.names(rnaseq) %in%
                           intersect(transx.coding$transcript_id,transx.ri$transcript_id),]

dim(rnaseq.coding) ## 149483
dim(rnaseq.ri) ## 27819
dim(rnaseq.coding.ri) ## 27042

## -----------------------------------------------------

## calculate frac
if(TRUE) {

  rnaseq.perSample = cbind(data.frame(total = apply(rnaseq, 2, sum, na.rm=T)),
                           data.frame(coding = apply(rnaseq.coding, 2, sum, na.rm=T)),
                           data.frame(ri = apply(rnaseq.ri, 2, sum, na.rm=T)),
                           data.frame(coding.ri = apply(rnaseq.coding.ri, 2, sum, na.rm=T)),
                           data.frame(total.tranx =
                                        apply(rnaseq, 2, function(x) sum(x>0, na.rm=T))),
                           data.frame(coding.tranx =
                                        apply(rnaseq.coding, 2, function(x) sum(x>0, na.rm=T))),
                           data.frame(ri.tranx =
                                        apply(rnaseq.ri, 2, function(x) sum(x>0, na.rm=T))),
                           data.frame(coding.ri.tranx =
                                        apply(rnaseq.coding.ri, 2,
                                              function(x) sum(x>0, na.rm=T))))

  dim(rnaseq.perSample) ## 109

  rnaseq.perSample$coding.frac = rnaseq.perSample$coding / rnaseq.perSample$total
  rnaseq.perSample$ri.frac = rnaseq.perSample$ri / rnaseq.perSample$total
  rnaseq.perSample$coding.ri.frac = rnaseq.perSample$coding.ri / rnaseq.perSample$coding
  rnaseq.perSample$coding.tranx.frac =
    rnaseq.perSample$coding.tranx / rnaseq.perSample$total.tranx
  rnaseq.perSample$ri.tranx.frac = rnaseq.perSample$ri.tranx / rnaseq.perSample$total.tranx
  rnaseq.perSample$coding.ri.tranx.frac =
    rnaseq.perSample$coding.ri.tranx / rnaseq.perSample$coding.tranx

  write.csv(rnaseq.perSample,
            file = paste0(rnaseq.file,sample.ex.string,'.perSample.stats.csv'))
}

## -----------------------------------------------------

row.names(sampleinfo) = sampleinfo$Sample
sum(duplicated(sampleinfo$Patient))
row.names(sampleinfo)[!row.names(sampleinfo) %in% colnames(rnaseq)] ## 0

## -----------------------------------------------------

## calculate relative enrichment of ri TPMs over all transcript TPMs
if(TRUE) {

  dim(transx.coding) ## 149520
  dim(transx.ri) ## 27819
  dim(rnaseq.coding) ## 149483
  dim(rnaseq.coding.ri) ## 27042

  ## protein_coding genes that have retained intron transcripts
  genes.coding.ri = unique(transx$gene_id[transx$gene_type %in% c('protein_coding') &
                              transx$gene_name %in% c('retained_intron')])
  length(genes.coding.ri) ## 9164
  length(unique(transx$gene_id[transx$gene_type %in% c('protein_coding')])) ## 19912

  genes.coding.ri.frac = NULL
  genes.coding.ri.tpm.1 = NULL ## ri transcripts sum TPM
  genes.coding.ri.tpm.2 = NULL ## ri-containing genes sum TPM

  ## calculate frac
  for(i in 1:length(genes.coding.ri)) {
    print(i)
    my.gene = genes.coding.ri[i]
    my.tranx = transx$transcript_id[transx$gene_id==my.gene]
    my.tranx.ri = transx$transcript_id[transx$gene_id==my.gene &
                                         transx$gene_name %in% c('retained_intron')]

    if(length(my.tranx)==0 | length(my.tranx.ri)==0) { break }

    my.data = NULL
    my.data = t(rnaseq.coding[row.names(rnaseq.coding)%in%my.tranx,,drop=F])
    my.data = data.frame(coding = apply(my.data, 1, sum, na.rm=T))

    my.data.ri = NULL
    my.data.ri = t(rnaseq.coding[row.names(rnaseq.coding)%in%my.tranx.ri,,drop=F])
    my.data.ri = data.frame(coding.ri = apply(my.data.ri, 1, sum, na.rm=T))

    ## calculate ratio
    my.data.ri.frac = NULL
    my.data.ri.frac =  my.data.ri / my.data

    if(i==1) {
      genes.coding.ri.frac =  my.data.ri.frac
      genes.coding.ri.tpm.1 = my.data.ri
      genes.coding.ri.tpm.2 = my.data
    } else {
      genes.coding.ri.frac = cbind(genes.coding.ri.frac,my.data.ri.frac)
      genes.coding.ri.tpm.1 = cbind(genes.coding.ri.tpm.1, my.data.ri)
      genes.coding.ri.tpm.2 = cbind(genes.coding.ri.tpm.2, my.data)
    }

    colnames(genes.coding.ri.frac)[ncol(genes.coding.ri.frac)] = my.gene
    colnames(genes.coding.ri.tpm.1)[ncol(genes.coding.ri.tpm.1)] = my.gene
    colnames(genes.coding.ri.tpm.2)[ncol(genes.coding.ri.tpm.2)] = my.gene

  }

  dim(genes.coding.ri.frac) ## 91 samples, 9164 coding genes
  dim(genes.coding.ri.tpm.1) ## 91 samples, 9164 coding genes
  dim(genes.coding.ri.tpm.2) ## 91 samples, 9164 coding genes

  ## NAs in genes.coding.ri.frac are the genes with
  ## 0 TPM in genes.coding.ri.tpm.2, hence the fraction cannot be
  ## calculated (I have manuallly confirmed this!!)
  genes.coding.ri.frac = data.frame(t(genes.coding.ri.frac))
  genes.coding.ri.tpm.1 = data.frame(t(genes.coding.ri.tpm.1))
  genes.coding.ri.tpm.2 = data.frame(t(genes.coding.ri.tpm.2))

  write.csv(genes.coding.ri.frac,
            file = paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.frac.sm',
                          ncol(genes.coding.ri.frac),'.csv'))
  write.csv(genes.coding.ri.tpm.1,
            file = paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.tpm.1.sm',
                          ncol(genes.coding.ri.tpm.1),'.csv'))
  write.csv(genes.coding.ri.tpm.2,
            file = paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.tpm.2.sm',
                          ncol(genes.coding.ri.tpm.2),'.csv'))


  ## select only the samples present in group file
  genes.coding.ri.frac = genes.coding.ri.frac[,row.names(sampleinfo)]
  genes.coding.ri.tpm.1 = genes.coding.ri.tpm.1[,row.names(sampleinfo)]
  genes.coding.ri.tpm.2 = genes.coding.ri.tpm.2[,row.names(sampleinfo)]
  dim(genes.coding.ri.frac) ## 91

  write.csv(genes.coding.ri.frac,
            file = paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.frac.sm',
                          ncol(genes.coding.ri.frac),'.csv'))
  write.csv(genes.coding.ri.tpm.1,
            file = paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.tpm.1.sm',
                          ncol(genes.coding.ri.tpm.1),'.csv'))
  write.csv(genes.coding.ri.tpm.2,
            file = paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.tpm.2.sm',
                          ncol(genes.coding.ri.tpm.2),'.csv'))

}

## -----------------------------------------------------

genes.coding.ri.frac = read.csv(paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.frac.sm109.csv'),
                                stringsAsFactors = F,
                                row.names = 1)
genes.coding.ri.tpm.1 = read.csv(paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.tpm.1.sm109.csv'),
                                stringsAsFactors = F,
                                row.names = 1)
genes.coding.ri.tpm.2 = read.csv(paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.tpm.2.sm109.csv'),
                                stringsAsFactors = F,
                                row.names = 1)
dim(genes.coding.ri.frac) ## 9164   109

## -----------------------------------------------------

## run test on genes.coding.ri.frac matrix
sample.min = round(0.5 * min(table(sampleinfo$BORR)))
print(sample.min) ## 10

row.sum.min = 0.05

if(TRUE) {

  ## distribution
  if(TRUE) {
    data.plot = melt(genes.coding.ri.frac)
    data.plot$logit = logit(data.plot$value, adjust = 0)
    ## adjust will put artificial peaks at both ends of data distribution.
    ## turn it to 0

    p1 = ggplot(data.plot, aes(value)) +
      geom_density() +
      geom_vline(xintercept = row.sum.min)

    p2 = ggplot(data.plot, aes(logit)) +
      geom_density()
    ## Removed 142302 rows containing non-finite values (stat_density).
    ## cannot deal with 0 or 1 values - logit function
  }

  ## filtering
  if(TRUE) {

    ## remove lowly frac genes
    row.select = which(apply(genes.coding.ri.frac, 1,
                             function(x) sum(x>row.sum.min, na.rm=T))
                       >= sample.min)
    length(row.select) ## 6230

    if(length(genes.coding.ri.frac)>0) {
      genes.coding.ri.frac.flt = genes.coding.ri.frac[row.select,,drop=F]
    }

    write.csv(genes.coding.ri.frac.flt,
              file = paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.frac.sm',
                            ncol(genes.coding.ri.frac),'.flt.csv'))

    ## -----------------------------------------------------------

    ## deal with 0 or 1 values....
    ## this actually change data ... skip it
    # my.min = min(genes.coding.ri.frac.flt[genes.coding.ri.frac.flt>0],na.rm=T)
    # genes.coding.ri.frac.flt[genes.coding.ri.frac.flt==0] = my.min
    # genes.coding.ri.frac.flt[genes.coding.ri.frac.flt==1] = 0.999

    ## transform proportion numbers ...
    genes.coding.ri.frac.flt.logit =
      apply(genes.coding.ri.frac.flt, 1, logit, adjust=0)
    genes.coding.ri.frac.flt.logit = data.frame(t(genes.coding.ri.frac.flt.logit))

    write.csv(genes.coding.ri.frac.flt.logit,
              file = paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.frac.sm',
                            ncol(genes.coding.ri.frac),'.flt.logit.csv'))

    ## distribution
    if(TRUE) {
      data.plot = melt(genes.coding.ri.frac.flt)
      data.plot$logit = logit(data.plot$value, adjust = 0)

      p3 = ggplot(data.plot, aes(logit)) +
        geom_density()
      ## Removed 104878 rows containing non-finite values (stat_density).

    }

    ## deal with -Inf and Inf values produced during logit function....
    if(TRUE) {
      my.min = min(genes.coding.ri.frac.flt.logit[
        genes.coding.ri.frac.flt.logit != -Inf &
          genes.coding.ri.frac.flt.logit != Inf &
          !is.na(genes.coding.ri.frac.flt.logit)
        ])
      print(my.min) ## -37.42995

      my.max = max(genes.coding.ri.frac.flt.logit[
        genes.coding.ri.frac.flt.logit != -Inf &
          genes.coding.ri.frac.flt.logit != Inf &
          !is.na(genes.coding.ri.frac.flt.logit)
        ])
      print(my.max) ## 28.88281

      replace.min.max = function(x, min, max) {
        y = x
        row.select = NA
        row.select = which(x==-Inf)
        if(length(row.select)>0) { y[row.select] = min  }
        row.select = NA
        row.select = which(x==Inf)
        if(length(row.select)>0) { y[row.select] = max  }

        return(y)
      }

      # replace.min.max.per.line = function(x) {
      #   y = x
      #   min = min(x[x!=-Inf & x!=Inf & !is.na(x)])
      #   max = max(x[x!=-Inf & x!=Inf & !is.na(x)])
      #   row.select = NA
      #   row.select = which(x==-Inf)
      #   if(length(row.select)>0) { y[row.select] = min  }
      #   row.select = NA
      #   row.select = which(x==Inf)
      #   if(length(row.select)>0) { y[row.select] = max  }
      #
      #   return(y)
      # }

      ## replace -Inf with my.min; Inf with my.max...
      for(i in 1:nrow(genes.coding.ri.frac.flt.logit)) {
        ## this replacement still leads to artifacts of inconsistent
        ## test statistic vs mean frac before logit conversion ...
        genes.coding.ri.frac.flt.logit[i,] =
          replace.min.max(genes.coding.ri.frac.flt.logit[i,],
                          my.min, my.max)
      }

      min(genes.coding.ri.frac.flt.logit, na.rm=T) ## -24.62311
      max(genes.coding.ri.frac.flt.logit, na.rm=T) ## 22.92545

      write.csv(genes.coding.ri.frac.flt.logit,
                file = paste0(rnaseq.file,sample.ex.string,'.genes.coding.ri.frac.sm',
                              ncol(genes.coding.ri.frac),'.flt.logit.replaceInf.csv'))

      ## distribution
      if(TRUE) {
        data.plot = melt(genes.coding.ri.frac.flt.logit)

        p4 = ggplot(data.plot, aes(value)) +
          geom_density()
        ## Removed 1832 rows containing non-finite values (stat_density).

      }

    }


  }

  ## statistical test
  if(TRUE) {

    dim(genes.coding.ri.frac.flt.logit) ## 7072,  109 samples
    group1 = row.names(sampleinfo)[!is.na(sampleinfo$BORR)  &
                                     sampleinfo$BORR=='NR']
    group2 = row.names(sampleinfo)[!is.na(sampleinfo$BORR) &
                                     sampleinfo$BORR=='R']

    genes.coding.ri.frac.flt.logit.stats = NULL
    genes.coding.ri.frac.flt.fc = NULL
    my.row.names = NULL


    for(i in 1:nrow(genes.coding.ri.frac.flt.logit)) {
      print(i)

      my.gene = row.names(genes.coding.ri.frac.flt.logit)[i]
      my.stats = NULL
      my.kept = NULL
      my.ri.tpm = NULL
      my.ri.frac = NULL
      my.ratio = NULL
      my.fc = NULL

      ## require at least 50% samples carry the values
      if(sum(is.na(genes.coding.ri.frac.flt.logit[i,])) >
         round(ncol(genes.coding.ri.frac.flt.logit) * 0.5)) {
        next
      }

      ## 06/20/2018
      ## if data are consistent ... t test would fail!
      if(sd(genes.coding.ri.frac.flt.logit[i,group1], na.rm = T)==0 &
         sd(genes.coding.ri.frac.flt.logit[i,group2], na.rm = T)==0) {
        next
      }
      if(sum(!is.na(genes.coding.ri.frac.flt.logit[i,group1]))<3 |
         sum(!is.na(genes.coding.ri.frac.flt.logit[i,group2]))<3) {
        next
      }

      my.stats = tidy(t.test(genes.coding.ri.frac.flt.logit[i,group1],
                       genes.coding.ri.frac.flt.logit[i,group2]))

      genes.coding.ri.frac.flt.logit.stats = rbind(genes.coding.ri.frac.flt.logit.stats,
                                                   my.stats)
      my.row.names = c(my.row.names, my.gene)

      ## how many samples kept for stats
      my.kept = c(sum(!is.na(genes.coding.ri.frac.flt.logit[i,group1])),
                  sum(!is.na(genes.coding.ri.frac.flt.logit[i,group2])))

      ## calculate mean ri tpm, ratio, and fc
      my.ri.tpm = c(mean(
        unlist(genes.coding.ri.tpm.1[row.names(genes.coding.ri.tpm.1)==my.gene,group1]),
        na.rm=T),
        mean(
          unlist(genes.coding.ri.tpm.1[row.names(genes.coding.ri.tpm.1)==my.gene,group2]),
          na.rm=T)
      )

      my.ri.frac = c(mean(unlist(genes.coding.ri.frac.flt[
        row.names(genes.coding.ri.frac.flt)==my.gene,group1]), na.rm = T),
        mean(unlist(genes.coding.ri.frac.flt[
          row.names(genes.coding.ri.frac.flt)==my.gene,group2]), na.rm = T)
      )

      my.ratio = my.ri.frac[1] / my.ri.frac[2]
      my.fc = my.ratio
      if(!is.na(my.ratio) & my.ratio < 1) {
        my.fc = -1 * (1 / my.ratio)
      }

      genes.coding.ri.frac.flt.fc = rbind(genes.coding.ri.frac.flt.fc,
                        c(my.gene, my.kept, my.ri.tpm, my.ri.frac, my.ratio, my.fc))

    }

    row.names(genes.coding.ri.frac.flt.logit.stats) = my.row.names

    colnames(genes.coding.ri.frac.flt.fc) = c(
      'Gene',
      'Group1.sm.kept', 'Group2.sm.kept',
      'Group1.ri.meanTPM', 'Group2.ri.meanTPM',
      'Group1.ri.frac', 'Group2.ri.frac',
      'Group1vs2.ratio', 'Group1vs2.FC'
    )
    genes.coding.ri.frac.flt.fc = data.frame(genes.coding.ri.frac.flt.fc)
    row.names(genes.coding.ri.frac.flt.fc) = genes.coding.ri.frac.flt.fc$Gene

    genes.coding.ri.frac.flt.logit.stats$p.adj =
      p.adjust(genes.coding.ri.frac.flt.logit.stats$p.value, method = 'fdr')

    dim(genes.coding.ri.frac.flt.logit.stats) ## 6207
    genes.coding.ri.frac.flt.logit.stats = merge(
      genes.coding.ri.frac.flt.logit.stats,
      genes.coding.ri.frac.flt.fc,
      by = 'row.names'
    )
    dim(genes.coding.ri.frac.flt.logit.stats) ## 6207
    row.names(genes.coding.ri.frac.flt.logit.stats) =
      genes.coding.ri.frac.flt.logit.stats$Row.names
    genes.coding.ri.frac.flt.logit.stats = genes.coding.ri.frac.flt.logit.stats[,-1]

    sum(genes.coding.ri.frac.flt.logit.stats$p.value<0.05) ## 211
    sum(genes.coding.ri.frac.flt.logit.stats$p.adj<0.05) ## 0

    data.plot = genes.coding.ri.frac.flt.logit.stats
    p5 = ggplot(data.plot, aes(p.value)) +
      geom_density()

    gene.anno =read.delim('/gpfs/data/pitroda-lab/ReferenceData/kallisto_indexes/GRCh38.primary_Gencode28_slim_maskPAR/gencode.v28.primary_assembly.annotation.maskPAR.gtf.geneinfo',
                          stringsAsFactors = F)
    row.names(gene.anno) = gene.anno$gene_id

    dim(genes.coding.ri.frac.flt.logit.stats)
    genes.coding.ri.frac.flt.logit.stats =
      merge(gene.anno[,c('gene_id','gene_type','gene_name')],
            genes.coding.ri.frac.flt.logit.stats, by = 'row.names')
    dim(genes.coding.ri.frac.flt.logit.stats)

    row.names(genes.coding.ri.frac.flt.logit.stats) = genes.coding.ri.frac.flt.logit.stats[,1]
    genes.coding.ri.frac.flt.logit.stats = genes.coding.ri.frac.flt.logit.stats[,-1]

    # genes.coding.ri.frac.flt.logit.stats[
    #   genes.coding.ri.frac.flt.logit.stats$gene_name=='TAP1',]
    # x = genes.coding.ri.frac.flt[row.names(genes.coding.ri.frac.flt)=='ENSG00000168394.11',
    #                              group1]
    # y = genes.coding.ri.frac.flt[row.names(genes.coding.ri.frac.flt)=='ENSG00000168394.11',
    #                              group2]
    # mean(unlist(x), na.rm = T)

    genes.coding.ri.frac.flt.logit.stats$change.in.direction = NA
    genes.coding.ri.frac.flt.logit.stats$change.in.direction[
      which(genes.coding.ri.frac.flt.logit.stats$estimate1>
              genes.coding.ri.frac.flt.logit.stats$estimate2)
    ] = 'higher in NR'
    genes.coding.ri.frac.flt.logit.stats$change.in.direction[
      which(genes.coding.ri.frac.flt.logit.stats$estimate1<
              genes.coding.ri.frac.flt.logit.stats$estimate2)
      ] = 'higher in R'

    write.csv(genes.coding.ri.frac.flt.logit.stats,
              file = paste0(rnaseq.file,sample.ex.string, '.genes.coding.ri.frac.sm',
                            ncol(genes.coding.ri.frac),out.string,'.flt.logit.stats.csv'))

    write.table(genes.coding.ri.frac.flt.logit.stats,
              file = paste0(rnaseq.file,sample.ex.string, '.genes.coding.ri.frac.sm',
                            ncol(genes.coding.ri.frac),out.string,'.flt.logit.stats.tsv'),
              sep = '\t', col.names = T, row.names = F, quote = F)

  }

  table(sampleinfo$BORR)
  # NR  R
  # 78 20

}

## -----------------------------------------------------





