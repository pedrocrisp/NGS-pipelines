#!/usr/bin/Rscript
##########
library(fields)

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
Sample <- args[1]
Sample
beds_folder <- args[2]
beds_folder
outDir <- args[3]
outDir
library_layout <- args[4]
library_layout
#minLength <- as.numeric(args[4])
#minLength
#trimLength <- as.numeric(args[5])
#trimLength
#minCoverage <- as.numeric(args[6])
#minCoverage
#Sample <- c("alx8_277_7")
sPath <- paste0(beds_folder, "/", Sample, "/")
outFolder <- paste0(outDir, "/", Sample)
dir.create(outFolder, showWarnings = F, recursive = T)
#####

if (library_layout == "nonstranded") {
  #read in file
  #debug
  # Sample="Sample_alx8_277_9"
  # plus.input=read.delim("Sample_alx8_277_9/Sample_alx8_277_9.plus.dist.1k.bed", head=F)
  plus.input=read.delim(paste0(sPath, Sample, ".dist.1k.bed"),head=F)
  #develop strand directional positioning
  real.dist=matrix(ifelse(plus.input[,10]=='+',-1*plus.input[,17],plus.input[,17]),ncol=1)
  plus.input=cbind(plus.input,real.dist)

  #create relative distance measure from coverage data start and entire gene (not just exon or intron) start stop.
  rel.dist=matrix(ifelse(plus.input$real.dist==0,
  ifelse(plus.input[,10]=="-",
  ((plus.input[,15] - (plus.input[,2]))/(plus.input[,15] - plus.input[,14]))*1000,
  (((plus.input[,2]) - plus.input[,14])/(plus.input[,15] - plus.input[,14]))*1000),
  ifelse(plus.input$real.dist>0,
  plus.input$real.dist + 1000,
  plus.input$real.dist)),
  ncol=1)
  plus.input=cbind(plus.input,rel.dist)

  #subset to exons only
  plus.exon=subset(plus.input,plus.input$V12=='exon')

  #look to subset coverage based on primary strand and input coverage information. Primary = matched, offstrand should be non genic strand
  plus.exon.primary=subset(plus.exon,plus.exon$V10=='+')
  plus.exon.offstrand=subset(plus.exon,plus.exon$V10=='-')

  #stats bin in 300 bins for both promary and offstrand
  plus.exon.primary.bin=stats.bin(plus.exon.primary$rel.dist,log(abs(plus.exon.primary[,4])+1),N=300)
  pepb=cbind(matrix(plus.exon.primary.bin$centers,ncol=1),plus.exon.primary.bin$stats["mean",])

  plus.exon.offstrand.bin=stats.bin(plus.exon.offstrand$rel.dist,log(abs(plus.exon.offstrand[,4])+1),N=300)
  peob=cbind(matrix(plus.exon.offstrand.bin$centers,ncol=1),plus.exon.offstrand.bin$stats["mean",])

  peob[,2]=-peob[,2]

# i think this bit is wrong...
  sense=pepb
  sense[,2]=sense[,2]+peob[,2]

  antisense=peob
  antisense[,2]=antisense[,2]+pepb[,2]

  out.table <- data.frame(cbind(sense,antisense[,2]))
  colnames(out.table) <- c("Position", "Sense", "Antisense")
  out.table$SampleName <- Sample

  #write.csv(out.table, 'average_coverage.csv')
  write.csv(out.table, paste0(outFolder, "/",Sample, '_average_coverage.csv'))

  pdf(paste0(outFolder, "/",Sample, '_gene_coverge_plot_WC.pdf'),h=10,w=12)
  #pdf("test_strands.pdf",h=10,w=12)
  plot(x=NULL,y=NULL,xlim=c(-1000,2000),ylim=c(-5,5), main=Sample)
  lines(pepb,col=1,lwd=2)
  lines(peob,col=2,lwd=2)
  lines(mepb,col=3,lwd=2)
  lines(meob,col=4,lwd=2)

  abline(v=0,lty=2)
  abline(v=1000,lty=2)
  abline(h=0,lty=1,col='grey')
  legend('topright',c('plus_primary','plus_offstrand','minus_primary','minus_offstrand'),lty=1,col=c(1,2,3,4))
  dev.off()

  pdf(paste0(outFolder, "/",Sample, '_gene_coverge_plot.pdf'),h=10,w=12)
  #pdf("test_sense.pdf",h=10,w=12)
  plot(x=NULL,y=NULL,xlim=c(-1000,2000),ylim=c(-5,5),main=Sample)
  lines(sense,col=1,lwd=2)
  lines(antisense,col=2,lwd=2)

  abline(v=0,lty=2)
  abline(v=1000,lty=2)
  abline(h=0,lty=1,col='grey')
  legend('topright',c('sense','antisense'),lty=1,col=c(1,2))
  dev.off()

}else{

#read in file
#debug
# Sample="Sample_alx8_277_9"
# plus.input=read.delim("Sample_alx8_277_9/Sample_alx8_277_9.plus.dist.1k.bed", head=F)
plus.input=read.delim(paste0(sPath, Sample, ".plus.dist.1k.bed"),head=F)
#develop strand directional positioning
real.dist=matrix(ifelse(plus.input[,10]=='+',-1*plus.input[,17],plus.input[,17]),ncol=1)
plus.input=cbind(plus.input,real.dist)

#create relative distance measure from coverage data start and entire gene (not just exon or intron) start stop.
rel.dist=matrix(ifelse(plus.input$real.dist==0,
ifelse(plus.input[,10]=="-",
((plus.input[,15] - (plus.input[,2]))/(plus.input[,15] - plus.input[,14]))*1000,
(((plus.input[,2]) - plus.input[,14])/(plus.input[,15] - plus.input[,14]))*1000),
ifelse(plus.input$real.dist>0,
plus.input$real.dist + 1000,
plus.input$real.dist)),
ncol=1)
plus.input=cbind(plus.input,rel.dist)

#subset to exons only
plus.exon=subset(plus.input,plus.input$V12=='exon')

#look to subset coverage based on primary strand and input coverage information. Primary = matched, offstrand should be non genic strand
plus.exon.primary=subset(plus.exon,plus.exon$V10=='+')
plus.exon.offstrand=subset(plus.exon,plus.exon$V10=='-')

#stats bin in 300 bins for both promary and offstrand
plus.exon.primary.bin=stats.bin(plus.exon.primary$rel.dist,log(abs(plus.exon.primary[,4])+1),N=300)
pepb=cbind(matrix(plus.exon.primary.bin$centers,ncol=1),plus.exon.primary.bin$stats["mean",])

plus.exon.offstrand.bin=stats.bin(plus.exon.offstrand$rel.dist,log(abs(plus.exon.offstrand[,4])+1),N=300)
peob=cbind(matrix(plus.exon.offstrand.bin$centers,ncol=1),plus.exon.offstrand.bin$stats["mean",])


#read in file
#debug
# minus.input=read.delim("Sample_alx8_277_9/Sample_alx8_277_9.minus.dist.1k.bed", head=F)
minus.input=read.delim(paste0(sPath, Sample, ".minus.dist.1k.bed"),head=F)
#develop strand directional positioning
real.dist=matrix(ifelse(minus.input[,10]=='+',-1*minus.input[,17],minus.input[,17]),ncol=1)
minus.input=cbind(minus.input,real.dist)

#create relative distance measure from coverage data start and entire gene (not just exon or intron) start stop.
rel.dist=matrix(ifelse(minus.input$real.dist==0,
ifelse(minus.input[,10]=="-",
((minus.input[,15] - (minus.input[,2]))/(minus.input[,15] - minus.input[,14]))*1000,
(((minus.input[,2]) - minus.input[,14])/(minus.input[,15] - minus.input[,14]))*1000),
ifelse(minus.input$real.dist>0,
minus.input$real.dist + 1000,
minus.input$real.dist)),
ncol=1)
minus.input=cbind(minus.input,rel.dist)

#subset to exons only
minus.exon=subset(minus.input,minus.input$V12=='exon')

#look to subset coverage based on primary strand and input coverage information. Primary = matched, offstrand should be non genic strand
minus.exon.primary=subset(minus.exon,minus.exon$V10=='-')
minus.exon.offstrand=subset(minus.exon,minus.exon$V10=='+')

#stats bin in 300 bins for both promary and offstrand using log transformed values
minus.exon.primary.bin=stats.bin(minus.exon.primary$rel.dist,log(abs(minus.exon.primary[,4])+1),N=300)
mepb=cbind(matrix(minus.exon.primary.bin$centers,ncol=1),minus.exon.primary.bin$stats["mean",])

minus.exon.offstrand.bin=stats.bin(minus.exon.offstrand$rel.dist,log(abs(minus.exon.offstrand[,4])+1),N=300)
meob=cbind(matrix(minus.exon.offstrand.bin$centers,ncol=1),minus.exon.offstrand.bin$stats["mean",])

peob[,2]=-peob[,2]
meob[,2]=-meob[,2]

sense=pepb
sense[,2]=sense[,2]+mepb[,2]

antisense=peob
antisense[,2]=antisense[,2]+meob[,2]

out.table <- data.frame(cbind(sense,antisense[,2]))
colnames(out.table) <- c("Position", "Sense", "Antisense")
out.table$SampleName <- Sample

#write.csv(out.table, 'average_coverage.csv')
write.csv(out.table, paste0(outFolder, "/",Sample, '_average_coverage.csv'))

pdf(paste0(outFolder, "/",Sample, '_gene_coverge_plot_WC.pdf'),h=10,w=12)
#pdf("test_strands.pdf",h=10,w=12)
plot(x=NULL,y=NULL,xlim=c(-1000,2000),ylim=c(-5,5), main=Sample)
lines(pepb,col=1,lwd=2)
lines(peob,col=2,lwd=2)
lines(mepb,col=3,lwd=2)
lines(meob,col=4,lwd=2)

abline(v=0,lty=2)
abline(v=1000,lty=2)
abline(h=0,lty=1,col='grey')
legend('topright',c('plus_primary','plus_offstrand','minus_primary','minus_offstrand'),lty=1,col=c(1,2,3,4))
dev.off()

pdf(paste0(outFolder, "/",Sample, '_gene_coverge_plot.pdf'),h=10,w=12)
#pdf("test_sense.pdf",h=10,w=12)
plot(x=NULL,y=NULL,xlim=c(-1000,2000),ylim=c(-5,5),main=Sample)
lines(sense,col=1,lwd=2)
lines(antisense,col=2,lwd=2)

abline(v=0,lty=2)
abline(v=1000,lty=2)
abline(h=0,lty=1,col='grey')
legend('topright',c('sense','antisense'),lty=1,col=c(1,2))
dev.off()

}
#write.csv(out.transposed.percentage, paste0(outFolder, "/",Sample, '_mRNA_densities_2bins_percentage.csv'))
