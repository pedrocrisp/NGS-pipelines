# Created by Peter Crisp
# August 2014
# pedrocrisp at gmail.com

# For use in differential gene expression analysis.

##### venn pieces #####
  
  VennPieces_2_sample <- function(test.results, venn.results, comparison1, comparison2, out.dir){
    
    
    pdf(paste0(out.dir,"/Venn.pdf"))
    vennDiagram(venn.results, 
                include=c("up","down"), 
                counts.col=c('red', 'blue'), 
                show.include=T, cex = c(0.75))
    dev.off()
    #up
    c1 <- test.results[venn.results[,1]==1 & venn.results[,2]!=1,]
    write.csv(c1, paste0(out.dir,"/",comparison1,"_up.csv"))
    c1.c2 <- test.results[venn.results[,1]==1 & venn.results[,2]==1,]
    write.csv(c1.c2, paste0(out.dir,"/",comparison1,".",comparison2,"_up.csv"))
    c2 <- test.results[venn.results[,1]!=1 & venn.results[,2]==1,]
    write.csv(c2, paste0(out.dir,"/",comparison2,"_up.csv"))
    
    #down
    c1 <- test.results[venn.results[,1]==-1 & venn.results[,2]!=-1,]
    write.csv(c1, paste0(out.dir,"/",comparison1,"_down.csv"))
    c1.c2 <- test.results[venn.results[,1]==-1 & venn.results[,2]==-1,]
    write.csv(c1.c2, paste0(out.dir,"/",comparison1,".",comparison2,"_down.csv"))
    c2 <- test.results[venn.results[,1]!=-1 & venn.results[,2]==-1,]
    write.csv(c2, paste0(out.dir,"/",comparison2,"_down.csv"))
  }
