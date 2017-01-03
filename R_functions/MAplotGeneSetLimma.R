# Created by Peter Crisp
# August 2014
# pedrocrisp at gmail.com

# For use in differential gene expression analysis.

##### MAplotGeneSet function #####
  
  #for debugging
  # inputList = "xrn2xrn3.vs.control"
  # inputListName = "xrn2xrn3"
  # geneSetList = "alx8.vs.control"
  # geneSetListName = "alx8"
  # MArrayLMobject = fit2
  
  MAplotGeneSetLimma <- function(MArrayLMobject, resultsMatrix, inputList, inputListName, geneSetList, geneSetListName){ 
    plotMA(MArrayLMobject, coef=inputList, cex=0.5, status=resultsMatrix[,geneSetList], values=c(1, -1), col=c("red", "blue"), main=paste0(inputList), legend=F)
    abline(h=0,col="darkgrey")
    legend("bottomright", legend=(c(paste0("up in ", geneSetListName), paste0("down in ", geneSetListName))), bty =  "n", pch=c(16, 16), col = c("red", "blue"), text.col = c("red", "blue"))
  }

MAplotGeneSetGreenLimma <- function(MArrayLMobject, resultsMatrix, inputList, inputListName, geneSetList, geneSetListName){ 
  plotMA(MArrayLMobject, coef=inputList, cex=0.5, status=resultsMatrix[,geneSetList], values=c(1, -1), col=c("red", "green"), main=paste0(inputList), legend=F)
  abline(h=0,col="darkgrey")
  legend("bottomright", legend=(c(paste0("up in ", geneSetListName), paste0("down in ", geneSetListName))), bty =  "n", pch=c(16, 16), col = c("red", "green"), text.col = c("red", "green"))
}

MAplotGeneSet2Limma <- function(MArrayLMobject, geneSetResultsMatrix, DEcallNames, inputList, inputListName, geneSetListName){ 
  plotMA(MArrayLMobject, coef=inputList, cex=0.35, status=geneSetResultsMatrix, values=DEcallNames, col=c("#00FF00",
                                                                                                        "#FF0000"), main=paste0(inputListName), legend=F)
  abline(h=0,col="darkgrey")
  legend("bottomright", cex = 0.5, legend=(c(paste0("up in ", geneSetListName), paste0("down in ", geneSetListName))), bty =  "n", pch=c(16, 16), col = c("#00FF00",
                                                                                                                                                          "#FF0000"), text.col = c("black", "black"))
}

MAplotGeneSet2rLimma <- function(MArrayLMobject, geneSetResultsMatrix, DEcallNames, inputList, inputListName, geneSetListName){ 
  plotMA(MArrayLMobject, coef=inputList, cex=0.35, status=geneSetResultsMatrix, values=DEcallNames, col=c("red",
                                                                                                        "blue"), main=paste0(inputListName), legend=F)
  abline(h=0,col="darkgrey")
  legend("bottomright", cex = 0.5, legend=(c(paste0("up in ", geneSetListName), paste0("down in ", geneSetListName))), bty =  "n", pch=c(16, 16), col = c("red",
                                                                                                                                                          "blue"), text.col = c("black", "black"))
}

MAplotGeneSet2gLimma <- function(MArrayLMobject, geneSetResultsMatrix, DEcallNames, inputList, inputListName, geneSetListName){ 
  plotMA(MArrayLMobject, coef=inputList, cex=0.35, status=geneSetResultsMatrix, values=DEcallNames, col=c("red",
                                                                                                        "green"), main=paste0(inputListName), legend=F)
  abline(h=0,col="darkgrey")
  legend("bottomright", cex = 0.5, legend=(c(paste0("up in ", geneSetListName), paste0("down in ", geneSetListName))), bty =  "n", pch=c(16, 16), col = c("red",
                                                                                                                                                          "green"), text.col = c("black", "black"))
}

MAplotGeneSet3Limma <- function(MArrayLMobject, geneSetResultsMatrix, DEcallNames, inputList, inputListName, geneSetListName){ 
  plotMA(MArrayLMobject, coef=inputList, cex=0.75, status=geneSetResultsMatrix, values=DEcallNames, col=c("#FF0000", "#00FF00", "#507DC3"), main=paste0(inputListName), legend=F)
  abline(h=0,col="darkgrey")
  legend("bottomright", cex = 0.5, legend=(c(paste0("up in ", geneSetListName), c(paste0("unchanged in ", geneSetListName)), paste0("down in ", geneSetListName))), bty =  "n", pch=c(16, 16), col = c("#FF0000", "#00FF00", "#507DC3"), text.col = c("black", "black", "black"))
}
