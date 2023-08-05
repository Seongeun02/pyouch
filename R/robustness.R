## install ouch package
library(ouch)


# pruneTree --------------------------------------------------------------
# prune tree data to get the same terminal nodes with expression data
pruneTree = function(tree, subData) {
  species = unique(as.character(tree$species[!is.na(tree$species)]))
  curTree = tree
  ##edit tree##
  for (curSpecies in species) {
    if (is.na(subData[curSpecies])) {
      curNode = curTree[curTree$species==curSpecies & !is.na(curTree$species),1]
      curParent = curTree[curTree$species==curSpecies & !is.na(curTree$species),3]
      curGrandparent = curTree[curTree$node==curParent & !is.na(curTree$node),3]
      curSib = curTree[curTree$ancestor==curParent & !is.na(curTree$ancestor) & curTree$node != curNode,1]     
      if (is.na(curTree[curTree$node == curParent,3])) {
        #remove opossum
        curTree = curTree[-which(curTree$node ==curNode),]
        #remove 1
        curTree = curTree[-which(curTree$node ==curParent),]
        #make opossum sib terminal node
        curTree[curTree$node == curSib,3] = NA
        curTree$time = curTree$time - curTree[curTree$node == curSib,4]
      } else {
        #remove terminal node
        curTree = curTree[-which(curTree$node ==curNode),]
        #connect sib to grandparent
        curTree[curTree$node == curSib,3] = curGrandparent
        #remove parent node
        curTree = curTree[-which(curTree$node ==curParent),]
      }
    }
  }
  tree.tr = seq(1, nrow(curTree))
  names(tree.tr) = curTree$node
  tree.tr = as.data.frame(tree.tr)
  curTree$node = tree.tr[as.character(curTree$node),]
  curTree$ancestor = tree.tr[as.character(curTree$ancestor),]
  rownames(curTree) = curTree$node
  return(curTree)
}



# robustness test -------------------------------------------------------------------
robustness = function(organ, ouch, dataPath, outPath,  num_sim = 10) {
  # tree
  mammals.tree = ouchtree(ouch$node, ouch$ancestor, ouch$time, labels = as.character(ouch$species))

  ## species matrix
  exp <- read.delim(paste0(dataPath, as.character(organ[1]),".txt"),
                    header=T, comment.char="", row.names=1, stringsAsFactors = F)
  species <- colnames(exp)
  sp_mat <- c()
  for (i in seq(1,10)){
    sp <- sample(species, length(species)-3)
    sp_mat <- cbind(sp_mat, sp)
  }
  write.table(sp_mat, file=paste0(outPath,"species.txt"), sep="\t", quote=F, row.names=T, col.names=T) 

  for (o in organ){
    out <- paste0(outPath,as.character(o),"/")
  
    exp <- read.delim(paste0(dataPath, as.character(o),".txt"),
                      header=T, comment.char="", row.names=1, stringsAsFactors = F)

    write.table(exp, file=paste0(outPath, as.character(o), "/original.txt"), sep="\t", quote=F, row.names=T, col.names=T) 
    
    ## fit model ---------------------------------------------------------
    for (m in seq(1,10)){
      sp <- sp_mat[,m]
      data_mean.norm <- exp
      
      for (r in seq(length(sp))) {
        data_mean.norm[,sp[r]] <- NA
        
        sigmas = rep(NA, nrow(data_mean.norm))
        alphas = rep(NA, nrow(data_mean.norm))
        logLik = rep(NA, nrow(data_mean.norm))
        pvalues = rep(NA, nrow(data_mean.norm))
        thetas = rep(NA, nrow(data_mean.norm))
        brownSigmas = rep(NA, nrow(data_mean.norm))
        
        cat("\nFitting...\n")
        for (i in seq(1, nrow(data_mean.norm))) {
          if (i %% 1000 == 0) cat(paste0(i, "\tgenes\r"))
          subData = data_mean.norm[i,intersect(colnames(data_mean.norm), as.character(ouch$species[!is.na(ouch$species)]))]
          
          if (sum(!is.na(subData)) < 3) next
          if (sd(subData, na.rm=T) == 0) next
          curTree = pruneTree(ouch, subData)
          curOuchTree = ouchtree(curTree$node, curTree$ancestor, curTree$time, labels = as.character(curTree$species))
          
          #rearrange subData
          subData = subData[as.character(curTree$species[!is.na(curTree$species)])]
          #add ancestral nodes
          curData = c(rep(NA, sum(is.na(curTree$species))), subData[!is.na(subData)])
          names(curData) = curTree$node
          curRegime = as.factor(rep("ns", nrow(curTree)))
          names(curRegime) = curTree$node
          
          h = hansen(data=curData, tree=curOuchTree, regimes=curRegime, sqrt.alpha = 1, sigma = 1, reltol=1e-5)
          b = brown(curData, curOuchTree)
          brownSigmas[i] = coef(b)$sigma.sq.matrix[,1]
          sigmas[i] = coef(h)$sigma.sq.matrix[,1]
          alphas[i] = coef(h)$alpha.matrix[,1]
          thetas[i] = coef(h)$theta$curData
          logLik[i] = logLik(h)
          pvalues[i] =  1 - pchisq((logLik(h) - logLik(b)) * 2, 1)
        }
        #print(head(sigmas))   
        names(sigmas) = rownames(data_mean.norm)
        names(alphas) = rownames(data_mean.norm)
        names(logLik) = rownames(data_mean.norm)
        names(pvalues) = rownames(data_mean.norm)  
        names(thetas) = rownames(data_mean.norm)
        names(brownSigmas) = rownames(data_mean.norm)
        qvalues = p.adjust(pvalues, method="fdr")
        var = sigmas / (alphas*2)
        gene_name = rownames(exp)
        
        stats = cbind(gene_name, qvalues, thetas, var, brownSigmas)
        write.table(stats, file=paste0(out,as.character(m),"_sim/", sp[r], ".txt"), sep="\t", quote=F, row.names=T, col.names=T) 
}
}}
}


# input data --------------------------------------------------------------
organ <- c("liver", "heart", "testis", "kidney", "brain")
ouch <- read.delim("/global/home/users/lizsarah/processed_chen/mammals.tree.txt")
dataPath <- "/global/home/users/lizsarah/processed_chen/mean/"
outPath <- "/global/home/users/lizsarah/rob_test/"


## run
robustness(organ, ouch, dataPath, outPath,  num_sim = 10) 
