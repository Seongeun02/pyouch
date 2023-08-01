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


# fitOUmodel --------------------------------------------------------------
fitOUmodel = function(exp, ouch, outPath){

mammals.tree = ouchtree(ouch$node, ouch$ancestor, ouch$time, labels = as.character(ouch$species))
#plot(mammals.tree)

# parameters
sigmas = rep(NA, nrow(exp))
alphas = rep(NA, nrow(exp))
logLik = rep(NA, nrow(exp))
pvalues = rep(NA, nrow(exp))
thetas = rep(NA, nrow(exp))
brownSigmas = rep(NA, nrow(exp))
    
cat("\nFitting...\n")
for (i in seq(1, nrow(exp))) {
  if (i %% 1000 == 0) cat(paste0(i, "\tgenes\r"))
  subData = exp[i,intersect(colnames(exp), as.character(ouch$species[!is.na(ouch$species)]))]
      
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
  
  tryCatch(
    {
      h = hansen(data = curData, tree = curOuchTree, regimes = curRegime, sqrt.alpha = 1, sigma = 1, reltol = 1e-5)
      b = brown(curData, curOuchTree)
      brownSigmas[i] = coef(b)$sigma.sq.matrix[, 1]
      sigmas[i] = coef(h)$sigma.sq.matrix[, 1]
      alphas[i] = coef(h)$alpha.matrix[, 1]
      thetas[i] = coef(h)$theta$curData
      logLik[i] = logLik(h)
      pvalues[i] = 1 - pchisq((logLik(h) - logLik(b)) * 2, 1)
    },
    error = function(e) {
      # Unsuccessful convergence or error during optimization
      # Handle the failure or try alternative approaches
      print("Optimization did not converge.")
      brownSigmas[i] = NA
      sigmas[i] = NA
      alphas[i] = NA
      thetas[i] = NA
      logLik[i] = NA
      pvalues[i] = NA
    }
  )
}

names(sigmas) = rownames(exp)
names(alphas) = rownames(exp)
names(logLik) = rownames(exp)
names(pvalues) = rownames(exp)  
names(thetas) = rownames(exp)
names(brownSigmas) = rownames(exp)
qvalues = p.adjust(pvalues, method="fdr")
var = sigmas / (alphas*2)
gene_name = rownames(exp)
    
stats = cbind(gene_name, qvalues, sigmas, alphas, thetas, var, brownSigmas)
write.table(stats, file=outPath, sep="\t", quote=F, row.names=T, col.names=T)
return(stats)
}


# input -------------------------------------------------------------------
## input > expression data should be mean data for each species
exp <- read.delim("/global/home/users/lizsarah/processed_chen/mean/liver.txt",
                  header=T, comment.char="", row.names=1, stringsAsFactors = F)
ouch <- read.delim("/global/home/users/lizsarah/processed_chen/mammals.tree.txt")
outPath <- "/global/home/users/lizsarah/test.txt"


# fitting
stats <- fitOUmodel(exp, ouch, outPath)

