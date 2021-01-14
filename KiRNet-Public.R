# Uncomment to set working directory
# setwd()

load("KiRNet_V3_Initializers")

##########################
# KiR Elastic Net models #
##########################
library(glmnet)
library(ggplot2)
library(tictoc)
library(readxl)
detach("package:dplyr")

drugTable <- as.data.frame(read_excel("drug_profiles_public.xlsx"))
kinaseIDs <- colnames(drugTable)[-1]
drugIDs <- drugTable$Compound
row.names(drugTable) <- drugIDs
drugTable <- drugTable[,-1]

kirData <- as.data.frame(read_excel("Huh7_allData_public.xlsx", sheet = "KiRdrugResponses"))
response <- kirData$`Huh7-Fzd2_Response`

controlIdx <- grep('^control|DMSO|empty|vehicle$', kirData[,1], ignore.case = TRUE) # Check for a control in data, used to normalize and/or include as full kinase activity
if(length(controlIdx)==1){
  control <- kirData[controlIdx,2]
  kirData[controlIdx,1] <- 'Control'
  response <- response / control}
inputs <- drugTable[match(kirData[,1],drugIDs),]
rm(controlIdx, control)

standardize <- FALSE        # Should x values (inhibitions) be standardized
controlWeight <- 1          # Relative weight to assign to control data point. Set to 0 to exclude, 1 for equal weighting.
alphaVals <- seq(from = 0.0, to = 1.0, by = 0.1)         # Sequence of alpha values (between 0 and 1) to use

numLambdas <- 300
lambdaMin <- 0.005
nfolds <- length(response)          # For LOOCV: length(response)     otherwise, usually 5 or 10
obsWeights <- ifelse(kirData[,1] == 'Control', controlWeight, 1)
varPenalty <- rep(1, ncol(inputs))
keep <- TRUE
algorithmType <- 'naive'    # Can be 'naive' or 'covariance'


# Model Construction

alphaLabels <- paste0('alpha = ', alphaVals)

# Run elastic net with the specified parameters
glmobjs <- lapply(alphaVals, function(x) {
  glmobj <- cv.glmnet(x = as.matrix(inputs), y = response,
                      family = 'gaussian', weights = obsWeights, alpha = x,
                      nlambda = numLambdas, lambda.min.ratio = lambdaMin, 
                      standardize = standardize, intercept = TRUE, nfolds = nfolds, 
                      penalty.factor = varPenalty, type.gaussian=algorithmType, keep = keep, grouped = FALSE)
  glmobj$s <- glmobj$lambda.min
  glmobj
})
names(glmobjs) <- alphaVals

#Outputs are formatted as dataframes
model <- sapply(glmobjs, function(x) as.matrix(coef(x, s = x$s)))
dimnames(model) <- list(c('Intercept', kinaseIDs), alphaLabels)
model <- model[do.call(order, as.data.frame(-model[,ncol(model):1])),]

preds <- sapply(glmobjs, function(x) predict(x, as.matrix(drugTable), s = x$s))
preds <- cbind(preds, apply(preds, 1, mean))
dimnames(preds) <- list(drugIDs, c(alphaLabels, 'Mean'))
preds <- preds[order(preds[,'Mean']),]

modelProps <- sapply(glmobjs, function(x) {
  lambda = x$s
  df <- x$glmnet.fit$df[which(x$lambda == lambda)]
  MSE <- x$cvm[which(x$lambda == lambda)]
  RMSE <- sqrt(MSE)
  CVRMSE <- sqrt(MSE) / diff(range(response))
  c(x$lambda.min, lambda, df, MSE, RMSE, CVRMSE)
})
dimnames(modelProps) <- list(c('Lambda Min', 'Lambda', 'Number of Variables', 'CrossValMSE', 'RMSE', 'CVRMSE'), alphaLabels)

rm(keep, controlWeight, lambdaMin, numLambdas, obsWeights, standardize, varPenalty, nfolds, algorithmType)

#Identify key kinases
keyKinases <- row.names(model[rowMeans(model[,-c(1)]) > 0,])
keyKinases <- keyKinases[keyKinases != "Intercept"]
keyKinaseInfo <- kinaseInfo[kinaseInfo$KIR %in% keyKinases, c(4,1:2,5:7)]


#############
# KiRNet v3 #
#############
library(tictoc)
library(tidyr)
library(ggplot2)
library(purrr)
library(dplyr)
library(PCSF)
library(igraph)
library(readxl)
library(pheatmap)
library(colorspace)

#############
# Functions #
#############

# Function to plot current state of the network
plotCurrentNetwork <- function(DFset, layout = layout_nicely, nodeColors = NA,
                               directOnly = TRUE, expressedOnly = TRUE, printPlot = TRUE){
  # Generate appropriate graph from node and interaction dataframes
  graphObj <- graph_from_data_frame(d = DFset$Interactions, directed = TRUE, vertices = DFset$Nodes)
  if(directOnly){graphObj <- delete.edges(graphObj, E(graphObj)[(E(graphObj)$SUBTYPE1 %in% "indirect effect") | (E(graphObj)$SUBTYPE2 %in% "indirect effect")])}
  if(expressedOnly){graphObj <- delete.vertices(graphObj, V(graphObj)[(V(graphObj)$Expressed == -1)])}
  
  # Add vertex aesthetics based on gene/protein properties
  V(graphObj)$shape <- "circle"
  V(graphObj)$shape[V(graphObj)$Kinase] <- "square"
  V(graphObj)$color <- "gray75"
  V(graphObj)$color[V(graphObj)$PhosphoSig == 1] <- "gold"
  V(graphObj)$color[V(graphObj)$PhosphoSig == -1] <- "deepskyblue"
  V(graphObj)$color[V(graphObj)$InKeyFunctional] <- "darkorchid"

  if(length(nodeColors) == length(V(graphObj))) V(graphObj)$color <- nodeColors
  
  V(graphObj)$label.color <- "black"
  V(graphObj)$label.color[V(graphObj)$InKeyFunctional] <- "white"
  V(graphObj)$label.family <- "sans"
  V(graphObj)$label.font <- 1
  V(graphObj)$label.cex <- 0.667
  # V(graphObj)$label.cex[V(graphObj)$InKeyFunctional] <- 0.87
  # V(graphObj)$label <- ifelse(V(graphObj)$InKeyFunctional | V(graphObj)$PhosphoSig != 0, V(graphObj)$name, NA)
  V(graphObj)$label <- V(graphObj)$name
  # V(graphObj)$size <- ifelse(!V(graphObj)$InKeyFunctional & V(graphObj)$PhosphoSig==0, 15, 25)
  V(graphObj)$size <- 30
  
  # Add edge aesthetics based on interaction properties
  E(graphObj)$color <- "gray80"
  E(graphObj)$color[E(graphObj)$SUBTYPE1 %in% c("activation", "expression")] <- "#009E73"
  E(graphObj)$color[E(graphObj)$SUBTYPE1 %in% c("inhibition", "repression")] <- "brown1"
  E(graphObj)$color[E(graphObj)$TYPE == "STRING"] <- "blue"
  E(graphObj)$arrow.size <- 0.2
  E(graphObj)$width <- 1.5
  E(graphObj)$lty <- "solid"
  E(graphObj)$lty[E(graphObj)$TYPE == "GErel"] <- "longdash"
  E(graphObj)$lty[E(graphObj)$TYPE == "PPrel"] <- "solid"
  E(graphObj)$lty[E(graphObj)$TYPE == "STRING"] <- "dotted"
  E(graphObj)$lty[(E(graphObj)$SUBTYPE1 %in% c("indirect effect")) | (E(graphObj)$SUBTYPE2 %in% c("indirect effect"))] <- "dotted"
  
  if(printPlot){plot(graphObj, layout = layout, frame = TRUE, edge.loop.angle = pi)}
  return(graphObj)
}

# Function to reconstruct a list of vertex and edge dataframes from an igraph object
graphToSet <- function(graphObj){
  graphObjSet <- list("Nodes" = igraph::as_data_frame(graphObj, "vertices"), "Interactions" =  igraph::as_data_frame(graphObj, "edges"))
  # Fix variable names to preserve compatibility with other functions/code
  colnames(graphObjSet$Nodes)[1] <- "SYMBOL"
  colnames(graphObjSet$Interactions)[1:2] <- c("FROM.SYMBOL", "TO.SYMBOL")
  #Remove graphing parameters
  graphObjSet$Nodes <- graphObjSet$Nodes[,-which(colnames(graphObjSet$Nodes) %in% c("shape","color","label.family", "label.font", "label.color","size","label.cex", "label"))]
  graphObjSet$Interactions <- graphObjSet$Interactions[, -which(colnames(graphObjSet$Interactions) %in% c("color","arrow.size","width","lty"))]
  # Fix variable types
  graphObjSet$Nodes$Kinase <- as.logical(graphObjSet$Nodes$Kinase)
  graphObjSet$Nodes$InKiR <- as.logical(graphObjSet$Nodes$InKiR)
  graphObjSet$Nodes$InKeyFunctional <- as.logical(graphObjSet$Nodes$InKeyFunctional)
  graphObjSet$Nodes$Mutated <- as.logical(graphObjSet$Nodes$Mutated)
  return(graphObjSet)
}

# Function to expand network by one intermediate DIRECT relation step
expandNetwork <- function(beforeSet){
  beforeNodes <- beforeSet$Nodes$SYMBOL
  directInteractions <- KEGGInteractions[!(KEGGInteractions$SUBTYPE1 %in% "indirect effect") & !(KEGGInteractions$SUBTYPE2 %in% "indirect effect"),]
  toNodes <- directInteractions$TO.SYMBOL[directInteractions$FROM.SYMBOL %in% beforeNodes]
  fromNodes <- directInteractions$FROM.SYMBOL[directInteractions$TO.SYMBOL %in% beforeNodes]
  nodesOfInterest <- union(beforeNodes, union(toNodes, fromNodes))
  afterNodes <- allNodes[match(nodesOfInterest, allNodes$SYMBOL),]
  afterInteractions <- KEGGInteractions[(KEGGInteractions$FROM.SYMBOL %in% afterNodes$SYMBOL) & (KEGGInteractions$TO.SYMBOL %in% afterNodes$SYMBOL), ]
  return(list("Nodes" = afterNodes, "Interactions" =  afterInteractions))
}

# Function to prune nodes/relations disconnected from given nodes
pruneNetwork <- function(beforeSet, alwaysNodes, maxDistance = Inf, directOnly = TRUE, expressedOnly = TRUE){
  beforeNodes <- beforeSet$Nodes
  beforeInteractions <- beforeSet$Interactions
  
  # Remove unwanted nodes/edges
  if(expressedOnly){
    beforeNodes <- beforeNodes[!(beforeNodes$Expressed == -1),]
    beforeInteractions <- beforeInteractions[beforeInteractions$FROM.SYMBOL %in% beforeNodes$SYMBOL & beforeInteractions$TO.SYMBOL %in% beforeNodes$SYMBOL,]}
  if(directOnly){
    beforeInteractions <- beforeInteractions[!(beforeInteractions$SUBTYPE1 %in% "indirect effect" | beforeInteractions$SUBTYPE2 %in% "indirect effect"),]}
  graphObj <- graph_from_data_frame(d = beforeInteractions, directed = TRUE, vertices = beforeNodes)
  
  # Identify nodes connected to seed nodes
  tmp <- distances(graphObj, v = V(graphObj)[V(graphObj)$name %in% alwaysNodes], mode = "all")
  connectedNodes <- colnames(tmp)[apply(tmp, 2, min) <= (maxDistance)]
  
  # Compile final list of nodes to keep
  afterNodes <- beforeNodes[beforeNodes$SYMBOL %in% c(alwaysNodes, connectedNodes),]
  afterInteractions <- beforeInteractions[(beforeInteractions$FROM.SYMBOL %in% afterNodes$SYMBOL) & (beforeInteractions$TO.SYMBOL %in% afterNodes$SYMBOL),]
  return(list("Nodes" = afterNodes, "Interactions" =  afterInteractions))
}

# Function to iteratively apply expandNetwork n times or until saturated (up to 100 expansions)
expandMany <- function(beforeSet, n = 100, directOnly = TRUE, expressedOnly = TRUE){
  tic("Total expansion time")
  expandReport <- data.frame("Num.Expansions" = 0:n,
                             "Num.Nodes" = NA,
                             "Num.Interactions" = NA,
                             "Num.KeyFunctional.Nodes" = NA,
                             stringsAsFactors = FALSE)
  expandReport$Num.Nodes[1] <- nrow(beforeSet$Nodes)
  expandReport$Num.Interactions[1] <- nrow(beforeSet$Interactions)
  expandReport$Num.KeyFunctional.Nodes[1] <- sum(beforeSet$Nodes$InKeyFunctional)
  workingSet <- beforeSet
  for(i in 1:n){
    tic(paste0("Expanding/pruning step ", i))
    workingSet <- expandNetwork(workingSet)
    expandReport$Num.Nodes[i+1] <- nrow(workingSet$Nodes)
    expandReport$Num.Interactions[i+1] <- nrow(workingSet$Interactions)
    expandReport$Num.KeyFunctional.Nodes[i+1] <- sum(workingSet$Nodes$InKeyFunctional)
    toc()
    if((expandReport$Num.Nodes[i+1] == expandReport$Num.Nodes[i]) & 
       (expandReport$Num.Interactions[i+1] == expandReport$Num.Interactions[i])) break
  }
  toc()
  expandReport <- expandReport[!apply(expandReport[,-1], 1, function(x) all(is.na(x))),]
  return(list("Set" = workingSet, "Report" = expandReport))
}

# Function to weight edges based on node degrees, penalizing highly-connected nodes
addGraphWeights <- function(graphToWeight, scale.mean = NA){
  E(graphToWeight)$weight <- log(degree(graphToWeight, v = ends(graphToWeight, E(graphToWeight), names = FALSE)[,1], mode = "out")
                                 *degree(graphToWeight, v = ends(graphToWeight, E(graphToWeight), names = FALSE)[,2], mode = "in"))+1
  if(is.numeric(scale.mean)){E(graphToWeight)$weight <- E(graphToWeight)$weight * scale.mean / mean(E(graphToWeight)$weight)}
  return(graphToWeight)}

# Function to collapse sibling nodes into family nodes
collapseSiblings <- function(beforeGraph){
  # Note: Vertex attributes are chosen via 'max',
  # Edge attributes are chose via 'first'
  allNeighbors <- lapply(V(beforeGraph), function(x) V(beforeGraph)$name[neighbors(beforeGraph, x, mode = "all")]) 
  familyGraph <- contract(beforeGraph, mapping = as.numeric(factor(sapply(allNeighbors, paste0, collapse = '+'),
                                                                   levels = unique(sapply(allNeighbors, paste0, collapse = '+')),
                                                                   ordered = TRUE)),
                          vertex.attr.comb = list(name = function(x) paste0(x, collapse = ', '), "max"))
  familyGraph <- simplify(familyGraph, remove.multiple = TRUE, remove.loops = FALSE, edge.attr.comb = "first")
  return(familyGraph)
}

##################################
# Initiating general interactome #
##################################

load("KiRNet_V3_Initializers")

allNodes <- data.frame(SYMBOL = IDs$SYMBOL[IDs$SYMBOL %in% c(KEGGInteractions$FROM.SYMBOL, KEGGInteractions$TO.SYMBOL)],
                       Kinase = NA,
                       InKiR = NA,
                       InKeyFunctional = NA,
                       Mutated = NA,
                       RNALevel = NA,
                       Expressed = 0,
                       RNAsig = 0,
                       ProteinLevel = NA,
                       ProteinSig = 0,
                       PhosphoSig = 0,
                       stringsAsFactors = FALSE)
allNodes$Kinase <- allNodes$SYMBOL %in% kinaseInfo$SYMBOL
allNodes$InKiR <- allNodes$SYMBOL %in% kinaseInfo$SYMBOL[kinaseInfo$IN.KIR]
generalSet <- list("Nodes" = allNodes, "Interactions" = KEGGInteractions)
rm(allNodes)
generalGraph <- plotCurrentNetwork(generalSet, directOnly = TRUE, expressedOnly = TRUE, printPlot = FALSE)
generalSet <- graphToSet(generalGraph)

##########################################################################################
# Annotate nodes (key functional nodes, RNA levels, DNA mutations, protein levels, etc.) #
##########################################################################################

# Load key functional nodes and continue annotating protein data
keyFunctionalFile <- read_excel("Huh7_allData_public.xlsx", sheet = "keyKinases")
keyFunctionalNames <- keyFunctionalFile$Huh7_Fzd2OE
keyFunctionalNames <- keyFunctionalNames[keyFunctionalNames != "" & !is.na(keyFunctionalNames)]
keyFunctionals <- kinaseInfo$SYMBOL[match(keyFunctionalNames, kinaseInfo$KIR)]
generalSet$Nodes$InKeyFunctional <- generalSet$Nodes$SYMBOL %in% keyFunctionals
rm(keyFunctionalNames, keyFunctionals)

# Load genetic mutations and annotate
mutationsFile <- read_excel("Huh7_allData_public.xlsx", sheet = "geneMutations")
generalSet$Nodes$Mutated <- generalSet$Nodes$SYMBOL %in% mutationsFile$Symbol

# Load RNA Seq data and continue annotating protein data
RNA_cutoff <- 1     # RNA values AT OR BELOW this are considered not expressed.

RNAfile <- read_excel("Huh7_allData_public.xlsx", sheet = "RNAseq")
generalSet$Nodes$RNALevel <- log2(RNAfile$Huh7_Fzd2OE_Average[match(IDs$ENTREZ[match(generalSet$Nodes$SYMBOL, IDs$SYMBOL)], RNAfile$Entrez)] +1)
generalSet$Nodes$Expressed <- if_else(generalSet$Nodes$RNALevel > RNA_cutoff, 1, -1)
generalSet$Nodes$Expressed[is.na(generalSet$Nodes$Expressed)] <- 0
generalSet$Nodes$RNALevel <- (generalSet$Nodes$RNALevel - min(generalSet$Nodes$RNALevel, na.rm = TRUE))
generalSet$Nodes$RNALevel <- generalSet$Nodes$RNALevel / mean(generalSet$Nodes$RNALevel, na.rm = TRUE)
generalSet$Nodes$RNAsig <- if_else(RNAfile$sig[match(IDs$ENTREZ[match(generalSet$Nodes$SYMBOL, IDs$SYMBOL)], RNAfile$Entrez)] == 1,
                                   if_else(RNAfile$logFC[match(IDs$ENTREZ[match(generalSet$Nodes$SYMBOL, IDs$SYMBOL)], RNAfile$Entrez)] > 1, 1, -1), 0)
generalSet$Nodes$RNAsig[is.na(generalSet$Nodes$RNAsig)] <- 0
rm(RNA_cutoff)

# Load protein data
proteinFile <- read_excel("Huh7_allData_public.xlsx", sheet = "proteinSites")
proteinFile <- proteinFile[order(proteinFile$Symbol, -abs(proteinFile$Difference * proteinFile$sig)),]
proteinFile <- proteinFile[!duplicated(proteinFile$Symbol),]
generalSet$Nodes$ProteinLevel <- proteinFile$Huh7_Fzd2OE_average[match(generalSet$Nodes$SYMBOL, proteinFile$Symbol)]
generalSet$Nodes$ProteinLevel <- generalSet$Nodes$ProteinLevel - min(generalSet$Nodes$ProteinLevel, na.rm = TRUE)
generalSet$Nodes$ProteinLevel <- generalSet$Nodes$ProteinLevel / mean(generalSet$Nodes$ProteinLevel, na.rm = TRUE)
generalSet$Nodes$ProteinSig <- proteinFile$sig[match(generalSet$Nodes$SYMBOL, proteinFile$Symbol)]
generalSet$Nodes$ProteinSig[is.na(generalSet$Nodes$ProteinSig)] <- 0

# Load phospho data
phosphoFile <- read_excel("Huh7_allData_public.xlsx", sheet = "phosphoSites")
phosphoFile <- phosphoFile[order(phosphoFile$Symbol, -abs(phosphoFile$Difference * phosphoFile$sig)),]
phosphoFile <- phosphoFile[!duplicated(phosphoFile$Symbol),]
generalSet$Nodes$PhosphoSig <- phosphoFile$sig[match(generalSet$Nodes$SYMBOL, phosphoFile$Symbol)]
generalSet$Nodes$PhosphoSig[is.na(generalSet$Nodes$PhosphoSig)] <- 0

########################################################
# Finalize general network and measure node properties #
########################################################

# Generate system network graph, select only largest component, and weight edges
generalGraph <- plotCurrentNetwork(generalSet, directOnly = TRUE, expressedOnly = TRUE, printPlot = FALSE)
tmp <- components(generalGraph, mode = "weak")
generalGraph <- delete.vertices(generalGraph, V(generalGraph)[!(tmp$membership == which.max(tmp$csize))])
rm(tmp)
generalGraph <- addGraphWeights(generalGraph)
generalSet <- graphToSet(generalGraph)

# Measure node properties
propDF <- generalSet$Nodes
propDF$KeyFunctionalDist <- apply(distances(generalGraph, to = V(generalGraph)[V(generalGraph)$InKeyFunctional], mode = "all"), 1, min, na.rm = TRUE)
propDF$Betweenness <- log(betweenness(generalGraph, normalized = FALSE)+1)
propDF$Betweenness <- propDF$Betweenness - min(propDF$Betweenness, na.rm = TRUE)
propDF$Betweenness <- propDF$Betweenness / mean(propDF$Betweenness, na.rm = TRUE)
propDF$Closeness <- closeness(generalGraph, mode = "all", normalized = FALSE)
propDF$Closeness <- propDF$Closeness - min(propDF$Closeness, na.rm = TRUE)
propDF$Closeness <- propDF$Closeness / mean(propDF$Closeness, na.rm = TRUE)
propDF$ReachCloseness <- closeness(generalGraph, mode = "out", normalized = FALSE)
propDF$ReachCloseness <- propDF$ReachCloseness - min(propDF$ReachCloseness, na.rm = TRUE)
propDF$ReachCloseness <- propDF$ReachCloseness / mean(propDF$ReachCloseness, na.rm = TRUE)
propDF$SinkCloseness <- closeness(generalGraph, mode = "in", normalized = FALSE)
propDF$SinkCloseness <- propDF$SinkCloseness - min(propDF$SinkCloseness, na.rm = TRUE)
propDF$SinkCloseness <- propDF$SinkCloseness / mean(propDF$SinkCloseness, na.rm = TRUE)

######################################################################
# Find and generate optimal subnetwork using key functional distance #
######################################################################

# Set node property to use as validation for network refinement (e.g. differential phosphorylation)
target <- propDF$PhosphoSig != 0

# Iterate over max distances from key functional nodes and record subnetwork properties
refinedRecordDF <- data.frame("MaxDist" = 0:30)
refinedRecordDF$size <- sapply(refinedRecordDF$MaxDist, function(x) sum(propDF$KeyFunctionalDist <= x))
refinedRecordDF$nodesCovered <- refinedRecordDF$size / nrow(propDF)
refinedRecordDF$numKinases <- sapply(refinedRecordDF$MaxDist, function(x) sum(propDF$Kinase & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist <= x)))
refinedRecordDF$kinasesCovered <- refinedRecordDF$numKinases / sum(propDF$Kinase & !propDF$InKeyFunctional)
refinedRecordDF$kinaseEnrich <- refinedRecordDF$kinasesCovered / refinedRecordDF$nodesCovered
refinedRecordDF$kinasePval <- sapply(refinedRecordDF$MaxDist, function(x) fisher.test(cbind(c(sum(propDF$Kinase & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist <= x)),
                                                                                              sum(!propDF$Kinase & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist <= x))),
                                                                                            c(sum(propDF$Kinase & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist > x)),
                                                                                              sum(!propDF$Kinase & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist > x)))))$p.value)
refinedRecordDF$kinaseLogPval <- -log10(refinedRecordDF$kinasePval)
refinedRecordDF$numTarget <- sapply(refinedRecordDF$MaxDist, function(x) sum(target & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist <= x)))
refinedRecordDF$targetCovered <- refinedRecordDF$numTarget / sum(target & !propDF$InKeyFunctional)
refinedRecordDF$targetEnrich <- refinedRecordDF$targetCovered / refinedRecordDF$nodesCovered
refinedRecordDF$targetPval <- sapply(refinedRecordDF$MaxDist, function(x) fisher.test(cbind(c(sum(target & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist <= x)),
                                                                                              sum(!target & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist <= x))),
                                                                                            c(sum(target & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist > x)),
                                                                                              sum(!target & !propDF$InKeyFunctional & (propDF$KeyFunctionalDist > x)))))$p.value)
refinedRecordDF$targetLogPval <- -log10(refinedRecordDF$targetPval)
rm(target)

# Set a distance threshold for max distance from key functional nodes, generate refined subnetwork, and calculate more properties
max_distance <- 11 # Set the functional distance cutoff for inclusion in the refined subnetwork

refinedSet <- pruneNetwork(generalSet, alwaysNodes = generalSet$Nodes$SYMBOL[generalSet$Nodes$InKeyFunctional], maxDistance = max_distance)
rm(max_distance)
refinedGraph <- plotCurrentNetwork(refinedSet, layout = layout_components, printPlot = FALSE)

refinedDF <- propDF[propDF$SYMBOL %in% refinedSet$Nodes$SYMBOL & propDF$Expressed!=-1,]
refinedDF$Betweenness <- refinedDF$Betweenness - min(refinedDF$Betweenness, na.rm = TRUE)
refinedDF$Betweenness <- refinedDF$Betweenness / mean(refinedDF$Betweenness, na.rm = TRUE)
refinedDF$Closeness <- refinedDF$Closeness - min(refinedDF$Closeness, na.rm = TRUE)
refinedDF$Closeness <- refinedDF$Closeness / mean(refinedDF$Closeness, na.rm = TRUE)
refinedDF$ReachCloseness <- refinedDF$ReachCloseness - min(refinedDF$ReachCloseness, na.rm = TRUE)
refinedDF$ReachCloseness <- refinedDF$ReachCloseness / mean(refinedDF$ReachCloseness, na.rm = TRUE)
refinedDF$SinkCloseness <- refinedDF$SinkCloseness - min(refinedDF$SinkCloseness, na.rm = TRUE)
refinedDF$SinkCloseness <- refinedDF$SinkCloseness / mean(refinedDF$SinkCloseness, na.rm = TRUE)
refinedDF$RefinedBetweenness <- log(betweenness(refinedGraph)+1)
refinedDF$RefinedBetweenness <- refinedDF$RefinedBetweenness - min(refinedDF$RefinedBetweenness, na.rm = TRUE)
refinedDF$RefinedBetweenness <- refinedDF$RefinedBetweenness / mean(refinedDF$RefinedBetweenness, na.rm = TRUE)
refinedDF$RefinedCloseness <- closeness(refinedGraph, mode = "all")
refinedDF$RefinedCloseness <- refinedDF$RefinedCloseness - min(refinedDF$RefinedCloseness, na.rm = TRUE)
refinedDF$RefinedCloseness <- refinedDF$RefinedCloseness / mean(refinedDF$RefinedCloseness, na.rm = TRUE)
refinedDF$RefinedReachCloseness <- closeness(refinedGraph, mode = "out")
refinedDF$RefinedReachCloseness <- refinedDF$RefinedReachCloseness - min(refinedDF$RefinedReachCloseness, na.rm = TRUE)
refinedDF$RefinedReachCloseness <- refinedDF$RefinedReachCloseness / mean(refinedDF$RefinedReachCloseness, na.rm = TRUE)
refinedDF$RefinedSinkCloseness <- closeness(refinedGraph, mode = "in")
refinedDF$RefinedSinkCloseness <- refinedDF$RefinedSinkCloseness - min(refinedDF$RefinedSinkCloseness, na.rm = TRUE)
refinedDF$RefinedSinkCloseness <- refinedDF$RefinedSinkCloseness / mean(refinedDF$RefinedSinkCloseness, na.rm = TRUE)

###################################################################
# Find and generate high-value subnetwork using refined closeness #
###################################################################

# Set the graph metric (e.g. closeness) to use as a cutoff and the target property (e.g. differential phosphorylation) to use as validation
metric <- refinedDF$RefinedCloseness
target <- refinedDF$PhosphoSig != 0

# Iterate over all values of the metric (i.e. refined closenes) and record subnetwork properties
highValRecordDF <- data.frame("Cutoff" = sort(unique(metric)))
highValRecordDF$size <- sapply(highValRecordDF$Cutoff, FUN = function(x) sum(metric > x))
highValRecordDF$nodeCovered <- highValRecordDF$size/nrow(refinedDF)
highValRecordDF$targetCovered <- sapply(highValRecordDF$Cutoff, FUN = function(x) sum(metric > x & target)/sum(target))
highValRecordDF$targetEnrich <- highValRecordDF$targetCovered / highValRecordDF$nodeCovered
highValRecordDF$targetPval <- sapply(highValRecordDF$Cutoff, FUN = function(x) fisher.test(cbind(c(sum(metric>x & target),sum(metric>x & !target)),
                                                                                                 c(sum(metric<=x & target),sum(metric<=x & !target))))$p.value)
highValRecordDF$logPval <- -log10(highValRecordDF$targetPval)
rm(metric,target)

# Enforce a refined closeness cutoff and generate high-value set and network
closeness_cutoff <- 1.1149274 # Set the refined closeness cutoff for inclusion in the high-value subnetwork

highValueDF <- refinedDF[refinedDF$RefinedCloseness > closeness_cutoff,]
rm(closeness_cutoff)
highValueSet <- pruneNetwork(refinedSet, alwaysNodes = highValueDF$SYMBOL, maxDistance = 0)
highValueGraph <- plotCurrentNetwork(highValueSet, layout = layout_components, printPlot = TRUE)

collapsedGraph <- collapseSiblings(highValueGraph)
collapsedSet <- graphToSet(collapsedGraph)
collapsedGraph <- plotCurrentNetwork(collapsedSet, layout = layout_components, printPlot = TRUE)

#####################################
# Writing output dataframes to CSVs #
#####################################

write.csv(refinedDF, file = "refinedDF_tmp.csv", row.names = FALSE)
write.csv(refinedRecordDF, file = "refined_record_tmp.csv", row.names = FALSE)

write.csv(highValueDF, file = "highValDF_tmp.csv", row.names = FALSE)
write.csv(highValRecordDF, file = "highVal_record_tmp.csv", row.names = FALSE)
