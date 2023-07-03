# Installation ------------------------------------------------------------
# install.packages("devtools")
install.packages('Seurat')                                          #install Seurat

devtools::install_github("satijalab/seurat-data", ref = 'develop')  #install Seurat-data

if (!requireNamespace("remotes", quietly = TRUE)) {                #install Seurat-disk
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
install.packages("anndata")

devtools::install_github("sqjin/CellChat")                     #install cellchat

devtools::install_github("thomasp85/patchwork")                #install patchwork

setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")  #install signac

if (!require("BiocManager", quietly = TRUE))                  #install scater
  install.packages("BiocManager")

BiocManager::install("scater")

devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")    #install loomR


library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(anndata)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)
library(Signac)
library(loomR)
library(scater)

# Reading Data  -----------------------------------------------------------

#"Convert from Scanpy to Seurat...
adata <- read_h5ad("D:/Deep learning Single Cell Course 5-2023/LAB4 Cellchat/MelanomaAnnNoS.h5ad")
adata
dataseurat<- CreateSeuratObject(counts = t(adata$X), meta.data = adata$obs)
Melnoma <- dataseurat 
Melnoma

# Part I: Data input & processing and initialization of CellChat ----------

#1- Create a CellChat object
cellchat <- createCellChat(object = Melnoma , group.by = "cell_type") #, assay = "RNA", datatype = "spatial", scale.factors = scale.factors,  do.sparse = T
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

#2- Set the ligand-receptor interaction database(Human)
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Paracrine Autocrine Signaling interactions
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

#3-Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat) #with 20 warnings
cellchat <- identifyOverExpressedInteractions(cellchat) # with waranings 

# # Part II: Inference of cell-cell communication network -------- --------
#1-Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1) #with warnings
cellchat@options$parameter #Parameter_Values_Storage

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#2-Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)
View(df.net)

df.net1 <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
View(df.net1)

df.net2<- subsetCommunication(cellchat, signaling = c("VEGF", "BMP"))
View(df.net2)

#3-Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
View(cellchat@netP)
cellchat@netP$pathways
cellchat@netP[["pathways"]]

#4-Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat) #Counting Number of links 
cellchat

#We can also VISUALIZE the aggregated cell-cell communication network. 
#For example, showing the number of interactions or the total interaction strength (weights) 
#between any two cell groups using CIRCLE plot.
groupSize

par(mfrow = c(1,2), xpd=TRUE) # Display 1 row 2 column / xpd = true all plotting is clipped to the figure region,
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
# Just in case : Error in plot.new() : figure margins too large!!
#ERROR SOLVE: dev.off()
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")

#Due to the complicated cell-cell communication network, 
#we can examine the signaling sent from each cell group. 
#Here we also control the parameter edge.weight.max 
#so that we can compare edge weights between differet networks.
mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#Expected Error in plot.new() : figure margins too large , Solution : dev.off()


# # Part III: Visualization of cell-cell communication network  -----------
#Data Exploration, Analysis, and Visualization.
#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
cellchat@netP$pathways #All the signaling pathways showing signifi cantcommunications  
pathways.show <- c("VEGF")

# 1-Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to CAFs and the right portion shows signaling to other cells 
vertex.receiver = seq(1,2) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
# netVisual_individual(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")

#2- Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

#3-Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

#4-Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# *Chord diagram* grouping cell clusters into different cell types.
group.cellType <- c(rep("CAFs", 5), rep("Endo", 5), rep("Macro", 5), rep("Mel", 5), rep("T_cells/NK", 5)) # grouping cell clusters into CAFs, Endo and Macro cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

#Compute the contribution of each ligand-receptor pair to the overall signaling 
netAnalysis_contribution(cellchat, signaling = pathways.show) 

#visualize the cell-cell communication mediated by a single ligand-receptor pair. 
#We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.EGVF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.EGVF[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,2) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout = "hierarchy")
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

##Automatically save the plots of the all inferred network for quick exploration

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,2)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

##Visualize cell-cell communication mediated by MULTIPLE LIGAND-RECEPTORS or signaling pathways
#a-Bubble plot 
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use' Endo) to other cell groups (defined by 'targets.use' all)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathway (VEGF)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5), signaling = c("VEGF"), remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("PTN","VEGF", "APRIL"))
netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(2:5), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object
#b-Chord Digram
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use':ENDO) to other cell groups (defined by 'targets.use')
# show all the interactions sending from Endo
netVisual_chord_gene(cellchat, sources.use = 2, targets.use = c(1:5), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by CAFS from Endo
netVisual_chord_gene(cellchat, sources.use = 2, targets.use = 1, legend.pos.x = 15)  

#c-Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "VEGF")
plotGeneExpression(cellchat, signaling = "PTN")
plotGeneExpression(cellchat, signaling = "ANGPT")
plotGeneExpression(cellchat, signaling = "BMP")
plotGeneExpression(cellchat, signaling = "APRIL")

# # Part IV: Systems analysis of cell-cell communication network ----------

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 12, height = 5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space  
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("VEGF"))
gg2
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("VEGF"))
ht


#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
#Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
#To infer number of patterns
selectK(cellchat, pattern = "outgoing") # takes time!!
nPatterns = 2
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 2.
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

#Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming") #Takes time be patient:)
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function
#> # dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

#>Manifold and classification learning analysis of signaling networks
#>dentify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


# # Part V: Save the CellChat object- -------------------------------------

saveRDS(cellchat, file = "cellchat_Melnoma_LS.rds")
x = readRDS("cellchat_Melnoma_LS.rds")

x = readRDS("D:/CELLCHAT2/Dr Abdelrahman TASK/Editting/Melnoma.single.r")

# Subset specific columns and rows from Seurat object using subset function
subset_obj <- subset(Melanoma, select = y, cells =c())
