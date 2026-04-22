# Code for the case studies of "SOUPE: an R package to test and optimize the phylogenetic and ecological representativeness of a taxonomic sample."
# Nyniane Steinkampf--Pellecuer, Idriss Pelletan, Pauline Provini

#--------------
# Libraries
#--------------

library(SOUPE)
library(ape)
library(RColorBrewer)
library(phytools)


#------------------------------------------------------------------------------------------------------------------------
# CS1
# A posteriori process: evaluation of a sample used in a study on birds’ locomotion abilities (Provini & Höfling 2020)
#------------------------------------------------------------------------------------------------------------------------


# files are available in the SM2 folder


# --------------------------
# Phylogenetic relevance
# --------------------------

tested.sample <- read.table("SM2a.txt", header=FALSE, sep=";")
tested.sample <- tested.sample$V1

tree <- read.nexus("SM2c.nex")
tree

output1 <- compost(tree = tree, sample = tested.sample, list = tree$tip.label, n = 10000)
output1.1 <- result.compost(tree = tree, data = output1)


# ------------------------
# Ecological relevance
# ------------------------

traits.table.CS1 <- read.csv("SM2b.csv", header = TRUE, sep=";")
traits.table.CS1

output2 <- metrics(tree = tree, sample = tested.sample, traits_table = traits.table.CS1, faith = FALSE, percentage = FALSE, sub_phylo = FALSE)
output2$Traits


#---------------------------------------------------------------------------------------------------------------------
# CS2
# A priori process: creation of a sample for the digitization campaign of birds’ fluid collection in natural history
# collections (e_COL+ project).
#---------------------------------------------------------------------------------------------------------------------


# files are available in the SM3 folder


# ------------------------
# opti.traits function
# ------------------------
tree<-read.nexus("SM3d.nex")
tree

available.taxa <- read.table("SM3a.txt", header = FALSE, sep = ",")
available.taxa <- available.taxa$V1
available.taxa

keyspe.list <- read.table("SM3b/keypse.txt", header = FALSE, sep = ",")
keyspe.list <- keyspe.list$V1
keyspe.list

traits.table.CS2 <- read.csv("SM3c.csv", header = TRUE, sep=";")
traits.table.CS2

prop.ecolp <- list(Size = c("ok" = 1, "not_ok" = 0))
weight.ecolp <- c(1)

output3 <- opti.traits(tree = tree, n = 172, list = available.taxa, keyspe = keyspe.list, traits_table = traits.table.CS2, proportions = prop.ecolp, traits_weight = weight.ecolp, delta = 0.5, history = FALSE, sub_phylo = FALSE)

output3$Best_sample
output3$Faith
output3$Percentage


# check phylogenetic representativity
# -----------------------------------------

#to run compost, create the list of available species as the list of species in the collections under 40cm
species.good.size <- c()
for (i in 1:length(traits.table.CS2$Species)) {
  print(i)
  if (traits.table.CS2$Size[[i]] == c("ok")) {
    species.good.size <- c(species.good.size, traits.table.CS2$Species[[i]])
  }
}
species.good.size

output4 <- compost(tree = tree, sample = output3$Best_sample, list = species.good.size, n = 10000)
output4.1 <- result.compost(tree = tree, data = output4)


# visualization of the sample on the phylogeny
# -------------------------------------------------

sample1 <- output3$Best_sample

# Visualize the wanted, scanned and selected species on the phylogeny of available species
pruned.tree<-drop.tip(tree,tree$tip.label[-match(available.taxa, tree$tip.label)]) #subphylogeny of available species
plot(pruned.tree, type="fan")

n<-length(sample1)
n

pa<-data.frame((matrix(NA,ncol=2,nrow=length(pruned.tree$tip.label))))
colnames(pa) <- c("Species", "Presence")

wanted_spec <- read.table("SM3b/wanted_spec.txt", header = FALSE, sep = ",")
wanted_spec <- wanted_spec$V1
wanted_spec

scanned_spec <- read.table("SM3b/scanned_spec.txt", header = FALSE, sep = ",")
scanned_spec <- scanned_spec$V1
scanned_spec

presence <- pruned.tree$tip.label %in% sample1 #species in available taxa present in the created sample
presence

for (i in 1:length(pruned.tree$tip.label)){ #create a table with status for each of the available species
  pa$Species[[i]] <- pruned.tree$tip.label[i]
  if (presence[[i]] ==TRUE) {
    if((pruned.tree$tip.label[i] %in% scanned_spec)==TRUE) {
      pa$Presence[[i]] <-c("scanned")
    }
    else if((pruned.tree$tip.label[i] %in% wanted_spec)==TRUE) {
      pa$Presence[[i]] <-c("wanted")
    }
    else {pa$Presence[[i]] <- c("sample") }
  }
  else {
    pa$Presence[[i]] <- c("absent")
  }
}

trait <- pa$Presence
names(trait) <- pa$Species
trait
colors<-list(setNames(c("#da627d","#00ff00ff","#0f9ed5ff", "#f5ebe0"),c("scanned", "priority", "opti.traits_sample", "not_taken")))
colors

dev.off()
plotFanTree.wTraits(pruned.tree, trait, lwd=12, colors=colors, ftype="off")


# Visualize species size on the available species tree
trait2 <- traits.table.CS2$Size
names(trait2) <- traits.table.CS2$Species
trait2
colors2<-list(setNames(c("#fdd06bff","#ff0000ff"),c("ok","not_ok")))
colors2

dev.off()
plotFanTree.wTraits(pruned.tree, trait2, lwd=12, colors=colors2, ftype="off")




# ------------------------
# bestsample function
# ------------------------

tree<-read.nexus("SM3d.nex")
tree

available.taxa <- read.table("SM3a.txt", header = FALSE, sep = ",")
available.taxa <- available.taxa$V1
available.taxa

keyspe.list <- read.table("SM3b/keyspe.txt", header = FALSE, sep = ",")
keyspe.list <- keyspe.list$V1
keyspe.list

traits.table.CS2 <- read.csv("SM3c.csv", header = TRUE, sep=";")
traits.table.CS2

species.good.size <- c()
for (i in 1:length(traits.table.CS2$Species)) {
  print(i)
  if (traits.table.CS2$Size[[i]] == c("ok")) {
    species.good.size <- c(species.good.size, traits.table.CS2$Species[[i]])
  }
}
species.good.size

output5 <- bestsample(tree = tree, init = keyspe.list, type = "number", n = 172, list = species.good.size, sub_phylo = FALSE)
output5$Best_sample
output5$Faith
output5$PR_percentage


# Check phylogenetic representativity
# -----------------------------------------

output6 <- compost(tree = tree, sample = output5$best_sample, list = species.good.size, n = 10000)
output6.1 <-result.compost(tree = tree, data = output6)



# Visualization of the sample on the phylogeny
# -----------------------------------------------

sample2 <- output5$best_sample


# Visualize the wanted, scanned and selected species on the phylogeny of available species
pruned.tree<-drop.tip(tree,tree$tip.label[-match(available.taxa, tree$tip.label)]) #subphylogeny of available species
plot(pruned.tree, type="fan")

n<-length(sample2)
n

wanted_spec <- read.table("SM3b/wanted_spec.txt", header = FALSE, sep = ",")
wanted_spec <- wanted_spec$V1
wanted_spec

scanned_spec <- read.table("SM3b/scanned_spec.txt", header = FALSE, sep = ",")
scanned_spec <- scanned_spec$V1
scanned_spec

pa2<-data.frame((matrix(NA,ncol=2,nrow=length(pruned.tree$tip.label))))
colnames(pa2) <- c("Species", "Presence")
pa2

presence2 <- pruned.tree$tip.label %in% sample2 #species in available taxa present in the created sample
presence2

for (i in 1:length(pruned.tree$tip.label)){ #create a table with status for each of the available species
  pa2$Species[[i]] <- pruned.tree$tip.label[i]
  if (presence2[[i]] ==TRUE) {
    if((pruned.tree$tip.label[i] %in% scanned_spec)==TRUE) {
      pa2$Presence[[i]] <-c("scanned")
    }
    else if((pruned.tree$tip.label[i] %in% wanted_spec)==TRUE) {
      pa2$Presence[[i]] <-c("wanted")
    }
    else {pa2$Presence[[i]] <- c("sample") }
  }
  else {
    pa2$Presence[[i]] <- c("absent")
  }
}
pa2

trait3 <- pa2$Presence
names(trait3) <- pa2$Species
trait3
colors3<-list(setNames(c("#da627d","#00ff00ff","#0055ffff", "#f5ebe0"),c("scanned", "priority", "bestsample_sample", "not_taken")))
colors3

dev.off()
plotFanTree.wTraits(pruned.tree, trait3, lwd=12, colors=colors3, ftype="off")


# Visualize species size on the available species tree
trait4 <- traits.table.CS2$Size
names(trait2) <- traits.table.CS2$Species
trait4
colors4<-list(setNames(c("#fdd06bff","#ff0000ff"),c("ok","not_ok")))
colors4

dev.off()
plotFanTree.wTraits(pruned.tree, trait2, lwd=12, colors=colors2, ftype="off")



# -------------------
# -------------------
