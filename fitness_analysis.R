library(phytools); library(ape); library(Biostrings); library(seqinr);
library(lubridate); library(magrittr); library(dplyr); library(stringr); library(RColorBrewer)
library(ggplot2); library(data.table); library(nlme)
source("D:/Documents/ARPE/Koelle lab-archive/Koelle lab/Phylodynamic analysis/R_codes/functions.R")
pattern='12-19'

setwd("~/ARPE/Koelle lab-archive/Koelle lab/Phylodynamic analysis/world_12-19")

### Run treetime and load the tree ---------------------------------------------
system(paste('/anaconda3/bin/treetime --aln nt', prefix, '_subsamp.fasta --tree ', "tree",prefix, '.treefile --clock-rate 0.004 --dates dates_tt.csv --outdir treetime004', sep=''), wait=TRUE)
(tre.tt <- ape::read.nexus("./treetime_large/timetree.nexus"))

### Load and process metadata --------------------------------------------------
tre.tt <- drop.tip(tre.tt, c('EPI_ISL_278569','EPI_ISL_223588','EPI_ISL_225316','EPI_ISL_278719','EPI_ISL_270745','EPI_ISL_242425','EPI_ISL_278549','EPI_ISL_278570','EPI_ISL_242562','EPI_ISL_242560','EPI_ISL_242484','EPI_ISL_242498','EPI_ISL_205025'))
meta = read.csv(paste('data_',pattern,'_subsamp.csv', sep=''))
meta <- meta[meta$Isolate_Id %in% tre.tt$tip.label, ]
meta_nodes <- data.frame(Isolate_Id = tre.tt$node.label)
tip.order <- c(1:length(tre.tt$tip.label))
names(tip.order) <- meta$Isolate_Id
tip.order <- tip.order[tre.tt$tip.label]
meta_tree <- bind_rows(meta[tip.order, ], meta_nodes) #build a new metadata folder with all computed informations from all nodes, in the order of the absolute node numbers
meta_tree$Decimal_Date <- time_nodes(tre.tt, meta)
row.names(meta_tree) <- meta_tree$Isolate_Id
meta_tree$Isolate_Id <- NULL

rm(meta); rm(tip.order);rm(meta_nodes)

### compute the fitness on each node -------------------------------------------
Alignment <- readDNAStringSet (paste("./treetime_large/ancestral_sequences.fasta", sep=""))
AA_Prefs <- read.csv(file= 'AA_prefs_avg.csv', header = TRUE)
Alignment <- Alignment[c(tre.tt$tip.label, tre.tt$node.label)] #reorder in the order of the absolute node numbers
Alignment_AA <- Biostrings::translate(Alignment)
list_mut_AA = apply(tre.tt$edge, 1, function(edge) mismatches(Alignment_AA[ edge[2] ], Alignment_AA[ edge[1] ]))
rm(Alignment)

list_ME = lapply(list_mut_AA, function(ListMut) sum_ME(ListMut))
meta_tree$Fitness <- fitness_nodes(tre.tt) #compute the fitness of each node/tip, with 0 for the root

### identify the antigenic clusters --------------------------------------------
meta_tree <- define_ag_clusters(meta_tree, Alignment_AA)

# plot the antigenic clusters
library(network)
edge.col= meta_tree[tre.tt$edge[,2],'VaxStrain']
plot.phylo(tre.tt, show.tip.label = FALSE, edge.color = as.color(edge.col))

# dst_small<- tre.tt_small
# dst_small$edge.length[dst_small$edge.length==0]<-max(nodeHeights(tre.tt_small))*1e-6 # set zero-length branches to be 1/1000000 total tree length
# nsim=25

### infer the phylogeography ---------------------------------------------------
dst<- tre.tt
dst$edge.length[dst$edge.length==0]<-max(nodeHeights(tre.tt))*1e-6 # set zero-length branches to be 1/1000000 total tree length

# compute the transition matrix
region <- meta_tree$Region[1:length(tre.tt$tip.label)]
names(region) <- tre.tt$tip.label
fit <- fitMk(dst, region, model='ARD', pi='estimated')
fittedQ<-matrix(NA,length(fit$states),length(fit$states))
fittedQ[] <- c(0,fit$rates)[fit$index.matrix+1]
diag(fittedQ) <- 0
diag(fittedQ) <- -rowSums(fittedQ)
colnames(fittedQ) <- rownames(fittedQ) <- fit$states
write.csv(fittedQ,file='fittedQ_large.csv')

#if importing transition matrix from csv file
fittedQ.tmp <- read.csv(file='fittedQ_small.csv')
fittedQ <- data.matrix(fittedQ.tmp[,2:9])
rownames(fittedQ)<-colnames(fittedQ)
region <- meta_tree[1:length(tre.tt$tip.label),'Region']
names(region) <- rownames(meta_tree[1:length(tre.tt$tip.label),])

# run the stochastic mapping
library(parallel)
trees.sm<-mclapply(1:5,function(n,tree,x,fixedQ) make.simmap(tree,x,Q=fixedQ,nsim=3, pi=fit$pi),
                   tree=dst,x=region,fixedQ=fittedQ,mc.cores=if(.Platform$OS.type=="windows") 1L else 5L)

trees.sm.mem <- trees.sm
tmp <- c()
for(t in trees.sm) tmp <- c(tmp, t)
trees.sm <- tmp
rm(tmp)
class(trees.sm)<-c("multiSimmap","multiPhylo")
write.simmap(trees.sm, 'trees.sm15')

#if importing stochastic mappings from file
trees.sm <- read.simmap(file='trees.sm100', format='phylip', version=1.0)

#plot the stochastic mappings
pdf(file='stochastic_mappings_large.pdf')
png(file='stochastic_mappings_large6.png', width=1200, height=1000)
cols=setNames(c('firebrick2','forestgreen','deepskyblue3', 'chocolate1','darkmagenta','chocolate4'),
          c('Europe', 'North_America', 'Oceania', 'China', 'SE_Asia', 'Southern_Asia'))

cols=setNames(brewer.pal(6,'Set1'),sort(unique(meta_tree$Region)))
# for (t in 1:5){
  plotSimmap(trees.sm[[10]], ftype='off', lwd=2, colors=cols)
  add.simmap.legend(prompt=FALSE, x=0.5, y=0.8*par()$usr[4], colors=cols)
  axisPhylo(lwd=2, yaxp=c(par()$usr[4],5,1))
# }
dev.off()

# merge Europe and North America into a super-region
trees.sm.hem <- lapply(trees.sm, function(t) mergeMappedStates(t,old.states = c('North_America','Europe'), new.state = 'N_Hemisphere' ))

### Analyze each stochastic mapping --------------------------------------------
#Identify the migration events of each stochastic realization

trunk_tree <- trunk(tre.tt)
length(trunk_tree)
edges.col <-rep(1,nrow(tre.tt$edge))
names(edges.col)<-tre.tt$edge[,2]
edges.col[trunk_tree] <- 2
plot.phylo(tre.tt, show.tip.label = FALSE, edge.color = as.color(edges.col))

Migr_per_tree <- mclapply(trees.sm, function(t){
  get_members_succes(t, trunk.tree = trunk_tree, meta = meta_tree)
}, mc.cores=if(.Platform$OS.type=="windows") 1L else 5L)

summary_proba_obs('Europe', simmaps = trees.sm, meta = meta_tree, Migr_per_tree, trnk = trunk_tree, ncores = 6L)
summary_proba_obs('Oceania', simmaps = trees.sm, meta = meta_tree, Migr_per_tree, trnk = trunk_tree, ncores= 6L)
summary_proba_obs('North_America', simmaps = trees.sm, meta = meta_tree, Migr_per_tree, trnk = trunk_tree, ncores = 6L)

summary_slopes(simmaps = trees.sm[c(1,3)], meta = meta_tree, ncores = 1L, trnk = trunk_tree)
