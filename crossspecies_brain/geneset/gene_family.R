# gene family
gene_fam <- readRDS('~/Dropbox (Gladstone)/single_cell/cross-species/Analysis_lein_bdbag_2020_08_11-2023-02-21_11.48.28/data/Transcriptomics/sncell/10X/human/processed/analysis/analysis/M1/cross_species_integration/MetaNeighbor/hgnc_syngo.rds')
library(ggplot2)
# plot z-score
aa = read.csv("results/geneset/cellwalk_integrate_all_wtree_gf_info1_rand_zscore.csv", row.names=1) # 0
aa = t(aa)

library(GO.db)
# extract a named vector of all terms
goterms <- Term(grep('^GO:', colnames(aa), value = T))
colnames(aa)[grep('^GO:', colnames(aa))] = goterms

aa[is.na(aa)] = 0
#aa[aa < 50] = NA #qnorm(0.05/nrow(aa)/ncol(aa))

sel = which(colSums(aa>50)  >0)
aa = aa[, sel]
aa2 = t(aa)
aa2[aa2<50] = 0
ord = hclust(dist(aa2, method = 'manhattan')) # clustering similar TF motifs, zscore1
aa = aa[, ord$order]


tr1 = readRDS('human_Inh_phylo.rds') # order cell types by trees
tr2 = readRDS('marmoset_Inh_phylo.rds')
tr3 = readRDS('mouse_Inh_phylo.rds')

ns1 = 2*tr1$Nnode + 1
ns2 = 2*tr2$Nnode + 1
ns3 = 2*tr3$Nnode + 1



aa1 = aa[1:ns1, ] 
aa2 = aa[(ns1 + 1) :(ns1 + ns2), ] 
aa3 = aa[-1 :-(ns1 + ns2), ] 

df = read.csv('results/cellwalk_integrate_all_wtree_info_rand_zscore_max_mapping2human.csv', row.names = 1)# mapping between species
rownames(df) = gsub('[ |-]', '.', rownames(df))

# optional
aa1[aa1 < quantile(c(aa1), 0.95)] = NA
aa2[aa2 < quantile(c(aa2), 0.95)] = NA
aa3[aa3 < quantile(c(aa3), 0.95)] = NA


aa1 = aa1[tr1$edge[,2],]
aa1 = aa1[1:36, 1:24] #part1
id2 = which(colnames(aa1) == 'Calcium voltage-gated channel subunits')
id1 = which(rownames(aa1) == 'Inh.L3.5.SST.CDH3.Inh.L5.6.SST.KLHL1.7')
id3 = which(rownames(aa1) == 'Inh.L6.SST.TH') 
aa1 = aa1[id3:id1, 25:id2] #part2
id2 = which(colnames(aa1) == 'Zinc fingers')
id4 = which(colnames(aa1) == 'regulation of postsynaptic density assembly')-1
id3 = which(rownames(aa1) == 'Inh.L6.SST.TH')+1
id5 = which(colnames(aa1) == 'Phosphatases')
id6 = which(colnames(aa1) == 'postsynaptic actin cytoskeleton')
aa1 = aa1[id3:nrow(aa1), c(id2:id4, id5:id6)] #part3

id3 = which(rownames(aa1) == 'Inh.L3.5.VIP.TAC3')
id2 = which(colnames(aa1) == 'regulation of postsynaptic density assembly')
id4 = which(colnames(aa1) == 'Neurexins')+1
id5 = which(colnames(aa1) == 'ZF class homeoboxes and pseudogenes')
id6 = ncol(aa1)
aa1 = aa1[37:id3, c(id2:id4, id5:id6)] #part4

## for not matching cell types
aa2 = aa2[tr2$edge[,2], ]
aa2 = aa2[3:31, 1:24] #part1

aa3 = aa3[tr3$edge[,2], ]
aa3 = aa3[2:29, 1:24] #part1

## separate plots per species
aa1 = reshape2::melt(as.matrix(aa1))
colnames(aa1) = c("cell_types", "gene_set", "zscore") 
###aa1$zscore[aa1$zscore < quantile(aa1$zscore, 0.95)] = NA

aa2 = reshape2::melt(as.matrix(aa2))
colnames(aa2) = c("cell_types", "gene_set", "zscore") 

aa3 = reshape2::melt(as.matrix(aa3))
colnames(aa3) = c("cell_types", "gene_set", "zscore") 


#aa = aa[c(tr1$edge[,2], tr2$edge[, 2] + ns1, tr3$edge[, 2] + ns1 + ns2), ]
# aa = as.matrix(aa)
# aa = reshape2::melt(aa)
# colnames(aa) = c("gene_set", "cell types", "zscore") 

pdf("results/geneset/cellwalk_integrate_all_wtree_gf_info1_rand_zscore_mouse.pdf",width= 45, height = 22 ) #width= 20, height = 15 
# ggplot(aa, aes(x= human, y= marmoset, fill=zscore, group= human)) + 
#   geom_tile() + 
#   theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   scale_fill_gradient(low="white", high="blue") + xlab('') + ylab('') #, limit= c(0,0.01)
print(ggplot(aa1, aes(x= gene_set, y= cell_types, group=gene_set, size = zscore, color = zscore)) +
        geom_point(alpha = 0.8) + #xlab('') +
        theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1)) + #45, 30
        scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab")) #, limit = c(0, 1)) #+scale_size(range = c(0.1, 8)))

print(ggplot(aa2, aes(x= gene_set, y= cell_types, group=gene_set, size = zscore, color = zscore)) +
        geom_point(alpha = 0.8) + #xlab('') +
        theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1)) + #45, 30
        scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab")) #, limit = c(0, 1)) #+scale_size(range = c(0.1, 8)))

print(ggplot(aa3, aes(x= gene_set, y= cell_types, group=gene_set, size = zscore, color = zscore)) +
        geom_point(alpha = 0.8) + #xlab('') +
        theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1)) + #45, 30
        scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab")) #, limit = c(0, 1)) #+scale_size(range = c(0.1, 8)))
dev.off()

## combined plots & all zscore
all(colnames(aa1) == colnames(aa2)); all(colnames(aa3) == colnames(aa2))
aa1 = aa1[rownames(aa1) %in% gsub('[ |-]', '.', tr1$tip.label), ]# only tips
aa2 = aa2[df[rownames(aa1), 'marmoset'], c(id2:id4, id5:id6)] #part1: 1:24, part2: 25:id2, part3:c(id2:id4, id5:id6)
rownames(aa2) = make.names(rownames(aa2), unique = T)
aa3 = aa3[df[rownames(aa1), 'mouse'], c(id2:id4, id5:id6)]
rownames(aa3) = make.names(rownames(aa3), unique = T)

## for not matching cell types
aa2 = aa2[rownames(aa2) %in% gsub('[ |-]', '.', tr2$tip.label), ]# only tips
aa3 = aa3[rownames(aa3) %in% gsub('[ |-]', '.', tr3$tip.label), ]# only tips

aa1_m = reshape2::melt(as.matrix(aa1))
aa2_m = reshape2::melt(as.matrix(aa2))
aa3_m = reshape2::melt(as.matrix(aa3))
to_plot = rbind(aa1_m, aa2_m, aa3_m)
colnames(to_plot) = c("cell_types", "gene_set", "zscore") 
to_plot$species = c(rep('human', nrow(aa1_m)), rep('marmoset', nrow(aa2_m)), rep('mouse', nrow(aa3_m)))
to_plot$zscore[to_plot$zscore < 0] = 0
library(viridis)
gg<-ggplot(to_plot, aes(x= cell_types, y= gene_set, group=cell_types, fill = zscore)) +
        geom_tile(alpha = 0.8) + #xlab('') +
        theme(axis.text.x = element_text(angle = 90, hjust=1), strip.background = element_blank()) + scale_fill_viridis() + 
        #scale_fill_gradient(low = "mediumblue",  high = "red2", space = "Lab")  + 
        facet_wrap(~species, ncol=3, strip.position = 'top', scales = 'free_x')
gg <- gg + theme(axis.ticks=element_blank())
#gg <- gg + theme(axis.text=element_text(size=5))
gg <- gg + theme(panel.border=element_blank())
gg <- gg + theme(plot.title=element_text(hjust=0))
gg <- gg + theme(strip.text=element_text(hjust=0))
gg <- gg + theme(panel.margin.x=unit(0.5, "cm"))
gg <- gg + theme(panel.margin.y=unit(0.5, "cm"))
#gg <- gg + theme(legend.title=element_text(size=6))
gg <- gg + theme(legend.title.align=1)
#gg <- gg + theme(legend.text=element_text(size=6))
gg <- gg + theme(legend.position="bottom")
gg <- gg + theme(legend.key.size=unit(0.4, "cm"))
gg <- gg + theme(legend.key.width=unit(1, "cm"))

pdf('results/geneset/cellwalk_integrate_stack_wtree_gf_info1_rand_zscore_part4_map_celltype.pdf',width= 11, height = 7)
plot(gg)
dev.off()
# interactive plot
library("plotly")
fig <- plot_ly(
  x = colnames(aa1), y = rownames(aa1),
  z = aa1, type = "heatmap"
)

fig
