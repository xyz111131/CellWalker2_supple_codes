# tf occurence
pRE_motif <- read.csv('GSE162170/pRE_motif_signac.csv')
motifs = unique(motif$motif_alt_id)
cooccur = NULL
for(i in 1:(length(motifs)-1))
{
  aa = pRE_motif$sequence_name[pRE_motif$motif_alt_id == motifs[i]]
  print(motifs[i])
  scores = sapply((i+1):length(motifs), function(j){
    bb = pRE_motif$sequence_name[pRE_motif$motif_alt_id == motifs[j]]
    length(intersect(aa, bb))/ length(union(aa, bb))
  })
  cooccur = rbind(cooccur, data.frame("motif1" = rep(motifs[i], length(scores)), "motif2" = motifs[-1:-i], "scores" = scores))
}

write.csv(cooccur, file = 'pRE_motif_signac_cooccur_Jaccard.csv', row.names = F)

ggplot(cooccur, aes(x = motif1, y = motif2, fill = scores)) + geom_tile()
