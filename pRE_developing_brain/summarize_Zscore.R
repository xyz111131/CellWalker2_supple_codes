case = 'info1'
region = 'TF_motif_signac' #'diseases_nearest_distance'#'diseases_region_temporal' #'layer_pfc'#'region_bg_cortex' #'layer_pfc' #'TF_motif_signac' #'TF_motif'#'DDD_region' #'region_unique' #cellregion_region
distan = 'log_Cosine'
#for(distan in  c('Cosine', 'Jaccard', 'Lsi'))
#{zeroD1
num_atac = 1500
# DevBrainCortex_integrate_cellwalk_tree2_all_TF_motif_signac_info1_label_rand_log_Jaccard_100.csv
# DevBrainCortex_integrate_cellwalk_tree2_all_DDD_region_info1_rand_log_Jaccard_51.csv
#DevBrainCortex_multi_cellwalk_TF_motif_signac_signac_info1_single_rand_log_Jaccard_100.csv
#DevBrainCortex_integrate_cellwalk_tree2_noATAC_TF_motif_signac_info_rand_log_Cosine_0.csv
for(seed in 3) #1:10
{
  #info = read.csv(paste0('results_compare_ATAC/DevBrainCortex_integrate_cellwalk_tree2_zeroD_', region, '_', case, '_rand_', distan,'_', num_atac, '_', seed, '_0.csv'), row.names=1) #metrics, linkATAC
   info = read.csv(paste0('results_compare_regions/DevBrainCortex_integrate_cellwalk_tree2_noATAC_', region, '_', case, '_rand_', distan,'_0.csv'), row.names=1) #label_rand
  #info = read.csv(paste0('results_compare_regions2/DevBrainCortex_integrate_cellwalk_tree_all_filter2_0.5_', region, '_', case, '_rand_', distan,'_0.csv'), row.names=1) #metrics, linkATAC, single_rand science
#info = read.csv(paste0('results/DevBrainCortex_multi_cellwalk_', region, '_', case, '_single_rand_', distan,'_0.csv'), row.names=1) #metrics, linkATAC
  #info = read.csv(paste0('results/DevBrainCortex_both_cellwalk_tree2_0.8_ASD_region_', case, '_ss0.8_0.csv'), row.names=1)
  #info = read.csv(paste0('results/DevBrainCortex_both_cellwalk_tree2_TF_motif_ds_', case, '_rand_0.csv'), row.names=1)
  

  info_rand =list()
  info_mean = matrix(0, nrow(info), ncol(info))
  info_std = matrix(0, nrow(info), ncol(info))
  #ni = 0
  n = 100
  for(i in c(1:n)) #11:20 is wrong
  {
    #info_rand[[i]]  = read.csv(paste0("results_compare_ATAC/DevBrainCortex_integrate_cellwalk_tree2_zeroD_", region, "_", 
  #				    case, "_rand_",distan, '_', num_atac, '_', seed, '_',i, '.csv'), row.names=1) #both
    #info_rand[[i]]  = read.csv(paste0("results/DevBrainCortex_both_cellwalk_tree2_TF_motif_ds_", case, "_rand_",i, '.csv'), row.names=1) #both
    info_rand[[i]] = read.csv(paste0('results_compare_regions/DevBrainCortex_integrate_cellwalk_tree2_noATAC_', region, '_', case, '_rand_', distan,'_',i,'.csv'), row.names=1) #metrics, linkATAC
    i#info_rand[[i]] = read.csv(paste0('results_compare_regions2/DevBrainCortex_integrate_cellwalk_tree_all_filter2_0.5_', region, '_', case, '_rand_', distan,'_',i,'.csv'), row.names=1) #metrics, linkATAC
    #info_rand[[i]] = read.csv(paste0('results/DevBrainCortex_multi_cellwalk_', region, '_', case, '_single_rand_', distan,'_',i,'.csv'), row.names=1) #metrics, linkATAC
    stopifnot(all(rownames(info) %in% rownames(info_rand[[i]])) )
    stopifnot(all(rownames(info_rand[[i]]) %in% rownames(info)) )
    info_mean = info_mean + info_rand[[i]][rownames(info), ]
    info_std = info_std + info_rand[[i]][rownames(info), ]^2
  #  ni = ni +1

  }
  
  info_mean = info_mean/n
  info_std = info_std/n - info_mean^2 
  zscore = (info - info_mean)/sqrt(info_std) 
  
  #write.csv(zscore, file = paste0("results/DevBrainCortex_multi_cellwalk_", region, "_", 
 # 				case,"_single_zscore_", distan, ".csv")) #both
  write.csv(zscore, file = paste0("results_compare_regions2/DevBrainCortex_integrate_cellwalk_tree_noATAC_filter2_0.5_", region, "_", 
  				case,"_zscore_", distan, ".csv")) #both,label
  
  write.csv(zscore, file = paste0("results_compare_regions/DevBrainCortex_integrate_cellwalk_tree2_noATAC_", region, "_", 
  				case,"_zscore_", distan, ".csv")) #both,label
  
  write.csv(info_mean, file = paste0("results_compare_ATAC/DevBrainCortex_integrate_cellwalk_tree2_zeroD_", region, "_", 
  				case,"_mean_", distan, "_", num_atac, '_', seed, ".csv")) #both
  
  write.csv(sqrt(info_std), file = paste0("results_compare_ATAC/DevBrainCortex_integrate_cellwalk_tree2_zeroD_", region, "_", 
  				case,"_std_", distan, "_", num_atac, '_', seed, ".csv")) #both
  
  write.csv(zscore, file = paste0("results_compare_ATAC/DevBrainCortex_integrate_cellwalk_tree2_zeroD_", region, "_", 
  				case,"_zscore_", distan, "_", num_atac, '_', seed, ".csv")) #both
}
#write.csv(zscore, file = paste0("DevBrainCortex_both_cellwalk_tree2_TF_motif_ds_",case,"_zscore.csv")) #both



info_rand =list() 
temp_mean = matrix(0, nrow(info), ncol(info))
#info_std2 = matrix(0, nrow(info), ncol(info))
for(i in 1:20)
{
  info_rand[[i]]  = read.csv(paste0("results/DevBrainCortex_both_cellwalk_tree2_0.8_ASD_region_", case, "_ss0.8_",i, '.csv'), row.names=1) #quan
  #temp_mean = temp_mean + info_rand[[i]]
  if(i == 1)
  {
     temp_mean = reshape2::melt(as.matrix(info_rand[[i]]))
  }else{
     temp_mean = cbind(temp_mean, reshape2::melt(as.matrix(info_rand[[i]]))[,3])
  }
  #info_std2 = info_std2 + info_rand[[i]]^2
}

#temp_mean = temp_mean/20
##info_std2 = info_std2/20 - temp_mean^2 
info_std2 = cbind(temp_mean[, 1:2], apply(temp_mean[, -1:-2], 1, function(x) mad(x)))
temp_mean = cbind(temp_mean[, 1:2], apply(temp_mean[, -1:-2], 1, function(x) median(x)))
temp_mean = reshape2::acast(temp_mean, Var1~Var2)
info_std2 = reshape2::acast(info_std2, Var1~Var2)
temp_mean = temp_mean[rownames(info), colnames(info)]
info_std2 = info_std2[rownames(info), colnames(info)]
zscore = (info - info_mean)/info_std2  #sqrt(info_std2)

print('info:')
print(info[1:5, 1:5])
print('bootstrap mean:')
print(temp_mean[1:5, 1:5])
print('null mean')
print(info_mean[1:5,1:5])
print('null std')
print(sqrt(info_std[1:5,1:5]))
print('bootstrap std')
print(info_std2[1:5,1:5])

write.csv(zscore, file = paste0("DevBrainCortex_both_cellwalk_tree2_0.8_ASD_region_",case,"_zscore_ss0.8.csv")) 
