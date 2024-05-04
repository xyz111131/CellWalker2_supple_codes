args = commandArgs(trailingOnly=T)
case = args[1]
#info = read.csv(paste0("permute_Geshwind/cellwalk_integrate_all_3000_cor_wtree_",case, "_rand1_markers_combine_nEN_0.csv"), row.names=1) #rand, ss1, quan
#info = read.csv(paste0("permute_Geshwind/cellwalk_integrate_all_3000_cor_drop_all_nEN_wtree_",case, "_rand1_0.csv"), row.names=1) #rand, ss1, quan
#info = read.csv(paste0("permute_Geshwind/cellwalk_integrate_all_3000_cor_noise_3_wtree_",case, "_rand1_0.csv"), row.names=1) #rand, ss1, quan
info = read.csv(paste0("permute_Geshwind/cellwalk_integrate_all_3000_cor_",case, "_rand1_0.csv"), row.names=1) #rand, ss1, quan


info_rand =list()
info_mean = matrix(0, nrow(info), ncol(info))
info_std = matrix(0, nrow(info), ncol(info))
n = 60 #50
for(i in 1:n)
{
  #file = paste0("permute_Geshwind/cellwalk_integrate_all_3000_cor_wtree_", case, "_rand1_markers_combine_nEN_",i, ".csv")
  #file = paste0("permute_Geshwind/cellwalk_integrate_all_3000_cor_drop_all_nEN_wtree_", case, "_rand1_",i, ".csv")
  #file = paste0("permute_Geshwind/cellwalk_integrate_all_3000_cor_noise_3_wtree_", case, "_rand1_",i, ".csv")
  file = paste0("permute_Geshwind/cellwalk_integrate_all_3000_cor_", case, "_rand1_",i, ".csv")
  if(!file.exists(file)) { print(i); next }
  info_rand[[i]] = read.csv(file, row.names=1) #quan
  info_mean = info_mean + info_rand[[i]]
  info_std = info_std + info_rand[[i]]^2
}

info_mean = info_mean/n
info_std = info_std/n - info_mean^2 
zscore = (info - info_mean)/sqrt(info_std) 

#write.csv(zscore, file = paste0("cellwalk_integrate_all_3000_cor_drop_all_nEN_wtree_", case, "_Geshwind_rand1_zscore_1.csv")) # _1 only permute UCSC labels
write.csv(zscore, file = paste0("cellwalk_integrate_all_3000_cor_", case, "_Geshwind_rand1_zscore_1.csv")) # _1 only permute UCSC labels
write.csv(zscore, file = paste0("cellwalk_integrate_all_3000_cor_noise_3_wtree_", case, "_Geshwind_rand1_zscore_1.csv")) # _1 only permute UCSC labels
write.csv(zscore, file = paste0("cellwalk_integrate_all_3000_cor_wtree_", case, "_Geshwind_rand1_markers_combine_nEN_zscore_1.csv")) # _1 only permute UCSC labels


info_rand =list()
temp_mean = matrix(0, nrow(info), ncol(info))
info_std2 = matrix(0, nrow(info), ncol(info))
n = 20
for(i in 1:n)
{
  info_rand[[i]] = read.csv(paste0("bootstrap/cellwalk_integrate_all_3000_cor_wtree0.8_", case, "_ss0.8_",i, ".csv"), row.names=1)
  temp_mean = temp_mean + info_rand[[i]]
  info_std2 = info_std2 + info_rand[[i]]^2
}

temp_mean = temp_mean/n
info_std2 = info_std2/n - temp_mean^2 
zscore = (info - info_mean)/sqrt(info_std2)

print('info:')
print(info[1:5, 1:5])
print('bootstrap mean:')
print(temp_mean[1:5, 1:5])
print('null mean')
print(info_mean[1:5,1:5])
print('permute se')
print(sqrt(info_std[1:5,1:5]))
print('bs se')
print(sqrt(info_std2[1:5,1:5]))

write.csv(zscore, file = paste0("cellwalk_integrate_all_3000_cor_wtree0.8_",case,"_rand_zscore_ss0.8.csv")) #_2 edge weight between labels same as cell-label, ct_2 permutation and optimal weights using cutoff too
