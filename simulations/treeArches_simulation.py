#!/usr/bin/env python
# coding: utf-8

# # treeArches: learning and updating a cell-type hierarchy (basic tutorial)
# 
# In this tutorial, we explain the different functionalities of treeArches. We show how to:
# 
# - [Step 1](#Create-scVI-model-and-train-it-on-reference-dataset): Integrate reference datasets using scVI
# - [Step 2](#Construct-hierarchy-for-the-reference-using-scHPL): Match the cell-types in the reference datasets to learn the cell-type hierarchy of the reference datasets using scHPL
# - [Step 3](#Use-pretrained-reference-model-and-apply-surgery-with-a-new-query-dataset-to-get-a-bigger-reference-atlas): Apply architural surgery to extend the reference dataset using scArches
# - [Step 4a](#Updating-the-hierarchy-using-scHPL): Update the learned hierarchy with the cell-types from the query dataset using scHPL (useful when the query dataset is labeled)
# - [Step 4b](#Predicting-cell-type-labels-using-scHPL): Predict the labels of the cells in the query dataset using scHPL (useful when the query dataset is unlabeled)

# In[2]:

import sys
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)


# In[3]:


os.environ["OMP_NUM_THREADS"] = "8" # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] = "8" # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] = "8" # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] = "8" # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] = "8" # export NUMEXPR_NUM_THREADS=1


# In[ ]:



# In[4]:


import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown
import copy as cp
import seaborn as sns


# In[5]:


sc.settings.set_figure_params(dpi=1000, frameon=False)
sc.set_figure_params(dpi=1000)
sc.set_figure_params(figsize=(7,7))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42


# ### Convert raw Dataset to h5ad

# In[6]:


import pandas as pd
from anndata import AnnData


prefix = "/pollard/home/zhhu/cellwalk/simulation_new/"


# 4-1: cell type 5 is ancestor of 3 and 4; 4-0: cell type 5 is a new cell type; 4: cell type 5 is slightly different from cell type 4

# In[7]:


#dat = pd.read_csv('/pollard/home/zhhu/cellwalk/simulation_new/simulation4_merge_rawcount.csv', index_col=0)
#meta = pd.read_csv('/pollard/home/zhhu/cellwalk/simulation_new/simulation4_merge_metadata.csv', index_col=0)


dat = pd.read_csv(prefix + "data2/simulation" + sys.argv[1] + "_merge_rawcount.csv", index_col = 0) 
meta = pd.read_csv(prefix + "data2/simulation" + sys.argv[1] + "_merge_metadata.csv", index_col = 0)



# In[8]:


adata = AnnData(dat.values.transpose())
adata.obs_names = dat.columns
adata.var_names = dat.index
adata.obs['study'] =  meta['Batch']
adata.obs['ground_truth'] = pd.Categorical(meta['Group'])


# We now split the data into reference and query dataset to simulate the building process. Here we use the 'geschwind' batch as query data.

# In[ ]:


# target_conditions = ["geschwind"]
# source_adata = adata[~adata.obs.study.isin(target_conditions)].copy()
# target_adata = adata[adata.obs.study.isin(target_conditions)].copy()
# print(source_adata)
# print(target_adata)


# In[9]:


source_adata = adata


# For a better model performance it is necessary to select HVGs. We are doing this by applying the function `scanpy.pp.highly_variable_genes()`. The parameter `n_top_genes` is set to 2000 here. However, for more complicated datasets you might have to increase number of genes to capture more diversity in the data.

# In[10]:


source_adata.raw = source_adata




# In[11]:


sc.pp.normalize_total(source_adata)
sc.pp.log1p(source_adata)


# In[12]:


sc.pp.highly_variable_genes(
    source_adata,
    n_top_genes=1000,
    batch_key= "study",
    subset=True)


# For consistency we set adata.X to be raw counts. In other datasets that may be already the case

# In[13]:


source_adata.X = source_adata.raw[:, source_adata.var_names].X #1000 genes


# ### Create scVI model and train it on reference dataset

# 
# <div class="alert alert-warning">
# <b>Remember:</b> The adata object has to have count data in adata.X for scVI/scANVI if not further specified.
# 
# </div>
# 
# 

# In[14]:


sca.models.SCVI.setup_anndata(source_adata, batch_key="study") #None


# The scVI model uses the zero-inflated negative binomial (ZINB) loss by default. Insert `gene_likelihood='nb'` to change the reconstruction loss to negative binomial (NB) loss.

# In[15]:


vae = sca.models.SCVI(
    source_adata,
    n_layers=1,
    n_latent = 20,
    encode_covariates=True,
    deeply_inject_covariates=True,
    use_layer_norm="both",
    use_batch_norm="none",
)


# In[16]:


vae.train(max_epochs=100, early_stopping = True)




# In[17]:


#vae.get_reconstruction_error(indices = vae.validation_indices)
#vae.get_reconstruction_error(indices = vae.train_indices)


# In[18]:


#plt.plot(vae.history['train_loss_epoch'])


# The resulting latent representation of the data can then be visualized with UMAP

# In[19]:


reference_latent = sc.AnnData(vae.get_latent_representation())
reference_latent.obs["cell_type"] = source_adata.obs["ground_truth"].tolist()
#reference_latent.obs["batch"] = source_adata.obs["batch"].tolist()
reference_latent.obs["study"] = source_adata.obs["study"].tolist()


# In[ ]:


#reference_latent


# In[20]:


# sc.pp.neighbors(reference_latent, n_neighbors=8)
# sc.tl.leiden(reference_latent)
# sc.tl.umap(reference_latent)


# In[21]:


reference_latent.obs['study'] = reference_latent.obs['study'].astype('category')

# Reorder categories, so smallest dataset is plotted on top
#reference_latent.obs['study'].cat.reorder_categories(['Oetjen', 'Sun', 'Freytag'], inplace=True)


# In[22]:


# sc.pl.umap(reference_latent,
#            color=['cell_type'], #cell_type
#            frameon=False,
#            wspace=0.6, s=25, 
#            palette=sns.color_palette('pastel', as_cmap=True)[:5]
#            )




# ### Construct hierarchy for the reference using scHPL
# 
# First, we concatenate all cell type labels with the study labels. This way, we ensure that the cell types of the different studies are seen as unique.

# <div class="alert alert-warning">
# <b>Warning:</b> Always ensure that the cell type labels of each dataset are unique!
# 
# </div>

# In[23]:


reference_latent.obs['celltype_batch'] = np.char.add(np.char.add(np.array(reference_latent.obs['cell_type'], dtype= str), '-'),
                                             np.array(reference_latent.obs['study'], dtype=str))


# Now, we are ready to learn the cell-type hierarchy. In this example we use the `classifier='knn'`, this can be changed to either a linear SVM (`'svm'`) or a one-class SVM (`'svm_occ'`). We recommend to use the kNN classifier when the dimensionality is low since the cell-types are not linearly separable anymore. 
# 
# The option `dynamic_neighbors=True` implies that the number of neighbors changes depending on the number of cells in the dataset. If a cell-type is small, the number of neighbors used will also be lower. The number of neighbors can also be set manually using `n_neighbors`.
# 
# During each step of scHPL, a classifier is trained on the datasets we want to match and the labels are cross-predicted. If you're interested in the confusion matrices used for the matching, set `print_conf=True`. The confusion matrices are also saved to .csv files then.
# 
# For more details about other parameters, take a look at the scHPL [GitHub](https://github.com/lcmmichielsen/scHPL)

# In[28]:


# import sys
# orig_stdout = sys.stdout
# f = open('out.txt', 'w')
# sys.stdout = f

tree_ref, mp_ref = sca.classifiers.scHPL.learn_tree(data = reference_latent, 
                batch_key = 'study',
                batch_order = ['Batch1', 'Batch2'],
                cell_type_key='celltype_batch',
                classifier = 'knn', dynamic_neighbors=True,
                dimred = False, print_conf= True) # will output NC1.csv and NC2.csv

# sys.stdout = orig_stdout
# f.close()


# In[ ]:


# tree_ref, mp_ref = sca.classifiers.scHPL.learn_tree(data = reference_latent, 
#                 cell_type_key='celltype_batch',
#                  batch_key = 'celltype_batch',
#                 batch_order = ['Group1-Batch1','Group2-Batch1','Group3-Batch1','Group4-Batch1',
#                                'Group1-Batch2','Group2-Batch2','Group3-Batch2','Group5-Batch2'],
#                 classifier = 'knn', dynamic_neighbors=True,
#                 dimred = False, print_conf= True)




