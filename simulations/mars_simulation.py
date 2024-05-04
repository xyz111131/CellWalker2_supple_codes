#!/usr/bin/env python
# coding: utf-8

# In[1]:


import torch
import sys
import numpy as np
import pandas as pd
import scanpy.api as sc
from anndata import AnnData
from scipy.optimize import linear_sum_assignment 
#import h5py
#import anndata as ad

#get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib import pyplot as plt
import matplotlib as mpl

sys.path.append("../")
from args_parser import get_parser
from model.mars import MARS
from model.experiment_dataset import ExperimentDataset
from data.benchmarks import BenchmarkData
import warnings
warnings.filterwarnings('ignore')


# # Setting parameters

# Loading default parameters

# In[2]:


params, unknown = get_parser().parse_known_args()
prefix = "/pollard/home/zhhu/cellwalk/simulation_new/"

# In[3]:



# Checking if CUDA device is available

# In[4]:


if torch.cuda.is_available() and not params.cuda:
        print("WARNING: You have a CUDA device, so you should probably run with --cuda")
device = 'cuda:0' if torch.cuda.is_available() and params.cuda else 'cpu'
params.device = device


# In[5]:


def init_seed(opt):
    torch.cuda.cudnn_enabled = False
    np.random.seed(opt.manual_seed)
    torch.manual_seed(opt.manual_seed)
    torch.cuda.manual_seed(opt.manual_seed)
init_seed(params)



# # Prepare merged data

# In[6]:


dat = pd.read_csv(prefix + "data/simulation" + sys.argv[1] + "_integrated_scaled.csv", index_col = 0) #3000 × 7498


# In[7]:


adata = AnnData(dat.values.transpose())
adata.obs_names = dat.columns
adata.var_names = dat.index



# In[8]:


meta = pd.read_csv(prefix + "data/simulation" + sys.argv[1] + "_integrated_metadata.csv", index_col = 0) 



# In[10]:


adata.obs["ground_truth"] = meta.Group # Categoricals are preferred for efficiency
adata.obs["experiment"] = meta.Batch



# In[ ]:


#sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
#sc.pp.scale(adata, max_value=10, zero_center=True)


# In[11]:


# sc.pp.neighbors(adata, n_neighbors=30, use_rep='X')
# sc.pp.pca(adata, n_comps=50)
# sc.tl.tsne(adata, n_jobs = 5 )
# sc.pl.tsne(adata, color=['experiment','ground_truth'],size=50)



# In[13]:

#datasets = list(set(adata.obs['experiment']))
#datasets


# # Train and evaluate MARS

# ### Use science dataset as annotated, and geschwind as unannotated 

# Prepare annotated, unannotated and pretrain datasets

# In[14]:


science = adata[adata.obs['experiment'] == 'Batch1',:]
geschwind = adata[adata.obs['experiment'] == 'Batch2',:]


# In[15]:


def celltype_to_numeric(adata, annotation_type):
       """Adds ground truth clusters data."""
       annotations = list(adata.obs[annotation_type])
       annotations_set = sorted(set(annotations))
       
       mapping = {a:idx for idx,a in enumerate(annotations_set)}
       
       truth_labels = [mapping[a] for a in annotations]
       adata.obs['truth_labels'] = pd.Categorical(values=truth_labels)
        
       return mapping


# In[ ]:


# from sklearn.preprocessing import LabelEncoder
# le = LabelEncoder()
# label = le.fit_transform(kolod.obs['ground_truth'])
# label2 = le.fit_transform(pollen.obs['ground_truth'])
# 


# In[16]:


celltype_id_map_science = celltype_to_numeric(science, 'ground_truth')
IDs_to_celltypes_science = {v:k for k,v in celltype_id_map_science.items()}
celltype_id_map_geschwind = celltype_to_numeric(geschwind, 'ground_truth')
IDs_to_celltypes_geschwind = {v:k for k,v in celltype_id_map_geschwind.items()}


# In[17]:


y_science = np.array(science.obs['truth_labels'], dtype=np.int64)
annotated = ExperimentDataset(science.X, science.obs_names, science.var_names, 'science', y_science)


# In[18]:


y_geschwind = np.array(geschwind.obs['truth_labels'], dtype=np.int64) # ground truth annotations will be only used for evaluation
unannotated = ExperimentDataset(geschwind.X, geschwind.obs_names, geschwind.var_names, 'geschwind', y_geschwind)


# In[19]:


pretrain_data = ExperimentDataset(geschwind.X, geschwind.obs_names, geschwind.var_names, 'geschwind')


# In[20]:


n_clusters = len(np.unique(unannnotated.y))




# Initialize MARS

# In[39]:


params.epochs = 100
params.epochs_pretrain = 50
params.learning_rate = 1e-3
mars = MARS(n_clusters, params, [annotated], unannotated, pretrain_data, hid_dim_1=500, hid_dim_2=100)


# Run MARS in evaluation mode. Ground truth annotations will be used to evaluate MARS performance and scores will be returned

# In[40]:


# return both annotated and unannotated datasets with save_all_embeddings
adata, landmarks, scores = mars.train(evaluation_mode=True, save_all_embeddings=True) # evaluation mode


# In[25]:


print(scores)


# Check MARS performance

# In[26]:


def hungarian_match(y_true, y_pred):
    """Matches predicted labels to original using hungarian algorithm."""
    
    #y_true = adjust_range(y_true)
    #y_pred = adjust_range(y_pred)
    
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)
    # Confusion matrix.
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    ind = linear_sum_assignment(-w)
    ind = np.asarray(ind)
    ind = np.transpose(ind)
    d = {i:j for i, j in ind}
    y_pred = np.array([d[v] for v in y_pred])
    
    return y_true, y_pred, d



# In[29]:


geschwind = adata[adata.obs['experiment'] == 'geschwind',:]
y_true, y_pred, d = hungarian_match(geschwind.obs['truth_labels'], geschwind.obs['MARS_labels'])


# In[31]:


name_maps = mars.name_cell_types(adata, landmarks, IDs_to_celltypes_science)



# In[33]:


confusion = pd.DataFrame(data=np.zeros((4, 5)), columns = celltype_id_map_science.keys()) #4 5
for idx in name_maps:
    for x in name_maps[idx]:
        confusion.loc[idx,x[0]] = x[1]

confusion.index = [IDs_to_celltypes_geschwind[d[v]] for v in d.keys()]


# In[ ]:


confusion.to_csv(prefix + 'results/' +  sys.argv[1] + '_confusion.csv')


# In[34]:


print(confusion)


# Visualize in MARS embedding space



# In[ ]:

# #create anndata object using MARS embeddings as X
# adata_mars = AnnData(adata.obsm['MARS_embedding'])
# adata_mars.obs['MARS_labels'] = pd.Categorical(adata.obs['MARS_labels'])
# adata_mars.obs['ground_truth'] = pd.Categorical(adata.obs['truth_labels'])
# adata_mars.obs['experiment'] = pd.Categorical(adata.obs['experiment'])


# In[ ]:


# # visualize only unannotated dataset
# pollen_mars = adata_mars[adata_mars.obs['experiment'] == 'geschwind',:]
# sc.pp.neighbors(pollen_mars, n_neighbors=30, use_rep='X')
# sc.tl.umap(pollen_mars)
# sc.pl.umap(pollen_mars, color=['ground_truth','MARS_labels'],size=50)


# # In[ ]:


# #visualize annotated and unannotated datasets jointly
# sc.pp.neighbors(adata_mars, n_neighbors=30, use_rep='X')
# sc.tl.umap(adata_mars)
# sc.pl.umap(adata_mars, color=['experiment','ground_truth'], size=50)



