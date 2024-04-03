import os 
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
from IPython.display import Image
import scanpy as sc
from cnmf import cNMF


for i in ["CB_1","CB_2"]:
  cnmf_obj = cNMF(output_dir="./03_res1",name = i)
  cnmf_obj.prepare(counts_fn = i+'.count.txt',components = 10,n_iter=300,seed=14,num_highvar_genes=2000)
  cnmf_obj.factorize()
  cnmf_obj.combine()
  cnmf_obj.consensus(k=10,density_threshold=0.1,show_clustering=True)
