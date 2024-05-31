Tutorial 2: Slide-seq dataset
====================================

.. raw:: html

    <div style="font-size: 15px;">The source code package is freely available at https://github.com/cuiyaxuan/DiffusionST. The datasets used in this study can be found at https://drive.google.com/drive/folders/1qgn2UKpu4q14ysCoCKjWYVEHXIzHNoqq?usp=drive_link.</div>

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    <div style="font-size: 15px;"> First, cd /home/.../DiffusionST-main </div>

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    from DenoiseST import DenoiseST
    import os
    import torch
    import pandas as pd
    import numpy as np
    import scanpy as sc
    from sklearn import metrics
    import multiprocessing as mp

.. code:: ipython3

    device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
    
    n_clusters = 7
    file_fold = '/home/cuiyaxuan/spatialLIBD/151673'
    adata = sc.read_visium(file_fold, count_file='151673_filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()

.. code:: ipython3

    model = DenoiseST(adata,device=device,n_top_genes=4096)
    adata = model.train()

.. code:: ipython3

    from repair_model import main_repair
    df=pd.DataFrame(adata.obsm['emb'])
    main_repair(adata,df,device)

.. code:: ipython3

    csv_file = "example.csv"
    data_df = pd.read_csv(csv_file, header=None)
    data_df = data_df.values
    adata.obsm['emb'] = data_df

.. code:: ipython3

    from utils import clustering
    
    radius = 50
    tool = 'mclust' # mclust, leiden, and louvain
    if tool == 'mclust':
       clustering(adata, n_clusters, radius=radius, method=tool, refinement=True)
    elif tool in ['leiden', 'louvain']:
       clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)
    
    df=adata.obs['domain']
    df.to_csv("label_HP.csv")

.. code:: ipython3

    import matplotlib as mpl
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import warnings
    import visual_high
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams["font.sans-serif"] = "Arial"
    warnings.filterwarnings('ignore')
    file_fold = '/home/cuiyaxuan/spatialLIBD/6.Mouse_Hippocampus_Tissue/'
    #adata = sc.read_h5ad(file_fold + 'mouse_anterior_posterior_brain_merged.h5ad')
    adata = sc.read_h5ad(file_fold + 'filtered_feature_bc_matrix_200115_08.h5ad')
    df_label=pd.read_csv('/home/cuiyaxuan/DiffusionST/Diffusion_Hip_10/label_HP.csv', index_col=0) 
    #df_label=pd.read_csv('./label_5000.csv', index_col=0) ##If the dropout rate is less than 0.85, visualize the data using "label_5000.csv".
    visual_high.visual(adata,df_label)


.. parsed-literal::

    #cells after MT filter: 53208
    WARNING: saving figure to file figures/spatialHippocampus.pdf



.. image:: test3_files/test3_6_1.png
   :width: 315px
   :height: 277px






