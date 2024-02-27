# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 15:18:38 2023

@author: aosinuga2
"""

#%%
import os
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from adjustText import adjust_text
# conda install -c conda-forge adjusttext 
# pip install adjustText

#%%
sample = 'C' ;
#sample = 'D' ;

condition='1E6';  # WITH OUTLIERS
#condition='NoOutlier'   # WITHOUT OUTLIERS (NO OUTLIERS)
#condition='Outliers'    # ONLY OUTLIERS

# List of indices to remove
indices_to_remove = [39]
tmp_filename = '_File_'+condition+'.csv' 
dpi_save = 500

#%%
#%%
# Create the folder if it does not exist
folder_name = '../results/CodeFigures/'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
   
#%%
    
# Read heatmap_data from CSV
heatmap_data = pd.read_csv('../results/'+sample+'_FCdata_'+tmp_filename, header=None).values.T

enz_labels = pd.read_csv('../results/'+sample+'_Enzymelabel_'+tmp_filename, header=None).values

rxn_labels = pd.read_csv('../results/'+sample+'_Rxnlabel_'+tmp_filename, header=None).values

#%%
heatmap_data = np.delete(heatmap_data, indices_to_remove, axis=0)
enz_labels = np.delete(enz_labels, indices_to_remove, axis=0)
rxn_labels = np.delete(rxn_labels, indices_to_remove, axis=0)
#%%
n_samples = heatmap_data.shape[0]
n_features = heatmap_data.shape[1]
numPCAComponents = min(n_samples, n_features)
# Checking PCA before t-SNE
pca = PCA(n_components=numPCAComponents)
# heatmap_data_pca = pca.fit_transform(heatmap_data_std)
heatmap_data_pca = pca.fit_transform(heatmap_data)
# Plot the explained variance
explained_variance = pca.explained_variance_ratio_ * 100
plt.bar(range(1, numPCAComponents+1), explained_variance)
plt.xlabel('Principal Component')
plt.ylabel('Explained Variance (%)')
plt.title('Scree Plot of Explained Variance')
#plt.show()
plt.savefig(folder_name+sample+'_'+condition+'_Scree_plot.png', format='png', dpi = dpi_save) # Save as high

#%%
perplexity = np.arange(5, n_samples-1, 1)
divergence = []
for i in perplexity:
    model = TSNE(n_components=2, init="pca", perplexity=i)
    reduced = model.fit_transform(heatmap_data)
    divergence.append(model.kl_divergence_)
plt.plot(perplexity, divergence, marker='o', color="red", linewidth=1)
plt.xlabel("Perplexity Values")
plt.ylabel("Divergence")
plt.title("Perplexity vs Divergence")
plt.savefig(folder_name+sample+'_'+condition+'_Perplex_Divergence_plot.png', format='png', dpi = dpi_save) # Save as high
#plt.show()

#%%
perplexity = n_samples-1
randomState = 42 

#%
print("Step II_e: t-SNE for enzyme classes")
# 2-D t-SNE embedding
tsne = TSNE(n_components=2, perplexity=perplexity, method='barnes_hut', init="pca", random_state=randomState)
Ytsne = tsne.fit_transform(heatmap_data)
dftsne = pd.DataFrame(Ytsne)
dftsne['enzID'] = enz_labels
dftsne['rxnID'] = rxn_labels
dftsne.columns = ['y1','y2','enzID', 'rxnID']

# Scatter plot using Seaborn
plt.figure(num=None, figsize=None, dpi=800, facecolor=None, edgecolor=None, frameon=True)
sns.scatterplot(data=dftsne, x='y1', y='y2', hue='enzID', palette='tab10')
plt.xlabel('Component 1', fontweight='bold')
plt.ylabel('Component 2', fontweight='bold')
plt.title('Fold Changes Distribution with t-SNE: 2-D Embedding', fontweight='bold')
# # Annotate each point with the cluster label
# for index, row in dftsne.iterrows():
#     plt.annotate(row['enzID'], (row['y1'], row['y2']))

# Annotate each point with the cluster label using adjustText
texts = []
for index, row in dftsne.iterrows():
    texts.append(plt.text(row['y1'], row['y2'], str(row['rxnID']), fontsize=10))

# Adjust the position of labels to avoid overlap
adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
plt.legend(title='Enzyme', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(folder_name+'Fig_6_tsne2dPlot'+sample+'_'+condition+'.png', format='png', dpi = dpi_save) # Save as high-res PNG
# plt.show()


#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score

# Range of values for the hyperparameters
eps_values = np.arange(0.2, 0.001, -0.001)
min_samples_values = range(10, 40)

# Store the best score and corresponding hyperparameters
best_score = -1
best_eps = 0
best_min_samples = 0

for eps in eps_values:
    for min_samples in min_samples_values:
        # Perform DBSCAN
        db = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1)
        labels = db.fit_predict(Ytsne)
        
        # Calculate the silhouette score
        if len(np.unique(labels)) > 1:  # Silhouette score is only meaningful if there are more than one cluster.
            score = silhouette_score(Ytsne, labels)
            print(f'EPS: {eps}, min_samples: {min_samples}, Score: {score}')
            
            # Check if this is the best score so far for 3 clusters... 
            if (score > best_score) and max(labels)==2:
                best_score = score
                best_eps = eps
                best_min_samples = min_samples

# Output the best settings
print(f'Best EPS: {best_eps}')
print(f'Best min_samples: {best_min_samples}')
print(f'Best silhouette score: {best_score}')

# if can't fit into 3 clusters reduce the min_samples: 
if best_eps==0:
  # Range of values for the hyperparameters
  eps_values = np.arange(0.2, 0.001, -0.001)
  min_samples_values = range(1, 40)
  
  # Store the best score and corresponding hyperparameters
  best_score = -1
  best_eps = 0
  best_min_samples = 0
  
  for eps in eps_values:
      for min_samples in min_samples_values:
          # Perform DBSCAN
          db = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1)
          labels = db.fit_predict(Ytsne)
          
          # Calculate the silhouette score
          if len(np.unique(labels)) > 1:  # Silhouette score is only meaningful if there are more than one cluster.
              score = silhouette_score(Ytsne, labels)
              print(f'EPS: {eps}, min_samples: {min_samples}, Score: {score}')
              
              # Check if this is the best score so far for 3 clusters... 
              if (score > best_score) and max(labels)==2:
                  best_score = score
                  best_eps = eps
                  best_min_samples = min_samples
  
  # Output the best settings
  print(f'Best EPS: {best_eps}')
  print(f'Best min_samples: {best_min_samples}')
  print(f'Best silhouette score: {best_score}')

#%%

db = DBSCAN(eps=best_eps, min_samples=best_min_samples, n_jobs=-1)
dftsne['cluster'] = db.fit_predict(Ytsne)
print(dftsne)

import openpyxl  # Make sure to import openpyxl

file_name = 'DBSCANclusteringResults.xlsx'
full_path = '../results/' + file_name
sheet_name = sample + '_n_1E6_' + condition  # Ensure sample and condition are defined

# Use ExcelWriter to append or replace sheet in the Excel file
with pd.ExcelWriter(full_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
    # Load existing book
    if os.path.exists(full_path):
        writer.book = openpyxl.load_workbook(full_path)
        # Check if the sheet exists and remove it before appending/replacing
        if sheet_name in writer.book.sheetnames:
            # This step is now optional given 'if_sheet_exists' set to 'replace'
            std = writer.book[sheet_name]
            writer.book.remove(std)
    # Write DataFrame to an Excel sheet
    dftsne.to_excel(writer, index=False, sheet_name=sheet_name)

print("DataFrame is appended to Excel File on specified Sheet successfully.")
