# -*- coding: utf-8 -*-
"""
Created on Sat Aug  9 13:05:32 2025


Plot for the PCA analysis 

@author: Zach

"""
from re import search
from scipy.optimize import curve_fit
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pybaselines import Baseline
import pandas as pd
import numpy as np

def SNIP_Baseline_Correction(absorbances, wavenumbers):
    baseline_fitter = Baseline(x_data=wavenumbers)
    window_size = 10
    snip_baseline, snip_params = baseline_fitter.snip(absorbances, window_size, extrapolate_window=20)
    return absorbances - snip_baseline

#define a useful function for applying PCA to IR data, from Grossutti et al. 2022
def IR_PCA(input_data, PCs = 2):
    from sklearn.decomposition import PCA
    import pickle as pk
    
    pca = PCA(n_components = PCs)
    data_PCs = pca.fit_transform(input_data)
    columns_list = ['PC'+ str(x) for x in range(1, PCs+1)]
    PC_df = pd.DataFrame(data = data_PCs,
                                columns = columns_list)
    
    feature_weights = pca.components_
    loadings_df = pd.DataFrame(feature_weights)
    
    result = pca.fit_transform(input_data) # Assume X is having more than 2 dimensions    
    
    print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))
    
    return(PC_df, loadings_df)

def prep_dataframe(filepath):
    
    df = pd.read_csv(filepath, skiprows=[1,2])
    try:
        df['aged_time'] = df['Source ID'].str.strip().str.extract(r'(\d+)').astype(int)
        labels = df.loc[:,'Source ID'].values.copy()
        labels_number = np.array([int(search(r'\d+', label).group()) for label in labels])
    except:
        df['aged_time'] = df['Sample Name'].str.strip().str.extract(r'(\d+)').astype(int)
        labels = df.loc[:,'Sample Name'].values.copy()
        labels_number = np.array([int(search(r'\d+', label).group()) for label in labels])

    # Sort by the extracted numeric value
    sorted_df = df.sort_values('aged_time').drop(columns='aged_time')
    
    return sorted_df

def sort_df(df):
    try:
        df['aged_time'] = df['Source ID'].str.strip().str.extract(r'(\d+)').astype(int)
        labels = df.loc[:,'Source ID'].values.copy()
        labels_number = np.array([int(search(r'\d+', label).group()) for label in labels])
    except:
        df['aged_time'] = df['Sample Name'].str.strip().str.extract(r'(\d+)').astype(int)
        labels = df.loc[:,'Sample Name'].values.copy()
        labels_number = np.array([int(search(r'\d+', label).group()) for label in labels])
    
    # Sort by the extracted numeric value
    sorted_df = df.sort_values('aged_time').drop(columns='aged_time')
    
    return sorted_df

path = "./"
path += "Purelink+ 0-4000hrs TE-MCT ATR PCA Data.csv"

df = prep_dataframe(path)

wavenumbers = df.columns.values[3:].astype(float)
wavenumbers.sort()

spectra = df[wavenumbers.astype(str)]

spectra['Sample Name'] = df['Sample Name']

cmap = plt.get_cmap('viridis')
colors = cmap(np.linspace(0, 0.8, len(spectra['Sample Name'].unique())))

fig, ax = plt.subplots(figsize=(10, 5))
ax.set_ylabel('Mean Spectra', fontsize=32)
#ax.set(yticklabels=[])
ax.set_xlabel(r'Wavenumber (cm$^-$$^1$)', fontsize=28)
ax.set_xlim(wavenumbers[-1], wavenumbers[0])

ax.tick_params(which='major', axis='x', labelsize=24)
ax.tick_params(which='minor', axis='x', length=10)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
#ax.yaxis.set_major_locator(plt.NullLocator())

i=0
for label, group in spectra.groupby('Sample Name'):
    baselined_absorbances = SNIP_Baseline_Correction(group.iloc[:,:-1].mean(), wavenumbers)
    ax.plot(wavenumbers, baselined_absorbances.values, label=label, color=colors[i])
    i+=1

ax.tick_params(axis='y', labelsize=24)
ax.legend(fontsize=24, loc='upper left', bbox_to_anchor=(0.1, 0.52, 0.5, 0.5)) 
ax.grid(True, which='major', linestyle='-', linewidth=1, color='grey', alpha=0.8)
plt.grid(True, which='minor', linestyle='--', linewidth=0.5, color='grey', alpha=0.5)
ax.minorticks_on()
#plt.show()


"""
Below is the Script to plot the principle component analysis method of the spectra above
"""

PC_df , loading_df = IR_PCA(spectra[wavenumbers.astype(str)])

PC_df['Sample Name'] = df['Sample Name']

PC_df = sort_df(PC_df)

grouped_PC_dfs = PC_df.groupby(['Sample Name'])

PC_dfs = {key: grouped_PC_dfs.get_group(key) for key in grouped_PC_dfs.groups.keys()}

legend_labels = PC_df['Sample Name'].unique()

PC_dfs = [PC_dfs[key] for key in legend_labels]

cmap = plt.get_cmap('viridis')
colors = cmap(np.linspace(0, 0.8, len(legend_labels)))

fig, ax = plt.subplots(figsize=(10, 5))
ax.set_ylabel('PC1', fontsize=28)
#ax.set(yticklabels=[])
ax.set_xlabel('PC2', fontsize=28)
for i, label in enumerate(legend_labels):
    ax.scatter(PC_dfs[i]['PC2'], PC_dfs[i]['PC1'], label=label, color=colors[i])

ax.tick_params(which='major', axis='x', labelsize=24)
ax.tick_params(which='minor', axis='x', length=10)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
#ax.yaxis.set_major_locator(plt.NullLocator())
ax.tick_params(axis='y', labelsize=24)
ax.legend(fontsize=24, loc='best') 
ax.grid(True, which='major', linestyle='-', linewidth=1, color='grey', alpha=0.8)
plt.grid(True, which='minor', linestyle='--', linewidth=0.5, color='grey', alpha=0.5)
ax.minorticks_on()
plt.show()

fig, ax = plt.subplots(figsize=(10, 5))
ax.set_ylabel('PC1', fontsize=28)
#ax.set(yticklabels=[])
ax.set_xlabel('PC2', fontsize=28)
for i, label in enumerate(legend_labels):
    ax.scatter(PC_dfs[i]['Sample Name'], PC_dfs[i]['PC1'], label=label, color=colors[i])

ax.tick_params(which='major', axis='x', labelsize=24)
ax.tick_params(which='minor', axis='x', length=10)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
#ax.yaxis.set_major_locator(plt.NullLocator())
ax.tick_params(axis='y', labelsize=24)
ax.legend(fontsize=24, loc='best') 
ax.grid(True, which='major', linestyle='-', linewidth=1, color='grey', alpha=0.8)
plt.grid(True, which='minor', linestyle='--', linewidth=0.5, color='grey', alpha=0.5)
ax.minorticks_on()
plt.show()