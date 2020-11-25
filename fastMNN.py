#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:12:43 2020

@author: derrickliu
"""

import pandas as pd
import sklearn # scikit learn package will be very helpful
import numpy as np

from sklearn.decomposition import PCA  # for PCA
import matplotlib.pyplot as plt  # for plotting variance thresholds to pick the number of principal components to do the PCA with
%matplotlib inline

from sklearn.neighbors import NearestNeighbors  # for finding nearest neighbors between two data sets

# reading in the first pancreas data set

pancreas_data = pd.read_table("GSE81076_D2_3_7_10_17.txt")  
pancreas_data

# setting the gene names as the row names

pancreas_data = pancreas_data.rename(columns = {"Unnamed: 0":"Genes"}) 
pancreas_data.set_index(["Genes"], inplace = True)
pancreas_data

# function adapted from HW1 - cosine normalization 

def normalizing_cells_MaxAbsScaler(data):
    
    """
        input data: data frame with gene expression data  
            columns are the cells and rows are genes
        output data_normalized_data_values: normalized (scaled) data
        function:  dividing the value of each gene for each cell 
                by the maximum value of that cell.
    """
    import sklearn.preprocessing
    
    # creating the list of patients
    cells = data.columns.values
    
    #using Sklearn to scale the data 
    scaler = sklearn.preprocessing.MaxAbsScaler()
    data_scaled = scaler.fit_transform(data)
    
    #creating the dataframe, the output of sklearn MaxAbsScaler is an array
    data_normalized_data_values = pd.DataFrame(data_scaled, columns= cells, index = data.index)
      
    return data_normalized_data_values

# cosine normalization here keeps the gene and cell names in the rows and columns

normalized_pancreas = normalizing_cells_MaxAbsScaler(pancreas_data)
normalized_pancreas

# dropping missing values - looks like there aren't any here

normalized_pancreas.dropna(axis=1) 
normalized_pancreas

# graph to determine how many principal components should be selected for PCA of first data set

#Fitting the PCA algorithm with our Data
pca = PCA().fit(normalized_pancreas)
#Plotting the Cumulative Summation of the Explained Variance
plt.figure()
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('Number of Components')
plt.xticks([0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800])
plt.ylabel('Variance (%)') #for each component
plt.title('DR Dataset Explained Variance')
plt.show()

# Plot shows that we could do PCA with 100 principal components (could explain 95% of variance)
# It says the number of PCs should be = min(n_samples, n_features) in the class notes though.. --> I asked the prof about this

# MAY NOT NEED TO DO THIS BECAUSE FUNCTION IN NEXT CELL APPEARS TO DO SCALING
# PCA is effected by scale so you need to scale the features in your data before applying PCA. 
# ^ is from https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60 

from sklearn.preprocessing import StandardScaler
normalized_pancreas_stand = StandardScaler().fit_transform(normalized_pancreas)

# adapted from HW1 - trying to figure out number of principal components for PCA

from sklearn.feature_selection import VarianceThreshold

def variance_threshold_selector(data, threshold=0.04): #removes features with a variation below a cutoff
                                                       #we can play around with the threshold 
    selector = VarianceThreshold(threshold)
    selector.fit_transform(data)
    return data[data.columns[selector.get_support(indices=True)]]

topgenes_pancreas_1 = variance_threshold_selector(normalized_pancreas, 0.0005) #contains the features with a variation above the cutoff
topgenes_pancreas_1

# reading in the second pancreas data set - gene names are already the row names here

pancreas_data_2 = pd.read_table("GSE85241_cellsystems_dataset_4donors_updated.csv")  
pancreas_data_2

# cosine normalization for the second pancreas data set - gets rid of row names and column names though...

transformer = Normalizer().fit(pancreas_data_2)
pancreas_data_2_cnorm = pd.DataFrame(transformer.transform(pancreas_data_2))
pancreas_data_2_cnorm

# cosine normalization here keeps the gene and cell names in the rows and columns

normalized_pancreas_2 = normalizing_cells_MaxAbsScaler(pancreas_data_2)
normalized_pancreas_2

# dropping missing values - looks like there aren't any here

normalized_pancreas_2.dropna(axis=1) 
normalized_pancreas_2

topgenes_pancreas_2 = variance_threshold_selector(normalized_pancreas_2, 0.0005) #contains the features with a variation above the cutoff
topgenes_pancreas_2

# graph to determine how many principal components should be selected for PCA of second data set

#Fitting the PCA algorithm with our Data
pca = PCA().fit(normalized_pancreas_2)
#Plotting the Cumulative Summation of the Explained Variance
plt.figure()
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('Number of Components')
plt.xticks([0,200,400,600,800,1000,1200,1400,1600,1800])
plt.ylabel('Variance (%)') #for each component
plt.title('DR Dataset Explained Variance')
plt.show()

# Plot shows that we could do PCA with 100 principal components (could explain 98% of variance)
# It says the number of PCs should be = min(n_samples, n_features) in the class notes though..

