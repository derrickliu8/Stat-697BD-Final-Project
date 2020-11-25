#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 14:36:53 2020

@author: derrickliu
"""

import pandas as pd
import sklearn # scikit learn package will be very helpful
import numpy as np

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

# reading in the second pancreas data set - gene names are already the row names here

pancreas_data_2 = pd.read_table("GSE85241_cellsystems_dataset_4donors_updated.csv")  
pancreas_data_2

# cosine normalization here keeps the gene and cell names in the rows and columns

normalized_pancreas_2 = normalizing_cells_MaxAbsScaler(pancreas_data_2)
normalized_pancreas_2

# dropping missing values - looks like there aren't any here

normalized_pancreas_2.dropna(axis=1) 
normalized_pancreas_2