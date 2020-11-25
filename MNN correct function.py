# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 16:41:24 2020

@author: annab
"""
import sklearn
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import Normalizer
import numpy as np

#data1 = pd.read_table("GSE81076_D2_3_7_10_17.txt")
#data2 = pd.read_csv("GSE85241_cellsystems_dataset_4donors_updated.csv")

df1 = pd.DataFrame(np.random.randint(0,100,size=(15, 4)), columns=list('ABCD'))
df2 = pd.DataFrame(np.random.randint(0,100,size=(15, 4)), columns=list('ABCD'))

def MNNcorrect(data1, data2):
    '''
    Takes in 2 or more datasets

    Returns one dataset that has been corrected using the Mutual Nearest Neighbors
    batch correction algorithm
    '''    
 #Checking the datasets
    if len(data1) == 1:
        raise TypeError("Dataset needs to have more than 1 row")
    if len(data1) == 0:
        raise TypeError("Dataset is empty")
 #Renaming the datasets so that the genes are row names       
    #data1 = data1.rename(columns = {"Unnamed: 0":"Genes"}) 
    #data1.set_index(["Genes"], inplace = True)
    #data2 = data2.rename(columns = {"Unnamed: 0":"Genes"}) 
    #data2.set_index(["Genes"], inplace = True)   
    
    #Cosine normalizing the data
    transformer = Normalizer().fit(data1)
    data1_cnorm = pd.DataFrame(transformer.transform(data1))
    
    transformer2 = Normalizer().fit(data2)
    data2_cnorm = pd.DataFrame(transformer2.transform(data2))
    
    #Performing the Nearest Neighbors algorthim
    NN1 = NearestNeighbors(n_neighbors = 2, algorithm = "kd_tree").fit(data1_cnorm)
    NN2 = NearestNeighbors(n_neighbors = 2, algorithm = "kd_tree").fit(data2_cnorm)
    
    #Putting the nearest neighbors into an array
    dist1, indices1 = NN1.kneighbors(data1_cnorm)
    dist2, indices2 = NN2.kneighbors(data2_cnorm)
    
    graph1 = NN1.kneighbors_graph(data1_cnorm).toarray()
    graph2 = NN2.kneighbors_graph(data2_cnorm).toarray()
    
    #transposing the datasets to multiply them and find the MNN
    data1T = np.transpose(graph1)
    data2T = np.transpose(graph2)
    
    #W matrix of the indices of MNN
    #Still working on this
    W1 = np.concatenate((data1T, data2T))
    WT = np.transpose(W1)
    W = W1*WT
    
    
    return(W1)
   
 
    
    
        
