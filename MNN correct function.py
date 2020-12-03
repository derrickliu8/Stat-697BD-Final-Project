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
import pdb

#data1 = pd.read_table("GSE81076_D2_3_7_10_17.txt")
#data2 = pd.read_csv("GSE85241_cellsystems_dataset_4donors_updated.csv")

df1 = pd.DataFrame(np.random.randint(0,100,size=(15, 4)), columns=list('ABCD'))
df2 = pd.DataFrame(np.random.randint(0,100,size=(15, 4)), columns=list('ABCD'))


def my_kernel(X,Y):
    K = np.zeros((X.shape[0],Y.shape[0]))
    for i,x in enumerate(X):
        for j,y in enumerate(Y):
            K[i,j] = np.exp(-1*np.linalg.norm(x-y)**2)
    return K
    

def MNNcorrect(data1, data2, n_neighbors = 3):
    '''
    Takes in 2 or more datasets
    First dataset that is put in will be assumed to be the reference batch

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
    NN1 = NearestNeighbors(n_neighbors = n_neighbors, algorithm = "ball_tree").fit(data1_cnorm)
    NN2 = NearestNeighbors(n_neighbors = n_neighbors, algorithm = "ball_tree").fit(data2_cnorm)
    
    #Putting the nearest neighbors into an array
    #dist1 is the bacth correction vectors
    dist1, indices1 = NN1.kneighbors(data2_cnorm)
    dist2, indices2 = NN2.kneighbors(data1_cnorm)
    
    #Saving the bacth correction vectors in a dictionary to make it easy to access later
    #the keys in the first dictionary is the cell in the reference batch (row in W)
    #The keys in the nested dictionary are the cells in the secon batch (columns in w)
    #The values are the euclidean distance for that MNN
    correction_vectors = {}
    for key in range(0,len(indices1)):
        #breakpoint()
        correction_vectors[key] = {}
        dist = list(dist1[key])
        indices = list(indices1[key])
        for i in range(0,len(dist)):
            correction_vectors[key][indices[i]] = dist[i]
        
    
    graph1 = NN1.kneighbors_graph(data2_cnorm).toarray()
    graph2 = NN2.kneighbors_graph(data1_cnorm).toarray()
    
    #transposing the datasets to multiply them and find the MNN
    data2T = np.transpose(graph2)
    
    #W matrix of the indices of MNN, Windices are the MNN indices for each batch
    #[cell in reference batch, cell in second batch]
    W1 = np.multiply(data2T,graph1)
    Windices = np.transpose(((W1>0)).nonzero())
    
    #getting the correction vector for each MNN pair
    for i in Windices:
        one = i[0]
        two = i[1]
        #breakpoint()
        dist = correction_vectors[one][two] 

    K = my_kernel(data1_cnorm,data2_cnorm)


    return(correction_vectors, Windices, K)
   

    
    
        
