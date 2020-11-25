# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 16:41:24 2020

@author: annab
"""
import sklearn
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import Normalizer

data1 = pd.read_table("GSE81076_D2_3_7_10_17.txt")
data2 = pd.read_csv("GSE85241_cellsystems_dataset_4donors_updated.csv")


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
    data1 = data1.rename(columns = {"Unnamed: 0":"Genes"}) 
    data1.set_index(["Genes"], inplace = True)
    data2 = data2.rename(columns = {"Unnamed: 0":"Genes"}) 
    data2.set_index(["Genes"], inplace = True)   
    
    transformer = Normalizer().fit(data1g)
    data1_cnorm = pd.DataFrame(transformer.transform(data1g))
    
    
    return(data1_cnorm,data2)
   
 
    
    
        
