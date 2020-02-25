##Version 0.1: Going to use all the genes in the dataset rather than the ones jessica gave us - still unsure what those pertain to at this point 
from math import log
from copy import copy, deepcopy
from operator import itemgetter
import pandas as pd
import numpy as np
import networkx as nx
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy import stats
default_path = "C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files"
os.chdir(default_path)

#Load and clean data
df = pd.read_csv('rld_Counts.csv')

def clean_data(data): 
    data = data.rename(columns = {"Unnamed: 0":"Gene_Name"})
    
    # get rid of rows where any experiment has 0 expression - this is a bit of an issue
    #data = data[~data.eq(0).any(1)]
    
    # drop NaNs
    data = data.dropna()
    
    #Lets reshape the data for easier analysis 
    data = pd.melt(data, id_vars = ["Gene_Name"])
    data = data.rename(columns = {"variable":"experiment", "value":"expression"})
    def label_treated(row):
        if row['experiment'].endswith('HCI_2509'):
            return "treated"
        if row['experiment'].endswith('DMSO'):
            return "control"
    data['condition'] = data.apply(label_treated, axis = 1)
    return(data)

cleaned_df = clean_data(data = df) 

##Z-score normalization
 
#def normalize_data(data): 
 #   grouped = data.groupby('Gene_Name')
 #   transformed = grouped.transform(lambda x: (x - x.mean()) / x.std())
 #   data['Z_score_abs'] = abs(transformed)
 #   data['Z_score'] = transformed
 #   return(data)

cleaned_df.loc[cleaned_df.expression < 0, 'expression'] = 0 #Here we are assuming that everything less than 0 in the rlog parlance is equal to 0. this is probably not ideal. It'll do for now tho
##break out by treated/control cases - just necessary because of how the downstream code was written
control_df = cleaned_df.loc[cleaned_df.condition == 'control']
treated_df = cleaned_df.loc[cleaned_df.condition == 'treated']
    
##Protein interaction network data
biogrid = pd.read_table('BIOGRID-ORGANISM-Homo_sapiens-3.5.171.tab2.txt',
                        usecols=['Official Symbol Interactor A', 'Official Symbol Interactor B'])


def clean_biogrid(expression_df, grid_df):
    gene_list = np.unique(np.concatenate((pd.unique(biogrid['Official Symbol Interactor A']), 
                                  pd.unique(biogrid['Official Symbol Interactor B']))))
    genes_with_expression = pd.unique(expression_df.Gene_Name)
    genes_in_both = set(genes_with_expression).intersection(set(gene_list))
    
    # drop unnecessary edges
    grid_df = grid_df[grid_df.isin(genes_in_both).all(1)]

    # construct the network
    G = nx.from_pandas_edgelist(grid_df, source='Official Symbol Interactor A', target='Official Symbol Interactor B')
    G = G.to_undirected()
    G.remove_edges_from(G.selfloop_edges())
    
    # drop unnecessary genes
    ##"""Note that the list of nodes will be smaller than genes_in_both because we restrict G to
    #contain only edges where BOTH nodes are in the list, not just one.""" 
    
    expression_df = expression_df[expression_df.Gene_Name.isin(list(G))]
    #print(len(df.index.values))
    return(G, expression_df)
network_control, control_df = clean_biogrid(control_df, biogrid)
network_treated, treated_df = clean_biogrid(treated_df, biogrid)