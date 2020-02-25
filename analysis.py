
from math import log
from copy import copy, deepcopy
from operator import itemgetter
from scipy import stats
import pandas as pd
import numpy as np
import networkx as nx
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import os
from random import sample
from import_and_clean import treated_df, control_df, network_control, network_treated# control_df_preabs
import time

#Old implementation with for loop instead of adding the dictionary all at once
#def add_expression(df, protein_network):
 #   for i in np.arange(0, len(df)):
  #      gene = df.Gene_Name.iloc[i]
   #     experiment = df.experiment.iloc[i]
    #    if gene in list(network):
     #       protein_network.node[gene][experiment] = float(df.expression.iloc[i])
   # return(protein_network)

##This section of the code attaches expression values for each gene to the relevant node of the graph for each cell_line

def add_expression(df, protein_network):
    df_small = df[['Gene_Name', 'experiment','expression']]
    df_wide = df_small.pivot(index = 'experiment', columns = 'Gene_Name', values = 'expression') #This is necessary to get the dictionary structured right for the graph
    dict_target = df_wide.to_dict()
    nx.set_node_attributes(protein_network, dict_target)
    return(protein_network)

network_control = add_expression(df = control_df, protein_network = network_control)
#network_treated = add_expression(df = treated_df, protein_network = network_treated)


def gibbs_calc(G, cell_line, df): 

    cell_line_df = df[df.experiment == cell_line] #isolate observations for the cell_line in question
    target_genes = list(cell_line_df.Gene_Name) 
    nbhd_conc = 'nbhd_conc-' + str(cell_line) # name the neighborhood concentration with the cell line name
    nx.set_node_attributes(G, name=nbhd_conc, values=0) # initialize the nbhd_conc variable for each node
    gibbs = 'gibbs-' + str(cell_line) # name the gibbs with the cell line name
    nx.set_node_attributes(G, name=gibbs, values=0) # initialize the gibbs variable for each node
    
    for node in list(G): # iterate through the nodes of G
        if G.node[node][cell_line] == 0:
            G.node[node][gibbs] == 0
        else:
            for neighbor in G.neighbors(node): # iterate through the neighbors of each node
                G.node[node][nbhd_conc] += G.node[neighbor][cell_line]
            try:
                G.node[node][gibbs] = G.node[node][cell_line] * log(G.node[node][cell_line] / (G.node[node][nbhd_conc] + G.node[node][cell_line]))
            except ValueError:
                print(cell_line, node, neighbor)
                print(G.node[node][cell_line], G.node[node][nbhd_conc])

#def gibbs_calc(G, cell_line, df): ##Okay so the other one only takes 15 seconds now? what the fuck. 
 #   cell_line_df = df.loc[df.experiment == cell_line]
  #  nbhd_conc = 'nbhd_conc-' + str(cell_line) # name the neighborhood concentration with the cell line name
   # nx.set_node_attributes(G, name=nbhd_conc, values=0) # initialize the nbhd_conc variable for each node
    #gibbs = 'gibbs-' + str(cell_line) # name the gibbs with the cell line name
    #nx.set_node_attributes(G, name=gibbs, values=0) # initialize the gibbs variable for each node - hoping to do this later
    
    #for node in list(G): # iterate through the nodes of G
     #   ##nbhd concentration calculation
      #  neighbors = list(G[node])
       # neighbors_df = cell_line_df.loc[cell_line_df.Gene_Name.isin(neighbors)]
        #G.node[node][nbhd_conc] = sum(neighbors_df.expression)
             
    #I know this looks unweildy but its all vectorised so it should happen quick 
   # nbhd_conc = nx.get_node_attributes(G, nbhd_conc)
   # nbhd_conc_df = pd.DataFrame.from_dict(nbhd_conc, orient = "index", columns = ["nbhd_conc"])
   # nbhd_conc_df.index.name = "Gene_Name"
   # nbhd_conc_df.reset_index(inplace = True)
   # cell_line_df = pd.merge(cell_line_df, nbhd_conc_df, on =("Gene_Name"))
   # cell_line_df['gibbs'] = list(cell_line_df['expression']) * np.log(cell_line_df['expression'] / (cell_line_df['nbhd_conc'] + cell_line_df['expression']))
   #return cell_line_df


##This section runs the function and then binds the dataframes together - the attempted faster version that ended up way slower somehow. 
#tic = time.process_time()
#gibbs_df = []
#for experiment in np.unique(control_df.experiment): 
 #   cell_line_df = gibbs_calc(network, experiment, control_df)
#    gibbs_df.append(cell_line_df)
#gibbs_df = pd.concat(gibbs_df)
#gibbs_df.reset_index(inplace = True)
#gibbs_df.drop('index', axis = 1, inplace = True)
#toc = time.process_time()

tic = time.process_time()

#Calculate gibbs for all nodes and then play with the data a bit.
for experiment in np.unique(control_df.experiment): 
    gibbs_calc(network_control, experiment, control_df)

#for experiment in np.unique(treated_df.experiment):
 #   gibbs_calc(network_treated, experiment, treated_df)
toc = time.process_time()   
print(toc-tic)
gibbs_dict_control = {}
#gibbs_dict_treated = {}
for experiment in np.unique(control_df.experiment):
    gibbses_control = nx.get_node_attributes(network_control, 'gibbs-' + str(experiment))
    
    gibbs_dict_control[experiment] = gibbses_control
    
#for experiment in np.unique(treated_df.experiment):
 #   gibbses_treated = nx.get_node_attributes(network_treated, 'gibbs-' + str(experiment))
  #  gibbs_dict_treated[experiment] = gibbses_treated

gibbs_df_control = pd.DataFrame.from_dict(gibbs_dict_control)
gibbs_df_control.reset_index(inplace = True)
gibbs_df_control = gibbs_df_control.rename(columns = {'index':'Gene_Name'})
gibbs_df_control = pd.melt(gibbs_df_control, id_vars= 'Gene_Name', var_name = 'experiment', value_name = 'gibbs')
#gibbs_df_treated = pd.DataFrame.from_dict(gibbs_dict_treated)
#gibbs_df_treated.reset_index(inplace = True)
#gibbs_df_treated = gibbs_df_treated.rename(columns = {'index':'Gene_Name'})
#gibbs_df_treated = pd.melt(gibbs_df_treated, id_vars= 'Gene_Name', var_name = 'experiment', value_name = 'gibbs')

#Need to make a subnetwork to highlight
nodes = sample(list(network['EGFR']), k = 10)
nodes.append('EGFR')
nodes2 = sample(list(network[nodes[0]]), k = 10)
nodes = nodes + nodes2
subnetwork = network.subgraph(nodes)
##This is for graphing the first network
nx.draw_networkx(subnetwork, nodelist = nodes) #view it by using plt.show()
#plt.clf() ##Clear this figure to allow the next one to be made

##This block of code draws a subnetwork with the 10 nodes with the largest gibbs values (max of all cell lines) -hmm, I think this is just for the sub_network.
#nodes = [x[0] for x in max_gibbs[:10]]
#sub_G = network.subgraph(nodes)
#nx.draw_networkx(sub_G)
#plt.clf() 

##I think this is the function I have to fuck with to evaluate the impact of removing multiple nodes at once a lah miRNA targeting.
##Dictionaries seem akin to vectors of factors in R
##Dictionary elements are accessed via keys
def Gibbs_Calc_Diff(G, cell_line, proportion_target = 0.05, df = gibbs_df_control):
    gibbs_diff_dict = {}
    df_cell_line = df[df.experiment == cell_line]
    target_genes = list(df_cell_line.Gene_Name)
    #start by calculating targets for this cell_line
    def calculate_targets(data = df_cell_line, proportion_to_target = proportion_target):
        #Take top X % of gibbs free energies
        df_sorted = data.sort_values(by = 'gibbs', ascending = True)
        df_targets = df_sorted.iloc[np.arange(0,proportion_to_target*len(df_sorted))]
        targets = np.unique(df_targets.Gene_Name)
        return targets
    targets = calculate_targets(data = df_cell_line, proportion_to_target = proportion_target)
    gibbs1 = sum(df_cell_line.gibbs) ##calculate gibbs 1
    #Iterate through targets, removing nodes
    for target in targets:
        ##First need to calculate gibbs free energy pre- removal 
        old_value = G.node[target][cell_line]  # bad things happen when we try to remove nodes from the network over which we're iterating - this step is slow
        G.node[target][cell_line] = 0
        #have to initialize these
        nbhd_conc = 'nbhd_conc-' + str(cell_line)
        gibbs = 'gibbs-' + str(cell_line)
        nx.set_node_attributes(G, name=nbhd_conc, values=0) # initialize the nbhd_conc 
        nx.set_node_attributes(G, name=gibbs, values=0) #initialize gibbs again
        for node in list(G): # iterate through the nodes of G
            if G.node[node][cell_line] == 0:
                G.node[node][gibbs] == 0
            else:
                for neighbor in G.neighbors(node): # iterate through the neighbors of each node
                    G.node[node][nbhd_conc] += G.node[neighbor][cell_line]
                #Calculate new gibbs
                G.node[node][gibbs] = G.node[node][cell_line] * np.log(G.node[node][cell_line] / (G.node[node][nbhd_conc] + G.node[node][cell_line]))
            
        ##Here is where we take the new network and calculate the new gibbs free energy    
        gibbses = nx.get_node_attributes(G, 'gibbs-' + str(cell_line))
        gibbs2= sum(gibbses.values())
        gibbs_diff_dict[target, cell_line] = gibbs2 - gibbs1
        G.node[target][cell_line] = old_value
    return gibbs_diff_dict

#debugging
#test = Gibbs_Calc_Diff(G = network, cell_line = control_df.experiment[0], df = gibbs_df, proportion_target = 0.0005)

gibb_diff_list = []
for experiment in np.unique(control_df.experiment):
   gibb_diff_list.append(Gibbs_Calc_Diff(G = network_control, cell_line = experiment, df = gibbs_df_control, proportion_target = 0.05))
gibb_df = []
for i in np.arange(0, len(gibb_diff_list)):
    gibb_diff = gibb_diff_list[i]
    gibb_df1 = pd.DataFrame.from_dict(gibb_diff, orient = 'index').reset_index()
    gibb_df.append(gibb_df1)
gibb_df = pd.concat(gibb_df)
gibb_df.to_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/gibbs_diff_values_10_23_2019.csv')
gibb_df = pd.read_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/gibbs_diff_values_10_23_2019.csv')
#string_obj = gibb_df['index'].str.replace(pat = "'", repl = "")
#string_obj = string_obj.str.replace(pat = "(", repl = "").str.replace(pat = ")", repl = "")
#string_obj = string_obj.str.split(pat = ", ")

#for i in np.arange(0,len(gibb_df)):
 #   gibb_df['Gene_Name'][i] = string_obj[i][0]
  #  gibb_df['experiment'][i] = string_obj[i][1]

#First need to convert gibb_df to a dataframe
gibb_df = pd.DataFrame(gibb_df) 
#Need to split cell line and Gene Name
#Get rid of the extra characters and split the character vector up at the comma that divides the cell line and gene name
string_obj = gibb_df['index'].str.replace(pat = "(", repl = "")
string_obj = string_obj.str.replace(pat = ")", repl = "").str.replace(pat = "\'", repl = "").str.split(pat = ", ")

experiment_list = []
gene_name_list = []
#Add the variables to the original dataframe 
for i in np.arange(0,len(gibb_df)):
    experiment_list.append(string_obj[i][1])
    gene_name_list.append(string_obj[i][0])
gibb_df['experiment'] = experiment_list
gibb_df['Gene_Name'] = gene_name_list

gibb_df = gibb_df.drop(["Unnamed: 0", "index"], axis = 1)
gibb_df = gibb_df.rename(columns ={"0":"gibbs_diff"})
#split experiment as well
string_obj = gibb_df['experiment'].str.replace(pat = "_DMSO", repl = "")
string_obj = string_obj.str.split(pat = "_")
experiment_list = []
cell_line_list = []
for i in np.arange(0,len(gibb_df)):
    experiment_list.append(string_obj[i][0])
    cell_line_list.append(string_obj[i][1])

gibb_df['experiment'] = experiment_list
gibb_df['cell_line'] = cell_line_list


#Prepare to merge by doing the same string work on control_df and gibbs_df
gibbs_df_control['condition'] = "control"
#gibbs_df_treated['condition'] = "treated"
#gibbs_df = gibbs_df_control.append(gibbs_df_treated)
#gibbs_df.reset_index(inplace = True)
gibbs_df = gibbs_df_control
string_obj = gibbs_df['experiment'].str.replace(pat = "_DMSO", repl = "")
string_obj = string_obj.str.split(pat = "_")
experiment_list = []
cell_line_list = []
for i in np.arange(0,len(string_obj)):
    experiment_list.append(string_obj[i][0])
    cell_line_list.append(string_obj[i][1])

gibbs_df['experiment'] = experiment_list
gibbs_df['cell_line'] = cell_line_list

gibbs_df = gibbs_df.drop("condition", axis = 1)

#Do the same text processing on the control dataset
control_df = control_df.reset_index()
control_df = control_df.drop(columns = ["index"])
string_obj = control_df['experiment'].str.replace(pat = "_DMSO", repl = "")
string_obj = string_obj.str.split(pat = "_")
experiment_list = []
cell_line_list = []
for i in np.arange(0,len(control_df)):
    experiment_list.append(string_obj[i][0])
    cell_line_list.append(string_obj[i][1])

control_df['experiment'] = experiment_list
control_df['cell_line'] = cell_line_list

control_df = control_df.drop("condition", axis= 1)


#And on the treated Dataset -  a lot of what is commented out is a holdover from when we were mucking around with both the control and treatment data. 
#treated_df = treated_df.reset_index()
#treated_df = treated_df.drop(columns = ["index"])
#string_obj = treated_df['experiment'].str.replace(pat = "_DMSO", repl = "")
#string_obj = string_obj.str.split(pat = "_")
#experiment_list = []
#cell_line_list = []
#for i in np.arange(0,len(treated_df)):
 #   experiment_list.append(string_obj[i][0])
  #  cell_line_list.append(string_obj[i][1])

#treated_df['experiment'] = experiment_list
#treated_df['cell_line'] = cell_line_list
#expression_df = control_df.append(treated_df)
#expression_df.reset_index(inplace = True)

#drop index column in Gibbs_df
#gibbs_df = gibbs_df.drop(columns = ['index'])
#expression_df = expression_df.drop(columns = ['index'])

#Merge it up - this step is messed up somehow
final_df = pd.merge(control_df, gibbs_df, on = ["Gene_Name", "experiment", "cell_line"])
final_df.to_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/EwingsDatasetFullClean_11_12_2019.csv')
final_df = pd.merge(gibb_df, final_df, on = ["Gene_Name", "experiment", "cell_line"])
final_df.to_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/EwingsDatasetFinal_11_12_2019.csv')
##what is the point of this 100
#number_target = 100 ## make this top 5% of all genes
#betti_dict = {}
#okay so we only want to try removing the top 100 in max gibbs

##Looks like we are only calculating betti number on the first network? is this just doing the max of the cell lines
#max_gibbs = dict(max_gibbs)
#target_nodes = list(max_gibbs.keys())[0:number_target]
#G_subbed = deepcopy(network)
#G_subbed.remove_nodes_from([n for n in network if n not in set(nodes)])
#bettis = bettis_calc(network, targets = target_nodes)
#sorted_bettis = sorted(bettis.items(), key=itemgetter(1), reverse=True)

#Okay so this is taking only genes that are in the upper 50% of all genes (because it is after z-score normalization) -it is also doing it for pre_abs, not pre_absdf

##Why does it have to be upregulated in every cell line
#only_upregulated = [] Need to compare to normal tissue for only upregulated - ask jessica
#for entry in sorted_bettis:
#    if (control_df_preabs.loc[entry[0]] > 0).all(0):
#        only_upregulated.append(entry)
#only_upregulated

#(pre_abs.loc['TP53'] > 0).any(0)

#with open('ranked_betti_list.csv', 'w', newline='') as csvfile:
   # writer = csv.writer(csvfile)
   # writer.writerows(only_upregulated)
