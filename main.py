#/usr/local/bin/python3.9
import pandas as pd
import networkx as nx
import numpy as np
from thinkstats2 import Pmf
from thinkstats2 import thinkplot
import matplotlib.pyplot as plt
import os
from scipy.stats import ttest_ind, ttest_ind_from_stats
from scipy.special import stdtr
print(os.chdir(os.path.dirname(__file__)))

def degrees(H):
    return [H.degree(u) for u in H]

# Uploading the PPI network file.
data1=pd.read_csv("Human-PPI.txt",delimiter='\t',header=0,names=['protein1','protein2'])
data1.to_csv("Human_network.csv")
interactions = data1[['protein1','protein2']]

# created the 1st network
G=nx.Graph(name='Protein Interaction Graph')

interactions = np.array(interactions)
for i in range(len(interactions)):
    interaction = interactions[i]
    a = interaction[0] # protein a node
    b = interaction[1] # protein b node
    G.add_edges_from([(a,b)])

print("Number of nodes=",G.number_of_nodes())# Number of nodes
print("Number of edges=",G.number_of_edges())# Number of edges
print("\nG.clustering",nx.clustering(G))# Clustering coefficient
print("\n G.average_clustering",nx.average_clustering(G))# Average clustering



# I have plotted two graph to show that degree distribution is scale free as nodes are lieing betweein 2 to 3.
degreespmf = Pmf(degrees(G))
thinkplot.Pdf(degreespmf, label='distribution')
plt.show()
H=G.number_of_nodes()
m=3 # M=3 because scale free graph values lies between 2 to 3
G = nx.barabasi_albert_graph(H, m)
degree_freq = nx.degree_histogram(G)
degrees = range(len(degree_freq))
plt.figure(figsize=(12, 8))
plt.loglog(degrees[m:], degree_freq[m:],'go-')
plt.xlabel('Degree')
plt.ylabel('Frequency')
plt.show()


# Uploaing 2 Proteins files.
data2=pd.read_csv("protein-list1.txt",delimiter='\t',header=None,names=['protein1'])
data3=pd.read_csv("protein-list2.txt",delimiter='\t',header=None,names=['protein2'])
data2_array= data2['protein1'].to_list()
data3_array=data3['protein2'].to_list()


# Created the second Graph
S=nx.Graph()
for i in data2_array:
    for j in data3_array:
      S.add_edges_from([(i,j)])

# Measuring the shortest path length of first graph G.
sp1=list()
length1 = dict(nx.all_pairs_shortest_path_length(G))
for key,value in length1.items():
    for x1,y1 in value.items():
       sp1.append(y1)
print(sp1)
# Measuring the shortest path lenth of second graph S
sp2=list()
length2= dict(nx.all_pairs_shortest_path_length(S))
for key1,value1 in length2.items():
    for x2,y2 in value1.items():
       sp2.append(y2)
print(sp2)

# Compared the Shortest path (SP1 from G graph) with Shortest path (Sp2 from S graph) value using t two sample test on unequal distribution.
t, p = ttest_ind(sp1, sp2, equal_var=False)
print("ttest_ind:t = %g  p = %g" % (t, p))





















