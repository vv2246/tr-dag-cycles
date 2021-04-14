#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 10:33:03 2020

@author: vvasiliau
"""


import networkx as nx
from collections import defaultdict
from scipy.linalg import null_space
import matplotlib.pyplot as plt
from alg_height import *
import itertools
import numpy as np
import random
import string
from numpy.linalg import lstsq
from scipy.linalg import orth
import scipy
from collections import defaultdict
from networkx.utils import not_implemented_for, pairwise
from networkx import minimum_cycle_basis
import math
import pandas as pd





def eigenspectrum(G,cycle_basis):
    
    C_vector = []
    edge_id = {}
    i=0
    id_edge = {}
    for e in G.edges():
        edge_id[e] =i 
        id_edge[i] = e
        i+=1
    E_0 = G.number_of_edges()
    for c in cycle_basis:
        gc = get_subgraph(G,c)
        vec = np.zeros(E_0)
        for e in gc.edges():
            vec[edge_id[e]] = 1
        C_vector.append(vec)
    M = np.array(C_vector)
    M_C = M.dot(M.T)
       
    return sorted(list(np.linalg.eigvals(M_C) ))

def _random_subset(seq,m):
    # """ Return m non-unique elements from seq.
    # This differs from random.sample which can return repeated
    # elements if seq holds repeated elements.
    # """
    targets=random.sample(seq,m)
    return targets



def complete_dag(n,c):
    
    #Create a DAG with n nodes in which an edge exists between i and j if i<j. 
    #Return a list of nodes in which each node i appears k_out(i) number of times.
    
    graph = nx.DiGraph()
    list_repeated_nodes  = []
    for i in range(n):
        graph.add_node(i)
        list_repeated_nodes.extend([i]*c)
        for j in range(i+1,n):
            graph.add_edge(i,j)
            list_repeated_nodes.append(i)
                       
    
    return list_repeated_nodes, graph

def lattice_dag_2D(N,L):
    #add nodes
    Gud = nx.grid_2d_graph(N,L)
    G = nx.DiGraph()
    
    for e1,e2 in Gud.edges():
        if e1[1]<=e2[1]:
            G.add_edge(e1,e2)
    return G



def russian_doll_dag(N):
    G = nx.DiGraph()
    G.add_edges_from([(1,2),(2,4),(1,3),(3,4),(7,1),(7,6),(6,5),(4,5)])
    n =8
    while G.number_of_nodes() <N:
        G.add_edge(n-3,n)
        G.add_edge(n+1,n)
        G.add_edge(n+2,n+1)
        G.add_edge(n+2,n-1)
        n +=3
        
    return G
        
        

def price_dag(n, m, c,delta ,seed=None):
    """Return random graph using Price cummulative advantage model.
    A graph of n nodes is grown by attaching new nodes each with m
    edges that are preferentially attached to existing nodes with high
    degree.
    
    Node t can only connect to nodes, that are t-delta or newer
    Parameters
    ----------
    n : int
        Number of nodes
    m : int
        Number of edges to attach from a new node to existing nodes
    c: number of times a node with out-degree=0 is added to the target list
        c = 1 Price original model
        c = m Directed Barabasi-Albert
    seed : int, optional
        Seed for random number generator (default=None).
    Returns
    -------
    G : Graph
    Notes
    -----
    The initialization is a graph with with m nodes and no edges.
    References
    ----------
    .. [1] de-Solla Price Network of Scientific Publications
    """

    if m < 1 or  m >=n:
        raise nx.NetworkXError(\
              "Price network must have m>=1 and m<n, m=%d,n=%d"%(m,n))
    if seed is not None:
        random.seed(seed)

    # Add m initial nodes (m0 in barabasi-speak)
    targets, G = complete_dag(m,c)
    
    Gnew = nx.DiGraph()
    for i, j in G.edges():
        Gnew.add_edge(i+1,j+1)
    G = Gnew.copy()
    targets = [i+1 for i in targets]
    # List of existing nodes, with nodes repeated once for each adjacent edge
    repeated_nodes=[]
    # Start adding the other n-m nodes. The first node is m.
    source=m+1
    while source<n+1:
        # Add edges to m nodes from the source.
        G.add_edges_from(zip(targets,[source]*m))
        # Add one node to the list for each new edge just created.
        repeated_nodes.extend(targets)
        # And the new node "source" has c times  to add to the list.
        
        repeated_nodes.extend([source]*c)
        # Now choose m unique nodes from the existing nodes
        # Pick uniformly from repeated_nodes (preferential attachement)
        targets = _random_subset(repeated_nodes,m)
        source += 1
    return G



def random_dag(N,P):
    nodes = [n for n in range(1,N+1)]
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    for n1,n2 in itertools.combinations(nodes,2):
        p = random.random()
        if p <= P:
            if n1 > n2:
                G.add_edge(n2,n1)
            else:
                G.add_edge(n1,n2)
    return G
       


def tr(DAG, output=False):
    # for printing progress
    E = DAG.number_of_edges()
    i = 0
    print_limit = 10
    print_counter = print_limit
    edges = list(DAG.edges())
    #########################
    for edge in edges:
        # check edge is necessary for causal structure
        [a, b] = edge
        DAG.remove_edge(a, b)
        if not nx.has_path(DAG, a, b):
            DAG.add_edge(a, b)
        
        if output:
            i += 1
            pc = (i/float(E))*100
            if pc > print_counter:
                #print ('Finished %s percent' % int(math.floor(pc)))
                print_counter += print_limit
                
    return DAG    




def print_cycle_statistics(graph,cycles=None, return_eigvals = False):
    G = graph.copy()
    E_0 = G.number_of_edges()
    tr(G)
    
    ###
    # Compute MCB
    ###
    if cycles == None :
        print("No MCB found, computing MCB")
        #E = G.number_of_edges()
        Gud = G.copy()
        Gud = Gud.to_undirected()
        C = minimum_cycle_basis(Gud)
        print("Betti 1= ",len(C), "Observed MCB=" ,G.number_of_edges()-G.number_of_nodes()+1)
    elif cycles != None:
        C = cycles
        
        
        
    ###
    # Find cycle vectors
    ###
    C_vector = []
    edge_id = {}
    i=0
    id_edge = {}
    for e in G.edges():
        edge_id[e] =i 
        id_edge[i] = e
        i+=1
    
    balance_norm = []
    balance = []
    cycle_size = []
    cycle_type = []
    for c in C:
        gc = get_subgraph(G,c)
        #print(gc.nodes())
        cycle_size.append(gc.number_of_nodes())
        vec = np.zeros(E_0)
        for e in gc.edges():
            vec[edge_id[e]] = 1
        C_vector.append(vec)#cycle_basis_vector(G,c,edge_id))
        
        
        ###
        ## Find subgraph longest, mean, shortest path
        ##
        sinks,sources = [n for n in gc if gc.out_degree(n)==0], [n for n in gc if gc.in_degree(n)==0]
        if len(sinks)>1:
            cycle_type.append("M")
        else:
            cycle_type.append("D")
        all_paths = []
        for src in sources:
            for snk in sinks:
                all_paths+=list(nx.all_simple_paths(gc,src,snk))
                
        all_paths = [len(x) for x in all_paths]
        balance.append(np.std(all_paths)/np.mean(all_paths))#2*(np.mean(all_paths)-min(all_paths))/(max(all_paths)-min(all_paths)))
        balance_norm.append(np.std(all_paths)/max(all_paths))#2*(np.mean(all_paths)-min(all_paths))/(max(all_paths)-min(all_paths)))
    ###
    # Compute edge laplacian and cycle correlation matrix
    ###
    balance = [x for x in balance if np.isnan(x)==False]
    
    balance_norm = [x for x in balance_norm if np.isnan(x)==False]
    M = np.array(C_vector)
    M_C = M.dot(M.T)
    M_E = M.T.dot(M)    
    M_C =np.absolute(M_C)
    S = sum([len(c) for c in C])
    edge_part = {e:0 for e in G.edges()}
    ###
    # Compute edge participation
    ###
    for e1,e2 in G.edges():
        for c in C:
            if e1 in c and e2 in c:
                edge_part[(e1,e2)] += 1
    
    ###
    # Compute largest effective cycle
    ###
    try:
        eigvals = np.linalg.eig(M_C)[0]
        x = [i/len(C) for i in range(len(C))]         
    except:
        pass
    
    A_C = M_C.copy()
    L_C = M_C.copy()
    for i in range(len(A_C)):
        A_C[i,i] = 0
    
    for i in range(len(A_C)):
        L_C[i,i] = -np.sum(A_C[i],axis =0)   
    L_C = L_C*-1
        
    print(np.mean(eigvals),np.mean([len(c) for c in C]))
    lambda_max_a =np.linalg.eigvals(A_C)[list(eigvals).index(max(eigvals))]
    
    ###
    # Compute cycle heights
    ###
    
    H  = set_heights(graph)
    cycle_height = []
    var_cycle_height = []
    for c in C:
        h = [H[n] for n in c]
        cycle_height.append((np.mean(h),np.std(h)))
        var_cycle_height.append(max(h)-min(h))
        
    result = {}
    result["Mean edge participation"] = np.mean(list(edge_part.values()))
    result["Std edge participation"] = np.std(list(edge_part.values()))
    result["Largest cycle size"] = max([len(c) for c in C])
    result["Mean cycle size"] = np.mean([len(c) for c in C])
    result["Std cycle size"] = np.std([len(c) for c in C])
    result["Number of cycles"] = len(C)
    result["Number cycle connected components"] = (null_space(L_C)).shape[1]
    result["Number of TR edges"] = G.number_of_edges()
    result["Longest path"] = nx.dag_longest_path_length(G)
    result["Number of nodes"] = G.number_of_nodes()
    #result["Number of edges"] = E_0
    result["Balance"] = np.mean(balance)
    #result["Balance_norm"] = np.mean(balance_norm)
    result["Number_diamonds"] = len([x for x in cycle_type if x == "D"])
    result["Number_mixers"] = len([x for x in cycle_type if x == "M"])
    result["Largest eigenvalue M_C"] = np.real(max(eigvals))
    #result["Largest eigenvalue A_C"] = np.real(max(np.linalg.eigvals(A_C)))
    result["Mean cycle height"] = np.mean([x[0] for x in cycle_height])#/ nx.dag_longest_path_length(G)
    result["Std cycle height"] = np.std([x[0] for x in cycle_height])#/ nx.dag_longest_path_length(G)    
    result["Mean stretch"] =np.mean(var_cycle_height)#max([x[0]/ nx.dag_longest_path_length(G) for x in cycle_height])-min([x[0]/ nx.dag_longest_path_length(G) for x in cycle_height])
    result["Std stretch"] =np.std(var_cycle_height)#max([x[0]/ nx.dag_longest_path_length(G) for x in cycle_height])-min([x[0]/ nx.dag_longest_path_length(G) for x in cycle_height])
    
    #result["eigenvalue_ratio"] = np.real(lambda_max_a)/np.real(max(eigvals))
    result["Largest M_C by corr size"] = np.real(max(eigvals))/[len(c) for c in C][list(eigvals).index(max(eigvals))]
    if return_eigvals ==True:
        return result,eigvals
    else:   
        return result



     
def get_subgraph(G,nodelist):
    #Find the subgraph that only contains nodes from the nodelist and edges, only adjacent to those nodes
    D = nx.DiGraph()
    D.add_nodes_from(nodelist)
    for n1,n2 in itertools.combinations(nodelist,2):
        if (n1,n2) in G.edges():
            D.add_edge(n1,n2)
        elif (n2,n1) in G.edges():
            D.add_edge(n2,n1)
    remove = []
    for n in D.nodes():
        if D.in_degree(n)+ D.out_degree(n)<2 :
            remove.append(n)
       
            
    D.remove_nodes_from(remove)
    return D



