# -*- coding: utf-8 -*-
"""
Created on Thu May  7 18:11:33 2020

@author: Vaiva
"""


import networkx as nx
from collections import defaultdict
#from mlxtend.evaluate import permutation_test
from scipy.linalg import null_space
#from cycle_utils import *
import matplotlib.pyplot as plt
#from run_cycles_random_graph import random_dag, price_dag
from directed_cycle_basis import dag_cycle_basis_horton
from directed_cycle_basis import dag_cycle_basis_depina
import sys
import matplotlib.pyplot as plt
#sys.path.insert(1,"C://Users/Vaiva//antichains")
from alg_height import *

price = True
citation_net=False
foodweb = False
run = False



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



def find_cycles(graph):
    G = graph.copy()
    E = G.number_of_edges()
    tr(G)

    Gud = G.copy()
    Gud = Gud.to_undirected()
    C = nx.minimum_cycle_basis(Gud)
    return C

def print_cycle_statistics(graph,cycles=None, return_eigvals = False,ctype = "DAG"):
    G = graph.copy()
    E_0 = G.number_of_edges()
    tr(G)
    
    ###
    # Compute MCB
    ###
    if cycles == None and ctype =="notdag":
        print("No cycles found, they will be computed")
        E = G.number_of_edges()
        
    
        Gud = G.copy()
        Gud = Gud.to_undirected()
        C = nx.minimum_cycle_basis(Gud)
    elif cycles != None:
        C = cycles
    elif ctype =="horton":
        C = dag_cycle_basis_horton(graph,False)
    elif ctype =="depina":
        print("type=depina")
        C = dag_cycle_basis_depina(graph,False)
        
        #if C == None:
        #    print("Horton failed.")
        
         #   Gud = G.copy()
         #   Gud = Gud.to_undirected()
         #   C = nx.minimum_cycle_basis(Gud)
    
        
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
        #gc_edges_remove= [(u,v) for (u,v) in gc.edges() if u in sources]+[(u,v) for (u,v) in gc.edges() if v in sinks]
        #gc.remove_edges_from(gc_edges_remove)
        #sinks,sources = [n for n in gc if gc.out_degree(n)==0], [n for n in gc if gc.in_degree(n)==0]
        #print(sinks, sources)
        all_paths = []
        for src in sources:
            for snk in sinks:
                all_paths+=list(nx.all_simple_paths(gc,src,snk))
                
        all_paths = [len(x) for x in all_paths]
        #print(all_paths)
        #all_path_res = (np.mean(all_paths)-min(all_paths))/(max(all_paths)-min(all_paths))
        #if np.isnan(all_paths) == False:
        balance.append(np.std(all_paths)/np.mean(all_paths))#2*(np.mean(all_paths)-min(all_paths))/(max(all_paths)-min(all_paths)))
        balance_norm.append(np.std(all_paths)/max(all_paths))#2*(np.mean(all_paths)-min(all_paths))/(max(all_paths)-min(all_paths)))
        #except:
        #    balance.append(0)
        

    ###
    # Compute edge laplacian and cycle correlation matrix
    ###
    balance = [x for x in balance if np.isnan(x)==False]
    
    balance_norm = [x for x in balance_norm if np.isnan(x)==False]
    #resprint(balance)
    #print(cycle_size)
    M = np.array(C_vector)
    M_C = M.dot(M.T)
    M_E = M.T.dot(M)    
    M_C =np.absolute(M_C)
    S = sum([len(c) for c in C])
    edge_part = {e:0 for e in G.edges()}
    
    #name = randomString()
    
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
        
        
    lambda_max_a = np.linalg.eigvals(A_C)[0][eigvals.index(max(eigvals))]
    lc_eigvals = [eigval for eigval in np.linalg.eigvals(L_C) if eigval!=0]
    lc_eigvals = sorted(lc_eigvals,reverse= True)
    #print(lc_eigvals)
    try:
        lambda_min_l = lc_eigvals[0]-lc_eigvals[1]
    except:
        lambda_min_l = 0
    
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
    result["Number of edges"] = E_0
    result["Balance"] = np.mean(balance)
    result["Balance_norm"] = np.mean(balance_norm)
    result["Number_diamonds"] = len([x for x in cycle_type if x == "D"])
    result["Number_mixers"] = len([x for x in cycle_type if x == "M"])
    #result["Largest height antichain"] = get_largest_height_antichain(G)
    result["Largest eigenvalue M_C"] = np.real(max(eigvals))
    result["Largest eigenvalue A_C"] = np.real(lambda_max_a)
    result["Ratio eigenvalues"] = np.real(lambda_max_a/max(eigvals))
    #result["Spectral gap of L_C"] = np.real(lambda_min_l)
    result["Mean cycle height"] = np.mean([x[0] for x in cycle_height])#/ nx.dag_longest_path_length(G)
    result["Std cycle height"] = np.std([x[0] for x in cycle_height])#/ nx.dag_longest_path_length(G)    
    result["Mean stretch"] =np.mean(var_cycle_height)#max([x[0]/ nx.dag_longest_path_length(G) for x in cycle_height])-min([x[0]/ nx.dag_longest_path_length(G) for x in cycle_height])
    if return_eigvals ==True:
        return result,eigvals
    else:   
        return result


def get_largest_height_antichain(graph):
    H = set_heights(graph)
    antichains = {key:0 for key in H.values()}
    for key,val in H.items():
        antichains[val] +=1
    antichains = list(antichains.values())
    return max(antichains)

def map_cycles(cycles,graph):
    H  = set_heights(graph)
    cycle_height = []
    for c in cycles:
        h = [H[n] for n in c]
        cycle_height.append((np.mean(h),np.std(h)))
    plt.figure()
    plt.errorbar([a for a in range(len(cycles))],[c[0] for c in cycle_height],[c[1] for c in cycle_height])
    plt.ylim(0,max(H.values()))
    
    
def cycle_stats_plot(cycles,graph):
    H  = set_heights(graph)
    cycle_height = []
    for c in cycles:
        h = [H[n] for n in c]
        cycle_height.append((np.mean(h),np.std(h)))
            
    C_vector = []    
    edge_id = {}    
    i=0
    for e in G.edges():
        edge_id[e] =i 
        i+=1       
    for c in cycles:
        C_vector.append(cycle_basis_vector(G,c,edge_id))
    
    M = np.array(C_vector)
    M_C = M.dot(M.T)    
    M_C = abs(M_C)
    
    np.fill_diagonal(M_C,0)
    c_degree = []
    row_sum = np.sum(M_C,1)
    for i in range(len(cycles)):
        c_degree.append(np.sum(M_C,1)[i])
    c_size = [2**len(c) for c in cycles]
    plt.figure()
    plt.scatter([i for i in range(len(cycles))],[c[0] for c in cycle_height],s=c_size,c=c_degree)
    plt.colorbar()
    #plt.clim(0,10)
    #plt.ylim(0,50)
    return M_C
    
#__name__="bla"    
if __name__=="__main__":  
    #for n in [100]:#,150,200]:
    G =price_dag(100,3,1,5)
    tr(G)
    nremove = []
    
    G_root = G.copy()
    for n in G:
        if n not in nx.descendants(G,0) and n !=0:
            nremove.append(n)
    G_root.remove_nodes_from(nremove)
    edges = list(G_root.edges())
    G = nx.DiGraph()
    for n1,n2 in edges:
        G.add_edge(n1+1,n2+1)
    
    #G = nx.DiGraph()
    #G.add_edges_from([(1,2),(1,3),(2,4),(4,5),(3,5)])
    tr(G)
    #C = find_cycles(G)
    
    #cycle_basis,edge_id,CB = dag_cycle_basis_depina(G)
    res = print_cycle_statistics(G,None,ctype="depina")
    #map_cycles(C,G)
    #A = cycle_stats_plot(C,G)
    
