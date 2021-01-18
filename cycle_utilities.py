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

import math





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
        
        
        
    


def minimum_cycle_basis(G, weight=None):
    """ Returns a minimum weight cycle basis for G

    Minimum weight means a cycle basis for which the total weight
    (length for unweighted graphs) of all the cycles is minimum.

    Parameters
    ----------
    G : NetworkX Graph
    weight: string
        name of the edge attribute to use for edge weights

    Returns
    -------
    A list of cycle lists.  Each cycle list is a list of nodes
    which forms a cycle (loop) in G. Note that the nodes are not
    necessarily returned in a order by which they appear in the cycle

    Examples
    --------
    >>> G=nx.Graph()
    >>> nx.add_cycle(G, [0,1,2,3])
    >>> nx.add_cycle(G, [0,3,4,5])
    >>> print([sorted(c) for c in nx.minimum_cycle_basis(G)])
    [[0, 1, 2, 3], [0, 3, 4, 5]]

    References:
        [1] Kavitha, Telikepalli, et al. "An O(m^2n) Algorithm for
        Minimum Cycle Basis of Graphs."
        http://link.springer.com/article/10.1007/s00453-007-9064-z
        [2] de Pina, J. 1995. Applications of shortest path methods.
        Ph.D. thesis, University of Amsterdam, Netherlands

    See Also
    --------
    simple_cycles, cycle_basis
    """
    # We first split the graph in commected subgraphs
    return sum((_min_cycle_basis(G.subgraph(c), weight) for c in
                nx.connected_components(G)), [])



def _min_cycle_basis(comp, weight):
    cb = []
    # We  extract the edges not in a spanning tree. We do not really need a
    # *minimum* spanning tree. That is why we call the next function with
    # weight=None. Depending on implementation, it may be faster as well
    spanning_tree_edges = list(nx.minimum_spanning_edges(comp, weight=None,
                                                         data=False))
    edges_excl = [frozenset(e) for e in comp.edges()
                  if e not in spanning_tree_edges]
    N = len(edges_excl)

    # We maintain a set of vectors orthogonal to sofar found cycles
    set_orth = [set([edge]) for edge in edges_excl]
    for k in range(N):
        # kth cycle is "parallel" to kth vector in set_orth
        new_cycle = _min_cycle(comp, set_orth[k], weight=weight)
        cb.append(list(set().union(*new_cycle)))
        # now update set_orth so that k+1,k+2... th elements are
        # orthogonal to the newly found cycle, as per [p. 336, 1]
        base = set_orth[k]
        set_orth[k + 1:] = [orth ^ base if len(orth & new_cycle) % 2 else orth
                            for orth in set_orth[k + 1:]]
    return cb


def _min_cycle(G, orth, weight=None):
    """
    Computes the minimum weight cycle in G,
    orthogonal to the vector orth as per [p. 338, 1]
    """
    T = nx.Graph()

    nodes_idx = {node: idx for idx, node in enumerate(G.nodes())}
    idx_nodes = {idx: node for node, idx in nodes_idx.items()}

    nnodes = len(nodes_idx)

    # Add 2 copies of each edge in G to T. If edge is in orth, add cross edge;
    # otherwise in-plane edge
    for u, v, data in G.edges(data=True):
        uidx, vidx = nodes_idx[u], nodes_idx[v]
        edge_w = data.get(weight, 1)
        if frozenset((u, v)) in orth:
            T.add_edges_from(
                [(uidx, nnodes + vidx), (nnodes + uidx, vidx)], weight=edge_w)
        else:
            T.add_edges_from(
                [(uidx, vidx), (nnodes + uidx, nnodes + vidx)], weight=edge_w)
    
    #all_shortest_pathlens = dict(nx.shortest_path_length(T, weight=weight))
    #print("size UD ASP_LEN",len(all_shortest_pathlens))
    #cross_paths_w_lens = {n: all_shortest_pathlens[n][nnodes + n]
     #                     for n in range(nnodes)}
    
    cross_paths_w_lens = {}
    for n in range(nnodes):
        cross_paths_w_lens[n] = len(nx.shortest_path(T,n,n+nnodes,weight="weight"))
    #print(all_shortest_pathlens)
    # Now compute shortest paths in T, which translates to cyles in G
    start = min(cross_paths_w_lens, key=cross_paths_w_lens.get)
    end = nnodes + start
    min_path = nx.shortest_path(T, source=start, target=end, weight='weight')
    #print(min_path)
    # Now we obtain the actual path, re-map nodes in T to those in G
    min_path_nodes = [node if node < nnodes else node - nnodes
                      for node in min_path]
    # Now remove the edges that occur two times
    mcycle_pruned = _path_to_cycle(min_path_nodes)

    return {frozenset((idx_nodes[u], idx_nodes[v])) for u, v in mcycle_pruned}


def _path_to_cycle(path):
    """
    Removes the edges from path that occur even number of times.
    Returns a set of edges
    """
    edges = set()
    for edge in pairwise(path):
        # Toggle whether to keep the current edge.
        edges ^= {edge}
    return edges

def roundup(x):
    if x < 60:
        return 50
    return int(math.ceil(x / 100.0)) * 100

def gf2_sum(vec1,vec2):
    vec3= np.zeros(len(vec1))
    for i in range(len(vec1)):
        if vec1[i]!=vec2[i]:
            vec3[i] = 1
    return vec3


def gq_sum(vec1,vec2):
    vec3= np.zeros(len(vec1))
    for i in range(len(vec1)):
        if vec1[i]==vec2[i]:
            vec3[i] = 0
        elif (vec1[i]==0 and vec2[i]==1) or (vec2[i]==0 and vec1[i]==1) :
            vec3[i] =1
        elif (vec1[i]==0 and vec2[i]==-1) or (vec2[i]==0 and vec1[i]==-1):
            vec3[i]=-1
        else:
            vec3[i] = 0
    return vec3

def gq_multi(vec1,vec2):
    vec3= np.zeros(len(vec1))
    for i in range(len(vec1)):
        if (vec1[i]==0 )or (vec2[i] == 0):
            vec3[i] = 0
        
        elif (vec1[i]==-1 and vec2[i]==-1)  :
            vec3[i] =1
        elif (vec1[i]==1 and vec2[i]==1)  :
            vec3[i] =1
        elif (vec1[i]==1 and vec2[i]==-1) or (vec2[i]==1 and vec1[i]==-1):
            vec3[i]=-1
        else:
            vec3[i] = 0
    return vec3


def linear_independence(vectors):
    N = len(vectors)
    
    for v in vectors:
        independent = True
        vvectors = [vv for vv in vectors if vv != v]
        for n in range(2,N):
            for SV in itertools.combinations(vvectors,n):
                pass
            
    


"""
def dag_cycle_basis_depina(graph,return_vector_edgeid = True):
    G = graph.copy()
    print(G.number_of_edges()-G.number_of_nodes()+1)
    T = nx.dfs_tree(G,[n for n in G if G.in_degree(n)==0][0])
    back_edges = []
    for e in G.edges():
        if e not in T.edges():
            back_edges.append(e)
    N= G.number_of_edges()-G.number_of_nodes()+1
    #print(N)
    edge_id = {}
    count = 0
    for e in back_edges:
        edge_id[e] = count
        count+=1
    for e in T.edges():
        edge_id[e] = count
        count+=1
    witnesses= []#cycle_basis_vector(G_root,list(e),edge_id) for e in back_edges]
    for be in back_edges:
        wa = np.zeros(len(edge_id))  
        wa[edge_id[be]] = 1
        witnesses.append(wa)
    id_edge = {val:key for key,val in edge_id.items()}        
    
    pnodes = list(G.nodes())
    nnodes = list(G.nodes())
    nnodes = [-n for n in nnodes]
    CB = []
    id_edge = {val:key for key,val in edge_id.items()}
    used_back_edges = []
    for i in  range(N):
        #print(i)
        S = witnesses[i]
        G_i = nx.Graph()
        G_i.add_nodes_from(nnodes+pnodes)
        for (u,v) in G.edges():
            if S[edge_id[(u,v)]] !=1:
                G_i.add_edge(u,v)
                G_i.add_edge(-u,-v)
            else:
                G_i.add_edge(u,-v)
                G_i.add_edge(-u,v)
        current_cycles = []
        for n in pnodes:
            if nx.has_path(G_i,-n,n):
            
                current_cycles.append(nx.shortest_path(G_i,-n,n)) 
            #except:
            #    pass
            #print("dir sp size",len(current_cycles))
        for c in sorted(current_cycles,key=len):
            c = [abs(cc) for cc in c]#[:-1]
            #print(c)
            vec= np.zeros(G.number_of_edges())
            for m in range(len(c)-1):
                n1,n2 = c[m],c[m+1]
                try:
                    k=edge_id[(n1,n2)]
                except:
                    k = edge_id[(n2,n1)]
                vec[k] = 1
            if np.dot(vec,S.T) % 2 ==1 :#and has_new_back_edge(all_back_edges,vec,edge_id,used_back_edges,return_edge =False)==True:#and vec[i]!=0:#and sum(vec)>3:# and len(np.where(CB==vec))==0:
                ns = []
                for index in range(len(vec)):
                    if vec[index] == 1:
                        ns.append(id_edge[index][0])
                        ns.append(id_edge[index][1])
                CB.append(vec)
                break
        for j in range(i+1,N):
            S_j = witnesses[j]
            if np.dot(vec,S_j.T) % 2 ==1:
                S_j = gf2_sum(S,S_j)
                witnesses[j] = S_j
    
    cycle_basis = []
    id_edge = {val:key for key,val in edge_id.items()}
    for cv in CB:
        nodeset = []
        for i in range(len(cv)):
            if cv[i] ==1:
                u,v = id_edge[i]
                nodeset += [u,v]
        cycle_basis.append(list(set(nodeset)))
    if return_vector_edgeid==True:  
        return cycle_basis,edge_id,CB
    else:
        return cycle_basis
"""

def dag_cycle_basis_depina(graph,return_vector_edgeid = True):
    print("other depina")
    G = graph.copy()
    T = nx.dfs_tree(G,[n for n in G if G.in_degree(n)==0][0])
    back_edges = []
    for e in G.edges():
        if e not in T.edges():
            back_edges.append(e)
    N= G.number_of_edges()-G.number_of_nodes()+1
    #print(N)
    edge_id = {}
    count = 0
    for e in back_edges:
        edge_id[e] = count
        count+=1
    for e in T.edges():
        edge_id[e] = count
        count+=1
    witnesses= []#cycle_basis_vector(G_root,list(e),edge_id) for e in back_edges]
    for be in back_edges:
        wa = np.zeros(len(edge_id))  
        wa[edge_id[be]] = 1
        witnesses.append(wa)
    id_edge = {val:key for key,val in edge_id.items()}        
    
    pnodes = list(G.nodes())
    nnodes = list(G.nodes())
    nnodes = [-n for n in nnodes]
    CB = []
    id_edge = {val:key for key,val in edge_id.items()}
    used_back_edges = []
    for i in  range(N):
        #print(i)
        S = witnesses[i]
        #print("current_back edge",id_edge[i])
        cbe = id_edge[i]
        #print(S)
        #find the shortest cycle in G_root such that <S,C_i>=1
        G_i = nx.Graph()
        G_i.add_nodes_from(nnodes+pnodes)
        
        for (u,v) in G.edges():
            if S[edge_id[(u,v)]] !=1:
                G_i.add_edge(u,v)
                G_i.add_edge(-u,-v)
            else:
                G_i.add_edge(u,-v)
                G_i.add_edge(-u,v)
        current_cycles = []
        for n in nnodes:
            try:
                current_cycles.append(nx.shortest_path(G_i,-n,n)) 
            except:
                pass
        #print("dir sp size",len(current_cycles))
        for c in sorted(current_cycles,key=len):
            c = [abs(cc) for cc in c]#[:-1]
            #print(c)
            vec= np.zeros(G.number_of_edges())
            for m in range(len(c)-1):
                n1,n2 = c[m],c[m+1]
                try:
                    k=edge_id[(n1,n2)]
                except:
                    k = edge_id[(n2,n1)]
                vec[k] = 1
            if np.dot(vec,S.T) % 2 ==1 :#and has_new_back_edge(all_back_edges,vec,edge_id,used_back_edges,return_edge =False)==True:#and vec[i]!=0:#and sum(vec)>3:# and len(np.where(CB==vec))==0:
                ns = []
                for index in range(len(vec)):
                    if vec[index] == 1:
                        ns.append(id_edge[index][0])
                        ns.append(id_edge[index][1])
                CB.append(vec)
                break
        for j in range(i+1,N):
            S_j = witnesses[j]
            if np.dot(vec,S_j.T) % 2 ==1:
                S_j = gf2_sum(S,S_j)
                witnesses[j] = S_j
    
    cycle_basis = []
    id_edge = {val:key for key,val in edge_id.items()}
    for cv in CB:
        nodeset = []
        for i in range(len(cv)):
            if cv[i] ==1:
                u,v = id_edge[i]
                nodeset += [u,v]
        cycle_basis.append(list(set(nodeset)))
    if return_vector_edgeid==True:  
        return cycle_basis,edge_id,CB
    else:
        return cycle_basis


def _random_subset(seq,m):
    """ Return m non-unique elements from seq.
    This differs from random.sample which can return repeated
    elements if seq holds repeated elements.
    """
    targets=random.sample(seq,m)
    return targets





def dag_cycle_basis_horton(graph,return_vector_edgeid = True):
    G = graph.copy()
    T = nx.dfs_tree(G,[n for n in G if G.in_degree(n)==0][0])
    #print("tree from ",[n for n in G if G.in_degree(n)==0][0])
    back_edges = []
    for e in G.edges():
        if e not in T.edges():
            back_edges.append(e)
    N= G.number_of_edges()-G.number_of_nodes()+1
    #print(N,N)
    edge_id = {}
    count = 0
    for e in back_edges:
        edge_id[e] = count
        count+=1
    for e in T.edges():
        edge_id[e] = count
        count+=1
    witnesses= []
    for be in back_edges:
        wa = np.zeros(len(edge_id))  
        wa[edge_id[be]] = 1
        witnesses.append(wa)
    id_edge = {val:key for key,val in edge_id.items()}     
    
    horton_set = []
    for v in G:
        for u in nx.descendants(G,v):
            for w in G.successors(u):
                SP_vu =set(nx.shortest_path(G,v,u))
                SP_vw =set(nx.shortest_path(G,v,w))
                if SP_vu.intersection(SP_vw) == set([v]):
                    C = list(SP_vu)+list(SP_vw)+[u,v]
                    
                    horton_set.append(set(C))
        
    horton_set= sorted(horton_set,key=len)
    CB = []
    #print(horton_set)
    for i in  range(N):
        S = witnesses[i]
        for c in horton_set:
            vec= np.zeros(G.number_of_edges())
            gc = get_subgraph(G,c)
            for (n1,n2) in gc.edges():
                k = edge_id[(n1,n2)]
                vec[k] = 1
            #print(np.dot(S,vec.T))
            if np.dot(S,vec.T) % 2 ==1:# and np.linalg.matrix_rank(np.row_stack(CB+[vec]))==len(CB)+1:       
                CB.append(vec)
                #print(i)
                break
        
        for j in range(i+1,N):
            S_j = witnesses[j]
            if np.dot(S_j,vec.T)% 2==1:
                S_j = gf2_sum(S,S_j)
                witnesses[j] = S_j
    
    cycle_basis = []
    id_edge = {val:key for key,val in edge_id.items()}
    for cv in CB:
        nodeset = []
        for i in range(len(cv)):
            if cv[i] ==1:
                u,v = id_edge[i]
                nodeset += [u,v]
        cycle_basis.append(list(set(nodeset)))
    A  = np.row_stack(CB)
    #r = np.linalg.matrix_rank(A)
    #print(r,len(CB))
    #if G.number_of_edges()-G.number_of_nodes()+1!=scipy.linalg.null_space((A.dot(A.T))).shape[0]:
    #    return None
    if return_vector_edgeid==True:  
        return cycle_basis,edge_id,CB
    else:
        return cycle_basis
    
def print_outersection(lol1,lol2):
    lol1 = set(frozenset(l) for l in lol1)
    lol2 = set(frozenset(l) for l in lol2)
    return lol1.symmetric_difference(lol2)
    

def directed_cycle_vector(G,c,edge_id):
    
    n0 = c[0]
    edges_undirected = []
    for n0,n1 in itertools.combinations(c,2):
        if (n0,n1) in G.edges() or (n1,n0) in G.edges():
            edges_undirected.append((n0,n1))
    n0,n1 = edges_undirected[0]
    n_origin = n0
    m0,m1 = edges_undirected[0]
    current_direction = 1
    edges = [(n0,n1)]
    edges_undirected = [e for e in edges_undirected if e != (m0,m1)]
    edge_direction = [1]
    while m1 !=n_origin:
        print(n0,n1)
        e = [(u,v) for (u,v) in edges_undirected if n1 ==u]
        if e == []:
            e = [(u,v) for (u,v) in edges_undirected if n1 ==v]
            current_direction = current_direction*-1
            m1,m0 = e[0]
        else:
            m0,m1 =e[0]
        edge_direction.append(current_direction)
        edges.append((m0,m1))
        n0,n1 = m0,m1
        edges_undirected = [e for e in edges_undirected if e != (m0,m1)]
 

def print_cycle_statistics(graph,cycles=None, return_eigvals = False,ctype = "DAG"):
    G = graph.copy()
    E_0 = G.number_of_edges()
    tr(G)
    
    ###
    # Compute MCB
    ###
    if cycles == None and ctype =="DAG":
        print("No cycles found, they will be computed")
        #E = G.number_of_edges()
        Gud = G.copy()
        Gud = Gud.to_undirected()
        C = minimum_cycle_basis(Gud)
        print(len(C), G.number_of_edges()-G.number_of_nodes()+1)
    elif cycles != None:
        C = cycles
    elif ctype =="horton":
        C = dag_cycle_basis_horton(G,True)[0]
    elif ctype =="depina":
        print("type=depina")
        C = dag_cycle_basis_depina(G,True)[0]
        
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
    result["Number of edges"] = E_0
    result["Balance"] = np.mean(balance)
    result["Balance_norm"] = np.mean(balance_norm)
    result["Number_diamonds"] = len([x for x in cycle_type if x == "D"])
    result["Number_mixers"] = len([x for x in cycle_type if x == "M"])
    #result["Largest height antichain"] = get_largest_height_antichain(G)
    result["Largest eigenvalue M_C"] = np.real(max(eigvals))
    result["Largest eigenvalue A_C"] = np.real(max(np.linalg.eigvals(A_C)))
    #result["Spectral gap of L_C"] = np.real(lambda_min_l)
    result["Mean cycle height"] = np.mean([x[0] for x in cycle_height])#/ nx.dag_longest_path_length(G)
    result["Std cycle height"] = np.std([x[0] for x in cycle_height])#/ nx.dag_longest_path_length(G)    
    result["Mean stretch"] =np.mean(var_cycle_height)#max([x[0]/ nx.dag_longest_path_length(G) for x in cycle_height])-min([x[0]/ nx.dag_longest_path_length(G) for x in cycle_height])
    result["eigenvalue_ratio"] = np.real(lambda_max_a)/np.real(max(eigvals))
    result["Largest M_C by corr size"] = np.real(max(eigvals))/[len(c) for c in C][list(eigvals).index(max(eigvals))]
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
    

def complete_dag(n,c):
    """
    Create a DAG with n nodes in which an edge exists between i and j if i<j. 
    Return a list of nodes in which each node i appears k_out(i) number of times.
    """
    
    graph = nx.DiGraph()
    list_repeated_nodes  = []
    for i in range(n):
        graph.add_node(i)
        list_repeated_nodes.extend([i]*c)
        for j in range(i+1,n):
            graph.add_edge(i,j)
            list_repeated_nodes.append(i)
                       
    
    return list_repeated_nodes, graph

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


def add_dag_edges(graph,p_target):
    G = graph.copy()
    N = G.number_of_nodes()
    E = int(N*(N-1)*p_target/2)
    edges = set(G.edges())
    new_edges = []
    nodes = list(G.nodes())
    no_edges = G.number_of_edges()
    while no_edges<E:
        n1,n2 = random.sample(nodes,2)
        if n1>n2 and (n2,n1) not in edges:
            new_edges.append((n2,n1))
            no_edges +=1
        elif n1<n2 and (n1,n2) not in edges:
            new_edges.append((n1,n2))
            no_edges+=1
    G.add_edges_from(new_edges)
    return G

def lattice_dag_2D(N,L):
    
    #add nodes
    Gud = nx.grid_2d_graph(N,L)
    G = nx.DiGraph()
    
    for e1,e2 in Gud.edges():
        if e1[1]<=e2[1]:
            G.add_edge(e1,e2)
    return G


def multi_tr(DAG):
    edges = list(DAG.edges())
    #########################
    for edge in edges:
        # check edge is necessary for causal structure
        [u, v] = edge
        if (u,v) in DAG.edges():
            DAG.remove_edge(u, v)
            if not nx.has_path(DAG, u, v):
                DAG.add_edge(u, v)
            else:
                sp = nx.shortest_path(DAG, u, v)
                sp = [(sp[i],sp[i+1]) for i in range(len(sp)-1)]
                
                
                if len(sp) == 2:
                #check if the third shortest path is also of lenght two
                #that case delete the entire clique
                    for x,z in sp:
                        DAG.remove_edge(x,z)
                    if nx.has_path(DAG,u,v):
                        sp2 = nx.shortest_path(DAG, u, v)
                        sp2 = [(sp2[i],sp2[i+1]) for i in range(len(sp2)-1)]
                        if len(sp2)==2:
                              #create a supernode
                              nodes = set()
                              for (n1,n2) in sp2:
                                  nodes.add(n1)
                                  nodes.add(n2)
                              for (n1,n2) in sp:
                                  nodes.add(n1)
                                  nodes.add(n2)
                              in_stubs,out_stubs = set(),set()
                              for n in nodes:
                                  for s in nx.descendants(DAG,n):
                                      if s not in nodes:
                                          in_stubs.add(s)
                                  for s in DAG.successors(n):
                                      if s not in nodes:
                                          out_stubs.add(s)
                              DAG.remove_nodes_from(nodes)
                              nnew=DAG.number_of_nodes()+1
                              DAG.add_node(nnew)
                              for s in in_stubs:
                                  DAG.add_edge(s,nnew)
                              for s in out_stubs:
                                  DAG.add_edge(nnew,s)
                                      
                        if len(sp2)>2:
                            #return the shortest paths, but not the original edge u,v
                            for (n1,n2) in sp:
                                DAG.add_edge(n1,n2)
                            for (n1,n2) in sp2:
                                DAG.add_edge(n1,n2)
                                
                                
                    else:
                        for x,z in sp:
                            #because there were two holes
                            DAG.add_edge(x,z)
                        
                            
                            
                if len(sp) >2:
                    #xto preserve a hole
                    DAG.add_edge(u,v)
    return DAG
                    
                
            



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
        
            
        
def cycle_basis_vector_directed(graph_original,cycle,edge_id):
    """
    Find cycle basis vector for a cycle. A cycle consists of forward edges and 
    backward edges. Correspondingly, the forward edges are given "+1" entry in the edge
    vector, backward edges - "-1". Edges that don't feature in the cycle are given "0".
    
    Parameters
    ----------
    graph_original - networkx.DiGraph directed acyclic graph
    cycle - list of nodes in the cycle
    edge_id - dictionary of edges (keys) and their corresponding ids (values)
    
    Returns 
    -------
    vector - a cycle basis vector for the given cycle
    """
    
    gc = get_subgraph(graph_original,cycle)
    
    curr_e = list(gc.edges())[0]
    direction = {curr_e:1}
    visited = {curr_e}
    curr_dir = 1
    #print(c)
    while len(visited) <gc.number_of_edges():
        #print(len(visited))
        count =0
        #print(direction)
        #print(visited)
        for e in gc.edges():
            count +=1
            if e not in visited and (e[1],e[0]) not in visited and e[0]==curr_e[1]:
                #print("i go from",curr_e,"to",e)
                visited.add(e)
                direction[e] = 1*curr_dir
                curr_dir = 1*curr_dir
                curr_e = e
                break
            elif e not in visited and (e[1],e[0]) not in visited and e[1]==curr_e[1]:
                #print("i go from",curr_e,"to",e)
                visited.add(e)
                direction[e] = -1*curr_dir
                curr_dir = -1*curr_dir
                gc = gc.reverse()
                curr_e = (e[1],e[0])
                break
        if count == gc.number_of_edges():
            curr_dir = -1*curr_dir
            gc= gc.reverse()
    vector = [0 for e in range(len(graph_original.edges()))]
    for e,val in direction.items():
        try:
            vector[edge_id[e]] =val
        except:
            e = (e[1],e[0])
            vector[edge_id[e]]=val
    return vector#,direction



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

