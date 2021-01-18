# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 11:08:56 2019

@author: Vaiva
"""



from cycle_utilities import *
import pandas as pd


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

    
    
def random_dag_iterative(N,P):
    
    G = nx.DiGraph()
    i = 1
    for i in range(1,N+1):
        G.add_node(i)
        for j in range(1,i):
            p = random.random()
            if p <= P:
                G.add_edge(j,i)
                
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



            
        
def cycle_basis_vector(graph_original,cycle,edge_id):
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

        

task="plot"
if task=="__main__":
    
    ntype = "random"
    N =500
    ctype = "DAG"
    if ntype =="random":
        pvalues = [0.6,0.7,0.8,0.9]
    elif ntype =="price":
        pvalues = [2,3,4,5,6,7,8]
    elif ntype =="lattice" or ntype =="russian_doll":
        pvalues = [5, 10,15,20,25,30]#[int(N**0.5)]
    filename = ntype+"_dags_cycle_data_"+ctype +"_no_nodes_eq_"+str(N)+"_1005.txt"
    for p in pvalues:
        print(p)
        r = 0
        while r <20:
            if ntype =="random":
                G = random_dag(N,p)
            elif ntype =="price":
                G = price_dag(N,p,1,N+1)
            elif ntype =="lattice":
                G = lattice_dag_2D(p,p)
            elif ntype =="russian_doll":
                G = russian_doll_dag(p*p)
            #tr(G)
            if G :
                r+=1
                print(r)
                if nx.is_weakly_connected(G) ==True:    
                    try:
                        E =G.number_of_edges()
                        tr(G)
                        #C = dag_cycle_basis_depina(G,True)[0]
                        #print("number of cycles in depina",len(C))
                        res = print_cycle_statistics(G,ctype = "DAG")
                        #name = randomString()
                        with open(filename,"a+") as f:
                            strname = str(p)+"\t"+str(E)+"\t"
                            for key,val in res.items():
                                strname+=str(val)+"\t"
                            f.write(strname+"\n") 
                    except:
                        pass
                    #print(M_C.shape)
                    
                    

if task == "quasi_unicity":
    N=50
    G = random_dag(N,0.1)
    tr(G)
    edges  = list(G.edges())
    df = pd.DataFrame()# print_cycle_statistics(G,cycles = C),orient ="index").T
    x_ran ,x_err_ran ,x_name,x_obs= [], [],[], []
    for l in range(10):
        Gud = nx.to_undirected(G)
        #C =minimum_cycle_basis(Gud)# dag_cycle_basis_depina(G,True)[0]
        C =dag_cycle_basis_horton(G)[0]
        print(G.number_of_edges())
        res = print_cycle_statistics(G,cycles = C)
        res = {key:[val] for key,val in res.items()}
        df =df.append(pd.DataFrame(res))
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
        print(len(C))
    for c in df.columns:
        x_ran.append(df[c].mean())
        x_err_ran.append(df[c].std())
        x_name.append("Random p=0.1")
        x_obs.append(c)
    G = random_dag(N,0.8)
    tr(G)
    edges  = list(G.edges())
    df = pd.DataFrame()# print_cycle_statistics(G,cycles = C),orient ="index").T
    for l in range(10):
        
        Gud = nx.to_undirected(G)
        #C = minimum_cycle_basis(Gud)#[0]
        
        C =dag_cycle_basis_horton(G)[0]
        res = print_cycle_statistics(G,cycles = C)
        res = {key:[val] for key,val in res.items()}
        df =df.append(pd.DataFrame(res))
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
    for c in df.columns:
        x_ran.append(df[c].mean())
        x_err_ran.append(df[c].std())
        x_name.append("Random p=0.8")
        x_obs.append(c)
    G = price_dag(N,3,1,N+1)
    
    tr(G)
    edges  = list(G.edges())
    df = pd.DataFrame()# print_cycle_statistics(G,cycles = C),orient ="index").T
    for l in range(10):
        Gud = nx.to_undirected(G)
        #C = minimum_cycle_basis(Gud)#[0]
        
        C =dag_cycle_basis_horton(G)[0]
        res = print_cycle_statistics(G,cycles = C)
        res = {key:[val] for key,val in res.items()}
        df =df.append(pd.DataFrame(res))
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
    for c in df.columns:
        x_ran.append(df[c].mean())
        x_err_ran.append(df[c].std())
        x_name.append("Price m=3")
        x_obs.append(c)
    
    G = price_dag(N,5,1,N+1)
    tr(G)
    edges  = list(G.edges())
    df = pd.DataFrame()# print_cycle_statistics(G,cycles = C),orient ="index").T
    for l in range(10):
        Gud = nx.to_undirected(G)
        #C = minimum_cycle_basis(Gud)#[0][0]
        C =dag_cycle_basis_horton(G)
        res = print_cycle_statistics(G,cycles = C)
        res = {key:[val] for key,val in res.items()}
        df =df.append(pd.DataFrame(res))
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
    for c in df.columns:
        x_ran.append(df[c].mean())
        x_err_ran.append(df[c].std())
        x_name.append("Price m=5")
        x_obs.append(c)
        
    df = pd.DataFrame([x_ran,x_err_ran,x_name,x_obs]).T
    df.columns =["mean","std","network","observable"]
    groups = df.groupby("observable")
    for key,group in groups:
        print(key)                
        plt.figure(figsize =(4,5))
        (_, caps, _) = plt.errorbar(x=group["network"],y=group["mean"],yerr=group["std"],fmt = "o",color="b" ,marker = "o",capsize= 10  )        
        for cap in caps:
            cap.set_color('b')
            cap.set_markeredgewidth(1)
        plt.tick_params(axis= "both",labelsize=17)
        plt.tick_params(axis="x", rotation= 90)
        plt.title(key,fontsize = 15)
        plt.tight_layout()
        #plt.savefig(key+"_quasiunicity_N_eq_100_iter_eq_10.pdf")
                    
                    
if task == "plot_eigenspectrum":
    
    ntype = "random"
    N =500
    ctype = "depina"
    
    if ntype =="random":
        pvalues = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    elif ntype =="price":
        pvalues = [2,3,4,5,6,7,8]
    elif ntype =="lattice" or ntype =="russian_doll":
        pvalues = [5,10,15,20,25,30]#[int(N**0.5)]
    filename = ntype+"_dags_cycle_data_"+ctype +"_no_nodes_eq_"+str(N)+"_0819.txt"
    plt.figure(figsize=  (6,5))
    symbol = ["o","^",">","<","*","v","+","x","d","1"]
    egvals_at_p = {}
    for p in pvalues:
        print(p)
        r = 0
        sumegval = []
        while r <1:
            if ntype =="random":
                    G = random_dag(N,p)
            elif ntype =="price":
                G = price_dag(N,p,1,N+1)
            elif ntype =="lattice":
                G = lattice_dag_2D(p,p)
            elif ntype =="russian_doll":
                G = russian_doll_dag(p*p)
            tr(G)
            try:
                Gud = nx.to_undirected(G)
                C = minimum_cycle_basis(Gud)
                
                egvals = eigenspectrum(G,C)
                sumegval.append(sum(egvals)/len(egvals))
            except:
                pass
        egvals_at_p[p] = sumegval
        plt.scatter([i/len(C) for i in range(len(C))],egvals,marker = symbol[pvalues.index(p)],label="p="+str(p))
             
    #plt.scatter([a for a in egvals_at_p.keys()],[np.mean(x) for x in egvals_at_p.values()],[np.std(x) for x in egvals_at_p.values()],capsize = 5,fmt=".")
    plt.xlabel("$1-i/d$",fontsize = 17)
    plt.ylabel("$\\lambda^C_i$",fontsize = 17)
    plt.tick_params(axis= "both",labelsize=17)
    plt.legend(fontsize =14)
    #plt.ylim(4,7)
    plt.tight_layout()
    plt.savefig(ntype+"_"+str(N)+"_eigenspectrum.pdf")
     
               
if task=="plot":
    
    ##########
    #FOR PRICE
    ##########
    
    data = pd.read_csv("price_dags_cycle_data_depina_no_nodes_eq_500_0826.txt",sep = "\t")
    #data = pd.read_csv("random_dags_cycle_data_DAG_no_nodes_eq_500_1005.txt",sep = "\t")
    data["Density"] = data['Number of TR edges']/(250*499)
    
    xlabel = "m"
    filestart = "price"
    mean_list = []
    for m in set(data["m"]):
        mean_list.append(list(data[data["m"]==m].mean()))
    
    mean_df = pd.DataFrame(mean_list)
    mean_df.columns = data.columns[:]
    
    std_list = []
    for m in set(data["m"]):
        std_list.append([m]+list(data[data["m"]==m].std())[1:])    
    mean_df.index = mean_df.m
    std_df = pd.DataFrame(std_list)
    std_df.columns = data.columns[:]
    std_df.index = std_df.m
    
    
    """
    
    fig, ax = plt.subplots(1,1,figsize = (5,4))
    ax.errorbar(mean_df["m"],mean_df["Density"],std_df["Density"],fmt= "o",capsize=5,label = "Random",color="red")
    data = pd.read_csv("price_dags_cycle_data_depina_no_nodes_eq_500_0826.txt",sep = "\t")
    data["Density"] = data['Number of TR edges']/(250*499)
    xlabel = "m"
    filestart = "random"
    mean_list = []
    for m in set(data["m"]):
        mean_list.append(list(data[data["m"]==m].mean()))
    
    mean_df = pd.DataFrame(mean_list)
    mean_df.columns = data.columns[:]
    
    std_list = []
    for m in set(data["m"]):
        std_list.append(list(data[data["m"]==m].std()))    
    mean_df.index = mean_df.m
    std_df = pd.DataFrame(std_list)
    std_df.columns = data.columns[:]
    std_df.index = std_df.m
    
    ax2 = ax.twiny()  # instantiate a second axes that shares the same x-axis
    ax2.set_xlabel(xlabel,fontsize = 18,c="g")
    ax2.tick_params(axis= "y",labelsize=18,labelcolor="g")
    ax2.errorbar(mean_df["m"],mean_df["Density"],std_df["Density"],fmt= "o",capsize=5,label = "Price",color="green")
    ax.set_ylabel("$2E_{\\mathrm{TR}}/N(N-1)$",fontsize = 18)
    ax.set_xlabel("p",fontsize = 18,c="r")
    ax2.tick_params(axis= "y",labelsize=18)
    ax.tick_params(axis ="x", labelsize=18)
    ax.tick_params(axis= "both",labelsize=18)
    ax2.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig("density_price_random_N_eq_500.pdf")
    
    
    """
    
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df["m"],mean_df["Number of TR edges"],std_df["Number of TR edges"],fmt= "o",capsize=5,)#label = "Mixers")
    #plt.legend(fontsize = 14)
    plt.ylabel("$E_{\\mathrm{TR}}$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"_E_TR_N_eq_500.pdf")
    
    
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df["m"],mean_df["Mean cycle height"],std_df["Mean cycle height"],fmt= "o",capsize=5,)#label = "Mixers")
    #plt.legend(fontsize = 14)
    plt.ylabel("$\\langle h \\rangle$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"_h_N_eq_500.pdf")
    
 
 
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df["m"],mean_df["Mean edge participation"],std_df["Mean edge participation"],fmt= "o",capsize=5,)#label = "Mixers")
    #plt.legend(fontsize = 14)
    plt.ylabel("$E_P$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"_edge_participation_N_eq_500.pdf")
 
    
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df["m"],mean_df["Mean cycle size"],std_df["Mean cycle size"],fmt= "o",capsize=5,)#label = "Mixers")
    #plt.legend(fontsize = 14)
    plt.ylabel("$\\langle C \\rangle$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"_mean_cycle_size_N_eq_500.pdf")
    
        
 
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df["m"],mean_df['Number_diamonds'],std_df["Number_diamonds"],fmt= "o",capsize=5,label = "Diamonds")

    plt.errorbar(mean_df["m"],mean_df['Number_mixers'],std_df["Number_mixers"],fmt= "o",capsize=5,label = "Mixers")
    plt.legend(fontsize = 18)
    plt.ylabel("Number of cycles",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"_number_mixers_diamonds_N_500.pdf")
    
    fig, ax = plt.subplots(1,1,figsize = (5,4))
    ax.errorbar(mean_df["m"],mean_df['Largest eigenvalue M_C'],std_df['Largest eigenvalue M_C'],fmt= "o",color = "b",capsize=5)#label = "Mixers")
    #plt.legend(fontsize = 14)
    ax.set_ylabel("$\\lambda^C_{\\mathrm{max}}$",fontsize = 16, color = "b")
    ax.set_xlabel(xlabel,fontsize = 18)
    ax.tick_params(axis= "y",labelsize=18)#, color = "b")
    
    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.errorbar(mean_df["m"]+0.2,mean_df['Largest M_C by corr size'],std_df['Largest M_C by corr size'],fmt= "o",color="g",capsize=5)#label = "Mixers")
    #plt.legend(fontsize = 14)
    ax2.set_ylabel("$\\lambda^C_{\\mathrm{max}}/C$",fontsize = 16, color = "g")
    ax2.set_xlabel(xlabel,fontsize = 18)
    ax2.tick_params(axis= "y",labelsize=18)
    ax.tick_params(axis ="x", labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"_eigvals_N_eq_500.pdf")
    
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df['m'],mean_df["Number cycle connected components"],std_df["Number cycle connected components"], fmt= "o",capsize=5,)#label = "Mixers")
    #plt.ylim(0.7,1.4)
    #plt.legend(fontsize = 14)
    plt.ylabel("null$L^C$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"_null_L_N_eq_500.pdf")
 
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df['m'],mean_df["Std cycle size"],std_df["Std cycle size"],fmt= "o",capsize=5,)#label = "Mixers")
    #plt.ylim(0.9,1.1)
    #plt.legend(fontsize = 14)
    plt.ylabel("$\\sigma(C)$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"std_cycle_size_N_eq_500.pdf")
    
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df['m'],mean_df["Mean cycle size"],std_df["Mean cycle size"],fmt= "o",capsize=5,)#label = "Mixers")
    #plt.ylim(0.9,1.1)
    #plt.legend(fontsize = 14)
    plt.ylabel("$\\langle C \\rangle$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"mean_cycle_size_N_eq_500.pdf")
    
    
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df['m'],mean_df["Largest cycle size"],std_df["Largest cycle size"],fmt= "o",capsize=5,)#label = "Mixers")
    plt.ylim(5,10)
    #plt.legend(fontsize = 14)
    plt.ylabel("$ C_{\\mathrm{max}}$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"max_cycle_size_N_eq_500.pdf")

    
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df['m'],mean_df["Mean stretch"],std_df["Mean stretch"],fmt= "o",capsize=5,)#label = "Mixers")
    #plt.ylim(0.9,1.1)
    #plt.legend(fontsize = 14)
    plt.ylabel("$\\langle s \\rangle$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"mean_stretch_size_N_eq_500.pdf")
    
    plt.figure(figsize =(5,4))
    plt.errorbar(mean_df['m'],mean_df["Balance"],std_df["Balance"],fmt= "o",capsize=5,)#label = "Mixers")
    plt.ylim(0.0,0.2)
    #plt.legend(fontsize = 14)
    plt.ylabel("$\\langle b \\rangle$",fontsize = 18)
    plt.xlabel(xlabel,fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig(filestart+"balance_size_N_eq_500.pdf")
    
    
    
    

    plt.figure(figsize =(5,4))
    df = pd.read_csv("lattice_dags_cycle_data_horton_310520.txt",sep ="\t")
    df.index = df["number of nodes "]
    df["Largest eigenvalue M_C"].plot(style ="-o",label  = "Lattice")
    df = pd.read_csv("russian_doll_dags_cycle_data_horton_310520.txt",sep ="\t")
    df.index = df["number of nodes "]
    df["Largest eigenvalue M_C"].plot(style ="-o",label  = "Russian doll")
    plt.ylabel("$\\lambda^C_{\\mathrm{max}}$",fontsize = 18)
    plt.xlabel("Number of nodes",fontsize = 18)
    plt.legend(fontsize = 18)
    plt.tick_params(axis= "both",labelsize=18)
    plt.tight_layout()
    plt.savefig("eigenvalues_lattice_russian_doll.pdf")
    
    
    
    