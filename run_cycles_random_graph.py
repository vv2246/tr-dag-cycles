# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 11:08:56 2019

@author: Vaiva
"""



from cycle_utilities import *

    
        
def calc_mcb(ntype, N):
    """
    

    Parameters
    ----------
    ntype : type of DAG: random, Price, lattice, or russian doll.
    N : number of nodes.

    Returns
    -------
    None. Saves results into a text file.

    """
    N =50
    if ntype =="random":
        pvalues = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    elif ntype =="price":
        pvalues = [2,3,4,5,6,7,8]
    elif ntype =="lattice" or ntype =="russian_doll":
        pvalues = [5, 10]#,15,20,25,30]#[int(N**0.5)]
    filename = "dags_cycle_data_{}_no_nodes_eq_{}.txt".format(ntype,N)
    for p in pvalues:
        print(p)
        r = 0
        while r <10:
            if ntype =="random":
                G = random_dag(N,p)
            elif ntype =="price":
                G = price_dag(N,p,1,N+1)
            elif ntype =="lattice":
                G = lattice_dag_2D(p,p)
            elif ntype =="russian_doll":
                G = russian_doll_dag(p*p)
            if G :
                if nx.is_weakly_connected(G) ==True:    
                    try:
                        E =G.number_of_edges()
                        res = print_cycle_statistics(G)
                        with open(filename,"a+") as f:
                            strname = str(p)+"\t"+str(E)+"\t"
                            for key,val in res.items():
                                strname+=str(val)+"\t"
                            f.write(strname+"\n") 
                        r+=1
                        print(r)
                    except:
                        pass
                    
                    
                    
                    
                    
def plot(filename,fig_save = False):
    
    plt.rcParams.update({'font.size': 30})
    data = pd.read_csv(filename,sep="\t",header= None)
    data.columns = ["p","E","E_p","E_p_std","S_max","S_mean","S_std","C","null_L","E_TR","L_max","N","b","D","M","lambda_max","h_mean","h_std","s_mean","s_std","lambda_max_byC","Nan"]
    if "price" in filename:
        xlabel = "m"
        data.index = data.p
        x = data.groupby(data.index)["p"].mean()
    elif "random" in filename:
        data.index = data.p
        x = data.groupby(data.index)["p"].mean()
        xlabel = "p"
    elif "lattice" in filename or "russian" in filename:
        xlabel = "N"
        
        data.index = data.N
        x = data.groupby(data.index)["N"].mean()
            
    figure_filename = filename.rstrip(".txt")+"_{}.pdf"
    
    for c in data.columns:
        y =  data.groupby(data.index)[c].mean()
        yerr = data.groupby(data.index)[c].std()
        if "price" in filename or "random" in filename:
            plt.errorbar(x,y,yerr,fmt = "o",capsize= 5)
        else:
            plt.errorbar(x,y,fmt="o",capsize= 0)
        plt.xlabel(xlabel)
        plt.ylabel(c)
        if fig_save:
            plt.tight_layout()
            plt.savefig(figure_filename.format(c))
        
        plt.show()
        
        
def test_quasiunicity(N,fig_save = False):
    
    N=50
    
    G = random_dag(N,0.1)
    tr(G)
    edges  = list(G.edges())
    df = pd.DataFrame()
    df_list = []
    for l in range(10):
        Gud = nx.to_undirected(G)
        C =minimum_cycle_basis(Gud)
        res = print_cycle_statistics(G,cycles = C)
        res = ["random_0_1"] + list(res.values())
        df_list.append(res)
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
        
    G = random_dag(N,0.8)
    tr(G)
    edges  = list(G.edges())
    for l in range(10):
        Gud = nx.to_undirected(G)
        C =minimum_cycle_basis(Gud)
        res = print_cycle_statistics(G,cycles = C)
        res = ["random_0_8"] + list(res.values())
        df_list.append(res)
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
        
        
    G = price_dag(N,3,1,N+1)
    tr(G)
    edges  = list(G.edges())
    for l in range(10):
        Gud = nx.to_undirected(G)
        C =minimum_cycle_basis(Gud)
        res = print_cycle_statistics(G,cycles = C)
        res = ["Price_m_3"] + list(res.values())
        df_list.append(res)
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
        
    G = price_dag(N,5,1,N+1)
    tr(G)
    edges  = list(G.edges())
    for l in range(10):
        Gud = nx.to_undirected(G)
        C =minimum_cycle_basis(Gud)
        res = print_cycle_statistics(G,cycles = C)
        res = ["Price_m_5"] + list(res.values())
        df_list.append(res)
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
        
        
    data = pd.DataFrame(df_list)
    data.columns =["name","E_p","E_p_std","S_max","S_mean","S_std","C","null_L","E_TR","L_max","N","b","D","M","lambda_max","h_mean","h_std","s_mean","s_std","lambda_max_byC"]
    data.index = data.name
    x = range(len(set(data.index)))
    for c in data.columns:
        if c!="name":
            y =  data.groupby(data.index)[c].mean()
            yerr = data.groupby(data.index)[c].std()
            plt.errorbar(x,y,yerr,fmt = "o",capsize= 5)
            plt.xlabel("Network model")
            plt.ylabel(c)
            plt.xticks(x,data.groupby(data.index)[c].mean().index, rotation='vertical')
            
            if fig_save:
                plt.tight_layout()
                plt.savefig("quasiunicity_N_eq_{}_observable_{}.pdf".format(N,c))
                
            plt.show()
            
   
            
   
def plot_eigenspectrum(ntype, N,fig_save = False):
    N=50
    plt.figure(figsize=  (6,5))
    symbol = ["o","^",">","<","*","v","+","x","d","1"]
    if ntype =="random":
        pvalues = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    elif ntype =="price":
        pvalues = [2,3,4,5,6,7,8]
    for p in pvalues:
        print(p)
        r = 0
        while r<1:
            if ntype =="random":
                G = random_dag(N,p)
            elif ntype =="price":
                G = price_dag(N,p,1,N+1)
            if nx.is_weakly_connected(G) ==True:  
                tr(G)
                Gud = nx.to_undirected(G)
                C = minimum_cycle_basis(Gud)
                egvals = eigenspectrum(G,C)
                plt.scatter([i/len(C) for i in range(len(C))],egvals,marker = symbol[pvalues.index(p)],label="p="+str(p))
                r+=1
           
    plt.xlabel("$1-i/d$",fontsize = 17)
    plt.ylabel("$\\lambda^C_i$",fontsize = 17)
    plt.tick_params(axis= "both",labelsize=17)
    plt.legend(fontsize =14)
    plt.tight_layout()
    if fig_save:
        plt.savefig(ntype+"_"+str(N)+"_eigenspectrum.pdf")


            
if __name__ =="__main__":
    
    
    calc_mcb("price",50)
    plot("dags_cycle_data_price_no_nodes_eq_50.txt",True)
    
    
    plot_eigenspectrum("price",50)
    




    