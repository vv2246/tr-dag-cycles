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
    
    
    #calc_mcb("lattice",50)
    #plot("dags_cycle_data_lattice_no_nodes_eq_50.txt",True)
    
    
    plot_eigenspectrum("price",50)
    




                    
# if task == "plot_eigenspectrum":
    
#     ntype = "random"
#     N =500
#     ctype = "depina"
    
#     if ntype =="random":
#         pvalues = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
#     elif ntype =="price":
#         pvalues = [2,3,4,5,6,7,8]
#     elif ntype =="lattice" or ntype =="russian_doll":
#         pvalues = [5,10,15,20,25,30]#[int(N**0.5)]
#     filename = ntype+"_dags_cycle_data_"+ctype +"_no_nodes_eq_"+str(N)+"_0819.txt"
#     plt.figure(figsize=  (6,5))
#     symbol = ["o","^",">","<","*","v","+","x","d","1"]
#     egvals_at_p = {}
#     for p in pvalues:
#         print(p)
#         r = 0
#         sumegval = []
#         while r <1:
#             if ntype =="random":
#                     G = random_dag(N,p)
#             elif ntype =="price":
#                 G = price_dag(N,p,1,N+1)
#             elif ntype =="lattice":
#                 G = lattice_dag_2D(p,p)
#             elif ntype =="russian_doll":
#                 G = russian_doll_dag(p*p)
#             tr(G)
#             try:
#                 Gud = nx.to_undirected(G)
#                 C = minimum_cycle_basis(Gud)
                
#                 egvals = eigenspectrum(G,C)
#                 sumegval.append(sum(egvals)/len(egvals))
#             except:
#                 pass
#         egvals_at_p[p] = sumegval
#         plt.scatter([i/len(C) for i in range(len(C))],egvals,marker = symbol[pvalues.index(p)],label="p="+str(p))
             
#     #plt.scatter([a for a in egvals_at_p.keys()],[np.mean(x) for x in egvals_at_p.values()],[np.std(x) for x in egvals_at_p.values()],capsize = 5,fmt=".")
#     plt.xlabel("$1-i/d$",fontsize = 17)
#     plt.ylabel("$\\lambda^C_i$",fontsize = 17)
#     plt.tick_params(axis= "both",labelsize=17)
#     plt.legend(fontsize =14)
#     #plt.ylim(4,7)
#     plt.tight_layout()
#     plt.savefig(ntype+"_"+str(N)+"_eigenspectrum.pdf")
     
               
# if task=="plot":
    
#     ##########
#     #FOR PRICE
#     ##########
    
#     data = pd.read_csv("price_dags_cycle_data_depina_no_nodes_eq_500_0826.txt",sep = "\t")
#     #data = pd.read_csv("random_dags_cycle_data_DAG_no_nodes_eq_500_1005.txt",sep = "\t")
#     data["Density"] = data['Number of TR edges']/(250*499)
    
#     xlabel = "m"
#     filestart = "price"
#     mean_list = []
#     for m in set(data["m"]):
#         mean_list.append(list(data[data["m"]==m].mean()))
    
#     mean_df = pd.DataFrame(mean_list)
#     mean_df.columns = data.columns[:]
    
#     std_list = []
#     for m in set(data["m"]):
#         std_list.append([m]+list(data[data["m"]==m].std())[1:])    
#     mean_df.index = mean_df.m
#     std_df = pd.DataFrame(std_list)
#     std_df.columns = data.columns[:]
#     std_df.index = std_df.m
    
    
#     """
    
#     fig, ax = plt.subplots(1,1,figsize = (5,4))
#     ax.errorbar(mean_df["m"],mean_df["Density"],std_df["Density"],fmt= "o",capsize=5,label = "Random",color="red")
#     data = pd.read_csv("price_dags_cycle_data_depina_no_nodes_eq_500_0826.txt",sep = "\t")
#     data["Density"] = data['Number of TR edges']/(250*499)
#     xlabel = "m"
#     filestart = "random"
#     mean_list = []
#     for m in set(data["m"]):
#         mean_list.append(list(data[data["m"]==m].mean()))
    
#     mean_df = pd.DataFrame(mean_list)
#     mean_df.columns = data.columns[:]
    
#     std_list = []
#     for m in set(data["m"]):
#         std_list.append(list(data[data["m"]==m].std()))    
#     mean_df.index = mean_df.m
#     std_df = pd.DataFrame(std_list)
#     std_df.columns = data.columns[:]
#     std_df.index = std_df.m
    
#     ax2 = ax.twiny()  # instantiate a second axes that shares the same x-axis
#     ax2.set_xlabel(xlabel,fontsize = 18,c="g")
#     ax2.tick_params(axis= "y",labelsize=18,labelcolor="g")
#     ax2.errorbar(mean_df["m"],mean_df["Density"],std_df["Density"],fmt= "o",capsize=5,label = "Price",color="green")
#     ax.set_ylabel("$2E_{\\mathrm{TR}}/N(N-1)$",fontsize = 18)
#     ax.set_xlabel("p",fontsize = 18,c="r")
#     ax2.tick_params(axis= "y",labelsize=18)
#     ax.tick_params(axis ="x", labelsize=18)
#     ax.tick_params(axis= "both",labelsize=18)
#     ax2.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig("density_price_random_N_eq_500.pdf")
    
    
#     """
    
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df["m"],mean_df["Number of TR edges"],std_df["Number of TR edges"],fmt= "o",capsize=5,)#label = "Mixers")
#     #plt.legend(fontsize = 14)
#     plt.ylabel("$E_{\\mathrm{TR}}$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"_E_TR_N_eq_500.pdf")
    
    
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df["m"],mean_df["Mean cycle height"],std_df["Mean cycle height"],fmt= "o",capsize=5,)#label = "Mixers")
#     #plt.legend(fontsize = 14)
#     plt.ylabel("$\\langle h \\rangle$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"_h_N_eq_500.pdf")
    
 
 
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df["m"],mean_df["Mean edge participation"],std_df["Mean edge participation"],fmt= "o",capsize=5,)#label = "Mixers")
#     #plt.legend(fontsize = 14)
#     plt.ylabel("$E_P$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"_edge_participation_N_eq_500.pdf")
 
    
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df["m"],mean_df["Mean cycle size"],std_df["Mean cycle size"],fmt= "o",capsize=5,)#label = "Mixers")
#     #plt.legend(fontsize = 14)
#     plt.ylabel("$\\langle C \\rangle$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"_mean_cycle_size_N_eq_500.pdf")
    
        
 
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df["m"],mean_df['Number_diamonds'],std_df["Number_diamonds"],fmt= "o",capsize=5,label = "Diamonds")

#     plt.errorbar(mean_df["m"],mean_df['Number_mixers'],std_df["Number_mixers"],fmt= "o",capsize=5,label = "Mixers")
#     plt.legend(fontsize = 18)
#     plt.ylabel("Number of cycles",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"_number_mixers_diamonds_N_500.pdf")
    
#     fig, ax = plt.subplots(1,1,figsize = (5,4))
#     ax.errorbar(mean_df["m"],mean_df['Largest eigenvalue M_C'],std_df['Largest eigenvalue M_C'],fmt= "o",color = "b",capsize=5)#label = "Mixers")
#     #plt.legend(fontsize = 14)
#     ax.set_ylabel("$\\lambda^C_{\\mathrm{max}}$",fontsize = 16, color = "b")
#     ax.set_xlabel(xlabel,fontsize = 18)
#     ax.tick_params(axis= "y",labelsize=18)#, color = "b")
    
#     ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
#     ax2.errorbar(mean_df["m"]+0.2,mean_df['Largest M_C by corr size'],std_df['Largest M_C by corr size'],fmt= "o",color="g",capsize=5)#label = "Mixers")
#     #plt.legend(fontsize = 14)
#     ax2.set_ylabel("$\\lambda^C_{\\mathrm{max}}/C$",fontsize = 16, color = "g")
#     ax2.set_xlabel(xlabel,fontsize = 18)
#     ax2.tick_params(axis= "y",labelsize=18)
#     ax.tick_params(axis ="x", labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"_eigvals_N_eq_500.pdf")
    
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df['m'],mean_df["Number cycle connected components"],std_df["Number cycle connected components"], fmt= "o",capsize=5,)#label = "Mixers")
#     #plt.ylim(0.7,1.4)
#     #plt.legend(fontsize = 14)
#     plt.ylabel("null$L^C$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"_null_L_N_eq_500.pdf")
 
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df['m'],mean_df["Std cycle size"],std_df["Std cycle size"],fmt= "o",capsize=5,)#label = "Mixers")
#     #plt.ylim(0.9,1.1)
#     #plt.legend(fontsize = 14)
#     plt.ylabel("$\\sigma(C)$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"std_cycle_size_N_eq_500.pdf")
    
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df['m'],mean_df["Mean cycle size"],std_df["Mean cycle size"],fmt= "o",capsize=5,)#label = "Mixers")
#     #plt.ylim(0.9,1.1)
#     #plt.legend(fontsize = 14)
#     plt.ylabel("$\\langle C \\rangle$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"mean_cycle_size_N_eq_500.pdf")
    
    
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df['m'],mean_df["Largest cycle size"],std_df["Largest cycle size"],fmt= "o",capsize=5,)#label = "Mixers")
#     plt.ylim(5,10)
#     #plt.legend(fontsize = 14)
#     plt.ylabel("$ C_{\\mathrm{max}}$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"max_cycle_size_N_eq_500.pdf")

    
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df['m'],mean_df["Mean stretch"],std_df["Mean stretch"],fmt= "o",capsize=5,)#label = "Mixers")
#     #plt.ylim(0.9,1.1)
#     #plt.legend(fontsize = 14)
#     plt.ylabel("$\\langle s \\rangle$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"mean_stretch_size_N_eq_500.pdf")
    
#     plt.figure(figsize =(5,4))
#     plt.errorbar(mean_df['m'],mean_df["Balance"],std_df["Balance"],fmt= "o",capsize=5,)#label = "Mixers")
#     plt.ylim(0.0,0.2)
#     #plt.legend(fontsize = 14)
#     plt.ylabel("$\\langle b \\rangle$",fontsize = 18)
#     plt.xlabel(xlabel,fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig(filestart+"balance_size_N_eq_500.pdf")
    
    
    
    

#     plt.figure(figsize =(5,4))
#     df = pd.read_csv("lattice_dags_cycle_data_horton_310520.txt",sep ="\t")
#     df.index = df["number of nodes "]
#     df["Largest eigenvalue M_C"].plot(style ="-o",label  = "Lattice")
#     df = pd.read_csv("russian_doll_dags_cycle_data_horton_310520.txt",sep ="\t")
#     df.index = df["number of nodes "]
#     df["Largest eigenvalue M_C"].plot(style ="-o",label  = "Russian doll")
#     plt.ylabel("$\\lambda^C_{\\mathrm{max}}$",fontsize = 18)
#     plt.xlabel("Number of nodes",fontsize = 18)
#     plt.legend(fontsize = 18)
#     plt.tick_params(axis= "both",labelsize=18)
#     plt.tight_layout()
#     plt.savefig("eigenvalues_lattice_russian_doll.pdf")
    
    
    
    