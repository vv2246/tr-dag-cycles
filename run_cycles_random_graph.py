# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 11:08:56 2019

@author: Vaiva
"""



from cycle_utilities import *


plt.rcParams.update({'font.size': 20})
        
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
    #N =50
    if ntype =="random":
        pvalues = [0.1]#,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
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
                    
                    
                    
 
def plot_figure_7(filename,fig_save = False,nruns= 20, labels = {'p':"$p$",
    'E':"$E$",
     'Mean edge participation':"$\\langle E_p\\rangle$",
      'Std edge participation':"$\\sigma(E_p)$",
       'Largest cycle size': "$C_{\\mathrm{max}}$",
        'Mean cycle size':"$\\langle C\\rangle$", 
        'Std cycle size':"$\\sigma(C)$",
       'Number of cycles':"d",
        'Number cycle connected components':"null$L^C$",
       'Number of TR edges':"$E_{TR}$", 
       'Longest path':"L",
        'Number of nodes':"N",
       'Number of edges':"E", 
       'Balance':"$\\langle b\\rangle$",
        'Balance_norm':"$\\langle b_{norm}\\rangle$",
         'Number_diamonds':"Diamonds",
       'Number_mixers':"Mixers",
        'Largest eigenvalue M_C':"$\\lambda_{\\mathrm{max}}^c$", 
        'Largest eigenvalue A_C':"labda_A_c",
       'Mean cycle height':"$\\langle h\\rangle$", 
       'Std cycle height':"$\\sigma( h)$", 
       'Mean stretch': "$\\langle s\\rangle$",
       'eigenvalue_ratio':'eigenvalue_ratio', 
       'Largest M_C by corr size':'Largest M_C by corr size'},error_in_mean = True):
    
    data = pd.read_csv(filename,sep="\t")
    print(data.columns)
    #data.columns = ["p","E","E_p","E_p_std","S_max","S_mean","S_std","C","null_L","E_TR","L_max","N","b","D","M","lambda_max","h_mean","h_std","s_mean","s_std","lambda_max_byC","Nan"]
    if "price" in filename:
        name="price"
        xlabel = "$m$"
        data.index = data.p
        x = data.groupby(data.index)["p"].mean()
        shift = 0.1
    elif "random" in filename:
        name="random"
        data.index = data.p
        x = data.groupby(data.index)["p"].mean()
        xlabel = "$p$"
        shift= 0.01
            
    figure_filename = name+"_{}.pdf"
    if error_in_mean==True:
        denom = np.sqrt(nruns)
    else:
        denom = 1
    for c in data.columns:
        plt.figure(figsize=(6,5))
        y =  data.groupby(data.index)[c].mean()
        yerr = data.groupby(data.index)[c].std()/denom
        if "price" in filename or "random" in filename:
            plt.errorbar(x,y,yerr,fmt = "o",capsize= 5,color= "teal")
        else:
            plt.errorbar(x,y,fmt="o",capsize= 0)
        plt.xlabel(xlabel)
        plt.ylabel(labels[c])
        if fig_save:
            plt.tight_layout()
            plt.savefig(figure_filename.format(c))
        
        plt.show()
        
    if ("price" in filename) or ("random" in filename):
        
        plt.figure(figsize=(6,5))
        y =  data.groupby(data.index)["Number_diamonds"].mean()
        yerr = data.groupby(data.index)["Number_diamonds"].std()/denom
        plt.errorbar(x,y,yerr,fmt = "o",capsize= 5,label = "Diamonds",color="teal")
        y =  data.groupby(data.index)["Number_mixers"].mean()
        yerr = data.groupby(data.index)["Number_mixers"].std()/denom
        plt.errorbar(np.array(x)+shift,y,yerr,fmt = "o",capsize= 5,label = "Mixers",color="goldenrod")
        plt.ylabel("Number of cycles")
        plt.legend()
        plt.xlabel(xlabel)
    
        if fig_save:
            plt.tight_layout()
            if "price" in filename:
                plt.savefig("price_number_mixers_diamonds_N_500.pdf")
            if "random" in filename:
                
                plt.savefig("random_number_mixers_diamonds_N_500.pdf")
        
        plt.show()
        
        
        fig,ax = plt.subplots(figsize=(6,5))
        y =  data.groupby(data.index)['Largest eigenvalue M_C'].mean()
        yerr = data.groupby(data.index)['Largest eigenvalue M_C'].std()/denom
        ax.errorbar(x,y,yerr,fmt = "o",capsize= 5,color= "teal")
        y =  data.groupby(data.index)['Largest M_C by corr size'].mean()
        yerr = data.groupby(data.index)['Largest M_C by corr size'].std()/denom
        ax2 = ax.twinx()
        ax2.errorbar(np.array(x)+shift,y,yerr,fmt = "o",capsize= 5,color="goldenrod")
        ax.set_ylabel("$\\lambda^C_{\\mathrm{max}}$",color= "teal")
        ax2.set_ylabel("$\\lambda^C_{\\mathrm{max}}/C$",color = "goldenrod")
        ax.set_xlabel(xlabel)
        
        
    
        if fig_save:
            plt.tight_layout()
            if "price" in filename:
                plt.savefig("price_eigvals_N_eq_500.pdf")
            if "random" in filename:
                
                plt.savefig("random_eigvals_N_eq_500.pdf")
        
        plt.show()
            


def plot_figure_9(fn_random,fn_price,fig_save = False,nruns= 20 ,error_in_mean = True):
    
    
    if error_in_mean==True:
        denom = np.sqrt(nruns)
    else:
        denom = 1
    fig,ax = plt.subplots(figsize=(6,5))
    df = pd.read_csv(fn_random,sep="\t")
    print(df.columns)
    df.index = df.p
    xlabel = "$p$"
    x = df.groupby(df.index)["p"].mean()
    s =df['Number of TR edges']*2/(df['Number of nodes']*(df['Number of nodes']-1))
    y =  s.groupby(s.index).mean()
    yerr = s.groupby(s.index).std()/denom
    ax.errorbar(x,y,yerr,fmt = "o",capsize= 5,color= "teal")
    ax.set_xlabel(xlabel,color= "teal")
    
    df = df[df.p ==0.3]
    print(df.mean())
    print(df.std()/denom)
    
    
    ax2 = ax.twiny()
    df = pd.read_csv(fn_price,sep="\t")
    df.index = df.p
    x = df.groupby(df.index)["p"].mean()
    xlabel = "$m$"
    s =df['Number of TR edges']*2/(df['Number of nodes']*(df['Number of nodes']-1))
    y =  s.groupby(s.index).mean()
    yerr = s.groupby(s.index).std()/denom
    ax2.errorbar(x,y,yerr,fmt = "o",capsize= 5,color= "goldenrod")
    ax2.set_xlabel(xlabel,color= "goldenrod")
    ax.set_ylabel("$2E_{\\mathrm{TR}}/N(N-1)$")
    plt.tight_layout()
    if fig_save:
        plt.savefig("density_price_random_N_eq_500.pdf")
    
    df = df[df.p ==4]
    print(df.mean())
    print(df.std()/denom)
    
    
        
def test_quasiunicity(N,fig_save = False,nruns=20, labels = {'p':"$p$",
    'E':"$E$",
     'Mean edge participation':"$\\langle E_p\\rangle$",
      'Std edge participation':"$\\sigma(E_p)$",
       'Largest cycle size': "$C_{\\mathrm{max}}$",
        'Mean cycle size':"$\\langle C\\rangle$", 
        'Std cycle size':"$\\sigma(C)$",
       'Number of cycles':"d",
        'Number cycle connected components':"null$L$",
       'Number of TR edges':"$E_{TR}$", 
       'Longest path':"L",
        'Number of nodes':"N",
       'Number of edges':"E", 
       'Balance':"$\\langle b\\rangle$",
        'Balance_norm':"$\\langle b_{norm}\\rangle$",
         'Number_diamonds':"Diamonds",
       'Number_mixers':"Mixers",
        'Largest eigenvalue M_C':"$\\lambda_{\\mathrm{max}}^c$", 
        'Largest eigenvalue A_C':"labda_A_c",
       'Mean cycle height':"$\\langle h\\rangle$", 
       'Std cycle height':"$\\sigma( h)$", 
       'Mean stretch': "$\\langle s\\rangle$",
       'eigenvalue_ratio':'eigenvalue_ratio', 
       'Largest M_C by corr size':'Largest M_C by corr size'},error_in_mean = True):
    
    #N=50
    
    G = random_dag(N,0.1)
    tr(G)
    edges  = list(G.edges())
    df = pd.DataFrame()
    df_list = []
    for l in range(nruns):
        Gud = nx.to_undirected(G)
        C =minimum_cycle_basis(Gud)
        res = print_cycle_statistics(G,cycles = C)
        res = ["Random $p=0.1$"] + list(res.values())
        df_list.append(res)
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
        
    G = random_dag(N,0.8)
    tr(G)
    edges  = list(G.edges())
    for l in range(nruns):
        Gud = nx.to_undirected(G)
        C =minimum_cycle_basis(Gud)
        res = print_cycle_statistics(G,cycles = C)
        res = ["Random $p=0.8$"] + list(res.values())
        df_list.append(res)
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
        
        
    G = price_dag(N,3,1,N+1)
    tr(G)
    edges  = list(G.edges())
    for l in range(nruns):
        Gud = nx.to_undirected(G)
        C =minimum_cycle_basis(Gud)
        res = print_cycle_statistics(G,cycles = C)
        res = ["Price $m=3$"] + list(res.values())
        df_list.append(res)
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
        
    G = price_dag(N,5,1,N+1)
    tr(G)
    edges  = list(G.edges())
    for l in range(nruns):
        Gud = nx.to_undirected(G)
        C =minimum_cycle_basis(Gud)
        res = print_cycle_statistics(G,cycles = C)
        
        print(res.keys())
        res = ["Price $m=5$"] + list(res.values())
        df_list.append(res)
        G = nx.DiGraph()
        random.shuffle(edges)
        G.add_edges_from(edges)
        
    data = pd.DataFrame(df_list)
    data.columns =["name","$E_p$","$\\sigma (E_p)$","S_max","S_mean","S_std","C","null_L","E_TR","L_max","N","$\\langle b\\rangle$","D","M","$\\lambda_{\\mathrm{max}}^c$","h_mean","h_std","$\\langle s \\rangle$","s_std","lambda_max_byC"]
                     
    data.index = data.name
    if error_in_mean:
        denom = np.sqrt(nruns)
    else:
        denom =1
    x = range(len(set(data.index)))
    for c in data.columns:
        if c!="name":
            y =  data.groupby(data.index)[c].mean()
            yerr = data.groupby(data.index)[c].std()/denom
            plt.figure(figsize=(6,6))
            plt.errorbar(x,y,yerr,fmt = "o",capsize= 5,color="blue")
            plt.xlabel("Network model")
            plt.ylabel(c)
            plt.xticks(x,data.groupby(data.index)[c].mean().index, rotation='vertical')
            
            if fig_save:
                plt.tight_layout()
                plt.savefig("quasiunicity_N_eq_{}_observable_{}.pdf".format(N,c))
                
            plt.show()
            
   
            
   
    
   
def plot_eigenspectrum(ntype, N,fig_save = False):
    #N=50
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
    
    
    #calc_mcb("random",500)
    plot_figure_7("price_dags_cycle_data_depina_no_nodes_eq_500_1005.txt",True,20,error_in_mean=True)
    #plot_figure_9("random_dags_cycle_data_DAG_no_nodes_eq_500.txt","price_dags_cycle_data_depina_no_nodes_eq_500_1005.txt",True)
    
    #plot_eigenspectrum("price",50)
    
    #test_quasiunicity(100,fig_save = True,nruns=10,error_in_mean = True)



    