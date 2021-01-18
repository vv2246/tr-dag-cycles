"""This contains code to assign a height to every vertex in the DAG.

Height is defined to be the longest length of any longest path from a source node 
       to the target node in question.

This is a copy of the version of this file in the james_dag library, taken 5/3/18

Author: Tim Evans
"""
import networkx as nx

# this can be changed to import any suitable test DAG
#import model_test2 as test

def set_heights(DAG, ordered=False, add_labels=False,add_labels_fractional=False):
    """    Find the heights of all nodes.
    
       Depth first search.  Starts from source nodes, nodes with no predecessors (zero in-degree?).
       Height is defined to be the longest length of any longest path from a source node 
       to the target node in question.
    
       Keyword arguments:
                 DAG -- networkX digraph forming a DAG
             ordered -- if ordered then node1>node2 is necessary for a path to exist from node1 to node2
          add_labels -- True if want to add 'height' label to each node with height value
          add_labels_fractional=False -- True if want 'height_fraction' label added, 1= max height, 0 = source node
       Return
        height -- dictionary such that height[node] is height of node
    """
#    if start == end:
#        return [0, None]
    height = {} # uses a dictionary?
    # find source nodes and set all heights to be negative number to indicate unvisited
    source_nodes=[]
    for node in DAG.nodes():
        height[node] = -1  # i.e. not connected to start
        if DAG.in_degree(node) ==0:
            height[node] = 0
            source_nodes.append(node)   
    if ordered is True:
        source_nodes.sort(reverse=True) # largest first           
    for source in source_nodes:
        height_recursive(DAG, source, height, ordered)    
    if add_labels:
        for node in DAG:
            DAG.node[node]['height']=height[node]
    if add_labels_fractional:
        max_height=float(max(height))
        for node in DAG:
            DAG.node[node]['height_fraction']=height[node]/max_height
    return height


def height_recursive(DAG, current, height, ordered, depth=0):
    """Recursively searches through the DAG assigning heights (distances from current) to nodes.
       
       Assigns the longest path distance inherited from current height list to nodes connected to current.
       Depth first search.
       Adapted from JC code lp_recursive in alg_paths.py.
        
       Keyword arguments:
                 DAG -- networkX digraph forming a DAG
             current -- node currently under investigation, must have a non-negative height
              height -- dictionary where height[n] has the longest path distance for node n (negative value then unvisited)
             ordered -- if ordered then current>=end is necessary for a path to exist
               depth -- current distance from current nodes
    """
    startneighbours= DAG.successors(current)
    if ordered:
        startneighbours.sort(reverse=True) # largest first to favour greedy paths first
#        if current<end:
            # we have gone past the end point
            # so we aren't going to find a path here
            # we can only make this check if the DAG is ordered
#            return None
    for node in startneighbours:
        height_n = height[node]
        height_current = height[current]
        if height_n < height_current + 1:
            height[node] = height_current + 1
            height_recursive(DAG, node, height, ordered, depth+1)
    return None

__name__ = "notmail"
if __name__ == "__main__":
    """Tests to see if height routines are working"""
    print ('Testing alg_height' )
    DAG= test.make_DAG()
    height=set_heights(DAG, ordered=True,add_labels=True)
    print ('node height node[height]' )
    for node in DAG.nodes():
        print (node,height[node],DAG.node[node]['height'] )                 
            
    """   
            
    date = {n:str2date(n) for n in G}
    min_time = min(date.values())
    partTime = defaultdict(list)
    for key,val in Ppart.items():
        partTime[key] = [(date[n]-min_time) for n in val]
        time_variation = {}
        year = 365.25*(8.64*10**10)
    for p, timedeltas in partTime.items():
        if len(timedeltas) >1:
            microtimes = [timedelta_to_microtime(td) for td in timedeltas]
            std = calc_std(microtimes)
            time_variation[p] = [ np.mean(microtimes)/(year),std/(year)]
        else:
            microtimes = [timedelta_to_microtime(td) for td in timedeltas] 
            std = calc_std(microtimes)
            time_variation[p] =[microtimes[0]/year for a in range(1)]+[std]
    x,y = [],[]
    for key,val in time_variation.items(): 
        x.append(val[0]) 
        y.append(len(Ppart[key]))
    plt.scatter(x,y)#,label="$\mathcal{A}\geq 10$")#plt.colorbar(label = "$|\mathcal{A}|$")
    plt.xlabel("$<t_{\mathcal{A}}>-t_0$",fontsize =14)
    #plt.legend()
    plt.ylabel("$|\mathcal{A}|$",fontsize = 14)
    plt.xscale("log")
    plt.yscale("log")#plt.xlim(9,100)#plt.title(title)
    plt.tight_layout()
    plt.savefig("size_vs_averageT_predecessors_common_bibliographies_antichains_hepth"+".png")
    plt.show()
    
    """
        
    x,y,z = [],[],[]  
    for i in range(len(successor_diversity)):
        x.append(successor_diversity[i][0])
        y.append(successor_diversity[i][1])
        z.append(successor_diversity[i][2])    
        
    plt.errorbar(x,z,xerr=y,capsize =2,fmt=".",label="Successors siblinarity")
    x,y,z = [],[],[]  
    for i in height_div.keys():
        x.append(i)
        y.append(height_div[i])
        
    plt.errorbar(x,y,fmt=".",label="Height")
    plt.legend()
    
    predecessor_diversity = []
    for key,val in partition.items():
        ave_height = np.mean([H[a] for a in val])
        std_height =np.std([H[a] for a in val])
        predecessor_diversity.append((ave_height,std_height,calc_diversity(val,label)))
    
    
    x,y,z = [],[],[]  
    for i in range(len(predecessor_diversity)):
        x.append(predecessor_diversity[i][0])
        y.append(predecessor_diversity[i][1])
        z.append(predecessor_diversity[i][2])    
        
    plt.errorbar(x,z,xerr=y,capsize =2,fmt=".",label="Predecessor siblinarity")
    x,y,z = [],[],[]  
    for i in depth_div.keys():
        x.append(i)
        y.append(depth_div[i])
        
    plt.errorbar(x,y,fmt=".",label="Depth")
    plt.legend()