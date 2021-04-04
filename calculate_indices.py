import networkx
from networkx import algorithms
from networkx.algorithms import isomorphism
import networkx as nx
import math

def get_copy(G):
    H = G.__class__()
    H.add_nodes_from(G)
    H.add_edges_from(G.edges)
    return(H)


class Connectivity:
   
   
    def calculate_connectivities(G,f,var=1):
        return sum(f(G.degree(e[0]),G.degree(e[1]))**var for e in G.edges)

    
    def zagreb(G):
        m1 = sum((G.degree(i))**2 for i in range(len(G)))
        m2 = Connectivity.calculate_connectivities(G, lambda a, b : a * b)
        return m1, m2
    
    
    def modified_zagreb(G):
        m1 = sum(1/((G.degree(i))**2) for i in range(len(G)))
        m2 = Connectivity.calculate_connectivities(G, lambda a, b : 1/(a * b))
        return m1, m2


    def variable_zagreb(G,var):
        m1 = sum((G.degree(i)**var) for i in range(len(G)))
        m2 = Connectivity.calculate_connectivities(G, lambda a, b : (a * b), var)
        return m1, m2
    
    
    def randic(G):
        return Connectivity.calculate_connectivities(G, lambda a, b : 1/math.sqrt((a * b)))
    
    
    def sum_connectivity(G):
        return Connectivity.calculate_connectivities(G, lambda a, b : 1/math.sqrt((a + b)))
    
    
    def general_randic(G,var):
        return Connectivity.calculate_connectivities(G, lambda a, b : (a * b), var)
    
    
    def general_sum_connectivity(G,var):
        return Connectivity.calculate_connectivities(G, lambda a, b : (a + b), var)
    
    
class Path:
    
       
    def balaban(G):
        distances = [sum(nx.shortest_path_length(G, source=i).values()) for i in range(len(G))]
        balaban =  sum(math.sqrt(distances[e[0]]*distances[e[1]]) for e in G.edges)
        balaban *= (G.number_of_edges()/(G.number_of_edges()-len(G)+2))
        return balaban
            
        
    def szeged(G):
        szeged = 0
        for e in G.edges:
            l_u = 0
            l_v = 0
            u =  nx.shortest_path_length(G, source=e[0])
            v =  nx.shortest_path_length(G, source=e[1])
            for i in range (len(G)):
                if u[i] < v[i]:
                    l_u += 1
                elif v[i] < u[i]:
                    l_v += 1
            szeged += l_u * l_v
        return szeged
        
        
    def revised_szeged(G):
        szeged = 0
        for e in G.edges:
            l_u = 0
            l_v = 0
            u =  nx.shortest_path_length(G, source=e[0])
            v =  nx.shortest_path_length(G, source=e[1])
            for i in range (len(G)):
                if u[i] < v[i]:
                    l_u += 1
                elif v[i] < u[i]:
                    l_v += 1
                else:
                    l_u += 0.5
                    l_v += 0.5
            szeged += l_u * l_v
        return szeged
        
        
    def detour(G):
        detour = 0
        for i in range(len(G)):
            for j in range(i):
                detour += max(len(x) for x in (nx.all_simple_paths(G, source=i, target=j)))-1
        return detour


    def harary(G):
        return sum(sum(1/x for x in v_lenghts.values() if x>0) for v_lenghts in dict(nx.shortest_path_length(G)).values())/2

        
    def wiener(G):
        return nx.average_shortest_path_length(G, weight='weight')* len(G)*(len(G)-1) * 0.5


    def pisanski(G):
        GM = algorithms.isomorphism.GraphMatcher(G, G)
        paths_len = 0
        for isomorph in GM.isomorphisms_iter():
            for v1, v2 in isomorph.items():
                paths_len += algorithms.shortest_path_length(G, v1, v2)
        return((paths_len*len(G))/(len(list(GM.isomorphisms_iter())) * 2))
