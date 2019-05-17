#! /usr/bin/env python
'''
@author: Jie Hu
'''
import networkx as nx
import numpy as np
import copy
import scipy.special
from collections import defaultdict
import csv
import sys

#read file
def read_input(network_file,seed_file):

    sniffer = csv.Sniffer()
    line_delimiter = None
    for line in open(network_file,'r'):
        if line[0]=='#':
            continue
        else:
            dialect = sniffer.sniff(line)
            line_delimiter = dialect.delimiter
            break
    if line_delimiter == None:
        print ('network_file format not correct')
        sys.exit(0)
    #read network
    G = nx.Graph()
    with open(network_file,'r') as f1:
        for line in f1:
            if line[0]=='#':
                continue
            line_data = line.strip().split(line_delimiter)
            node1 = line_data[0]
            node2 = line_data[1]
            G.add_edge(node1,node2)
    #read seed genes
    seed_genes = set()
    for line in open(seed_file,'r'):
        if line[0]=='#':
            continue
        line_data = line.strip().split('\t')
        seed_gene = line_data[0]
        seed_genes.add(seed_gene)

    return G,seed_genes

#n!
def compute_all_gamma_ln(N):

    gamma_ln = {}
    for i in range(1,N+1):
        gamma_ln[i] = scipy.special.gammaln(i)

    return gamma_ln

#lg
def logchoose(n, k, gamma_ln):

    if n-k+1 <= 0:
        return scipy.infty
    lgn1  = gamma_ln[n+1]
    lgk1  = gamma_ln[k+1]
    lgnk1 = gamma_ln[n-k+1]

    return lgn1 - [lgnk1 + lgk1]

#lg+-
def gauss_hypergeom(x, r, b, n, gamma_ln):

    return np.exp(logchoose(r, x, gamma_ln)+logchoose(b, n-x, gamma_ln)-logchoose(r+b, n, gamma_ln))

#pvalue
def pvalue(kb, k, N, s, gamma_ln):

    p = 0.0
    for n in range(kb,k+1):
        if n > s:
            break
        prob = gauss_hypergeom(n, s, N-s, k, gamma_ln)
        p += prob

    if p > 1:
        return 1
    else:
        return p

#get neighbors and degrees
def get_neighbors_and_degrees(G):

    neighbors,all_degrees = {},{}
    for node in G.nodes():
        neighbors[node] = set(G.neighbors(node))
        all_degrees[node] = G.degree(node)

    return neighbors,all_degrees

#Reduce number of calculations
def reduce_not_in_cluster_nodes(nGsub,all_degrees,neighbors,cluster_nodes,not_in_cluster):

    reduced_not_in_cluster = {}
    kb2k = defaultdict(dict)
    for node in not_in_cluster:
        kd = all_degrees[node]
        kb = 0
        ks= 0
        temps = copy.copy(nGsub)
        for neighbor in neighbors[node]:
            if neighbor in cluster_nodes:
                ks += 1
            for temp in temps[:]:
                if neighbor in temp:
                    kb+= 1
                    temps.remove(temp)
                    continue
        k=kd-ks+kb
        kb2k[kb][k] =node
    k2kb = defaultdict(dict)
    for kb,k2node in kb2k.items():
        min_k = min(k2node.keys())
        node = k2node[min_k]
        k2kb[min_k][kb] = node
    for k,kb2node in k2kb.items():
        max_kb = max(kb2node.keys())
        node = kb2node[max_kb]
        reduced_not_in_cluster[node] =(max_kb,k)

    return reduced_not_in_cluster

def reduce_not_in_cluster_edges(nGsub,all_degrees,neighbors,cluster_nodes,bG):

    reduced_not_in_cluster1 = {}
    kb2k = defaultdict(dict)
    for edge in bG.edges():
        kd = all_degrees[edge[0]]+all_degrees[edge[1]]-2
        kb1 = []
        kb2 = []
        ks1=0
        ks2=0
        temps1 = copy.copy(nGsub)
        for neighbor in neighbors[edge[0]]:
            if neighbor in cluster_nodes:
                ks1 += 1
            for temp in temps1[:]:
                if neighbor in temp:
                    kb1.append(temp)
                    temps1.remove(temp)
                    continue
        temps2 = copy.copy(nGsub)
        for neighbor in neighbors[edge[1]]:
            if neighbor in cluster_nodes:
                ks2 += 1
            for temp in temps2[:]:
                if neighbor in temp:
                    kb2.append(temp)
                    temps2.remove(temp)
                    continue
        for i in kb2:
            if i not in kb1:
                kb1.append(i)
        kb=len(kb1)
        k=kd-ks1-ks2+kb
        kb2k[kb][k] =edge
    k2kb = defaultdict(dict)
    for kb,k2edge in kb2k.items():
        min_k = min(k2edge.keys())
        edge = k2edge[min_k]
        k2kb[min_k][kb] = edge
    for k,kb2edge in k2kb.items():
        max_kb = max(kb2edge.keys())
        edge = kb2edge[max_kb]
        reduced_not_in_cluster1[edge] =(max_kb,k)

    return reduced_not_in_cluster1

#iteration
def iteration(G,S):

    N = G.number_of_nodes()
    C3_proteins=[]
    neighbors,all_degrees = get_neighbors_and_degrees(G)
    cluster_nodes = set(S)
    not_in_cluster = set()
    s0 = len(cluster_nodes)
    gamma_ln = compute_all_gamma_ln(N+1)
    for node in cluster_nodes:
        not_in_cluster |= neighbors[node]
    not_in_cluster -= cluster_nodes
    nG = G.subgraph(cluster_nodes)
    nGsub=list(nx.connected_components(nG))
    n=m=len(nGsub)
    all_p = {}
    while(True):
        pmin = 10
        next_node = 'nix'
        reduced_not_in_cluster = reduce_not_in_cluster_nodes(nGsub,all_degrees,neighbors,cluster_nodes,not_in_cluster)
        for node,kbk in reduced_not_in_cluster.items():
            kb,k = kbk
            try:
                p = all_p[(k,kb,n)]
            except KeyError:
                p = pvalue(kb, k, N-s0+n, n, gamma_ln)
                all_p[(k,kb,n)] = p
            if p < pmin:
                pmin = p
                next_node = node
        if pmin < 0.05:
            C3_proteins.append(next_node)
            # Updating
            cluster_nodes.add(next_node)
            s0 = len(cluster_nodes)
            not_in_cluster1=not_in_cluster
            not_in_cluster = not_in_cluster|( neighbors[next_node] - cluster_nodes )
            not_in_cluster.remove(next_node)
            nG1 = G.subgraph(cluster_nodes)
            nGsub=list(nx.connected_components(nG1))
            m=len(nGsub)
            nG1c = max(nx.connected_component_subgraphs(nG1), key=len)
            overlap=len(nG1c.node()&S)
        if m==n :
            C3_proteins.remove(next_node)
            cluster_nodes.remove(next_node)
            s0 = len(cluster_nodes)
            nG2 = G.subgraph(cluster_nodes)
            nGsub=list(nx.connected_components(nG2))
            n=m=len(nGsub)
            bG = nx.Graph()
            for node in not_in_cluster1:
                for neighbor in neighbors[node]:
                    if neighbor not in cluster_nodes:
                        bG.add_edge(node,neighbor)
            pmin = 10
            next_edge1 = 'one'
            next_edge2='two'
            reduced_not_in_cluster1 = reduce_not_in_cluster_edges(nGsub,all_degrees,neighbors,cluster_nodes,bG)
            for edge,kbk in reduced_not_in_cluster1.items():
                kb,k = kbk
                try:
                    p = all_p[(k,kb,n)]
                except KeyError:
                    p = pvalue(kb, k, N-s0+n, n, gamma_ln)
                    all_p[(k,kb,n)] = p
                if p < pmin:
                    pmin = p
                    next_edge1 = edge[0]
                    next_edge2 = edge[1]
            if pmin < 0.05:
                C3_proteins.append(next_edge1)
                C3_proteins.append(next_edge2)
                # Updating
                cluster_nodes.add(next_edge1)
                cluster_nodes.add(next_edge2)
                s0 = len(cluster_nodes)
                not_in_cluster1 |= ( neighbors[next_edge1] - cluster_nodes )
                not_in_cluster1 |= ( neighbors[next_edge2] - cluster_nodes )
                nG3 = G.subgraph(cluster_nodes)
                nGsub=list(nx.connected_components(nG3))
                m=len(nGsub)
                if next_edge1 in not_in_cluster1:
                    not_in_cluster1.remove(next_edge1)
                if next_edge2 in not_in_cluster1:
                    not_in_cluster1.remove(next_edge2)
                not_in_cluster=not_in_cluster1
                nG3c = max(nx.connected_component_subgraphs(nG3), key=len)
                overlap=len(nG3c.node()&S)
            if m==n:
                cluster_nodes.remove(next_edge1)
                if next_edge2 in cluster_nodes:
                    cluster_nodes.remove(next_edge2)
                C3_proteins.remove(next_edge1)
                C3_proteins.remove(next_edge2)
                break
        n=m
    mG = G.subgraph(cluster_nodes)
    mGc = max(nx.connected_component_subgraphs(mG), key=len)
    overlap=list(mGc.node()&S)

    return C3_proteins,mG,mGc,overlap

def C3(G,seed_genes):

    #1.Locating  the seed genes that are in the network
    all_genes_in_network = set(G.nodes())
    seed_genes = set(seed_genes)
    disease_genes = seed_genes & all_genes_in_network
    print ('There are %s disease genes mapped to the network' %(len(disease_genes)))

    #2.algorithm
    C3_proteins,mG,mGc,overlap = iteration(G,disease_genes)

    #3.the results
    return C3_proteins,mG,mGc,overlap

if __name__ == '__main__':

    input_list = sys.argv
    network_file=input_list[1]
    seeds_file=input_list[2]
    outfile_name1='C3_genes.txt'
    outfile_name2='connected_disease_genes.txt'
    G,seed_genes= read_input(network_file,seeds_file)
    C3_proteins,mG,mGc, overlap=C3(G,seed_genes)
    for node in C3_proteins:
        file=open(outfile_name1,'a')
        file.write(str(node)+'\n')
    for node in overlap:
        file=open(outfile_name2,'a')
        file.write(str(node)+'\n')

