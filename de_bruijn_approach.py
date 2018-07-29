import dinopy
import collections as col
import networkx as nx
from functools import reduce
import numpy as np
import yaml
import sys

otus = sys.argv[1]

class DeBruijnGraph:
    
    @staticmethod
    def chop(st, k):
        for i in range(len(st)-(k-1)):
            yield (st[i:i+k], st[i:i+k-1], st[i+1:i+k])
    
    class Node:
        def __init__(self, km1mer):
            self.km1mer = km1mer
            self.nin = 0
            self.nout = 0
            self.abundance = 0
            self.otu = None
        
        def isSemiBalanced(self):
            return abs(self.nin - self.nout) == 1
        
        def isBalanced(self):
            return self.nin == self.nout
        
        def isHead(self):
            return self.nin == 0
        
        def isTail(self):
            return self.nout == 0
        
        def __hash__(self):
            return hash(self.km1mer)
        
        def __str__(self):
            return self.km1mer

    def __init__(self, strIter, k):
        self.G = nx.MultiDiGraph()    # multimap from nodes to neighbors
        self.nodes = {} # maps k-1-mers to Node objects
        for st in strIter:
            for kmer, km1L, km1R in self.chop(st[0], k):
                nodeL, nodeR = None, None
                if km1L in self.nodes:
                    nodeL = self.nodes[km1L]
                else:
                    nodeL = self.nodes[km1L] = self.Node(km1L)
                if km1R in self.nodes:
                    nodeR = self.nodes[km1R]
                else:
                    nodeR = self.nodes[km1R] = self.Node(km1R)
                #nodeL.abundance += 0.5
                #nodeR.abundance += 0.5
                nodeL.nout += 1
                nodeR.nin += 1
                nodeL.otu = st[1].split(';')[0]
                nodeR.otu = st[1].split(';')[0]
                nodeL.abundance += float(st[1].split('=')[1][:-1])/2
                nodeR.abundance += float(st[1].split('=')[1][:-1])/2
                self.G.add_node(nodeL)
                self.G.add_node(nodeR)
                #self.G.setdefault(nodeL, []).append(nodeR)
        if len(self.G.nodes()) > 0:
            [self.G.add_edge(L, R) for L in self.G.nodes() for R 
                    in self.G.nodes() if L.km1mer[1:] == R.km1mer[:-1]]
        # Iterate over nodes; count # balanced, semi-balanced, neither
        self.nsemi, self.nbal, self.nneither = 0, 0, 0
        # Keep track of head and tail nodes in the case of a graph with
        # Eularian walk (not cycle)
        self.head, self.tail = [], []
        for node in iter(self.nodes.values()):
            if node.isBalanced():
                self.nbal += 1
            
            elif node.isHead():
                node.abundance = node.abundance*2
                self.head.append(node)
                if node.nin == node.nout - 1:
                    self.nsemi
            
            elif node.isTail():
                node.abundance = node.abundance*2
                self.tail.append(node)
                if node.nin == node.nout + 1:
                    self.nsemi += 1
                    
            # If a head or a tail is semi balanced it doesn't count
            # as semi balanced anymore, this has to be fixed
            elif node.isSemiBalanced():
                if node.nin == node.nout + 1 or node.nin == node.nout - 1:
                    self.nsemi += 1
            else:
                self.nneither += 1
    
    def nnodes(self):
        return len(self.nodes)
    
    def nedges(self):
        return len(self.G)

class kmer_filter:
    def __init__(self, otu_file, k):
        self.reads = dinopy.FastaReader(otu_file)
        self.de_bruijn_dict = self.make_de_bruijn_file(self.reads, k)
        self.high_divergence = self.high_divergence_dict(self.de_bruijn_dict)
        #self.reads_containing_query = self.reads_contianing_query(query_node)
        
    def make_de_bruijn_file(self, reads, k):
        bru_dict = col.defaultdict()
        nodes = []
        for f in reads.entries():
            seq = f.sequence.decode()
            name = f.name.decode()
            for i in range(len(seq)-k+1):
                node = seq[i:i+k-1]
                nodes.append(node)
                if not node in bru_dict.keys():
                    bru_dict[node] = {}
                    bru_dict[node]['name'] = []
                    bru_dict[node]['abu'] = 1
                else:
                    bru_dict[node]['abu'] += 1
                if not name in bru_dict[node]['name']:
                    bru_dict[node]['name'].append(int(name.split('=')[1][:-1]))
        return bru_dict
    
    def reads_contianing_query(self, query_node):
        read_list = []
        for seq in self.reads.entries():
            if query_node in seq.sequence.decode():
                read_list.append((seq.sequence.decode(), seq.name.decode()))
        return read_list
    
    def high_divergence_dict(self, de_bruijn_dict, cutoff=2):
        high_divergence = {}
        for key in de_bruijn_dict.keys():
            if len(de_bruijn_dict[key]['name']) >= cutoff:
                high_divergence[key] = de_bruijn_dict[key]
        return high_divergence
    
    def abundance_filter(self, reads, high_div):
        num_reads = sum(1 for read in reads.entries())
        low_abu = []
        for kmer in high_div:
            if high_div[kmer]['abu'] < num_reads/100:
                low_abu.append(kmer)
        return low_abu
    
class unitig_creator:
    def __init__(self, graph):
        self.unitig_graph = self.maximal_non_branching_paths(graph)
        
    # Creates a unitig from k-mers that have only 1 out degree,
    # as there is no loss in information when combining those
    # it creates a "non-branching-path" out of these k-mers by walking
    # along each node, if it encounters a node with more than one in-degree,
    # it wont add it to the non-branching path, insted,
    # it will break beforehand and will start a new nbp later
    # in the maximal_non_branching_paths algorithm
    # (to circumvent duplicate subsequences and differences in abundances).
    # If the out-degree is higher than 1, the node will be added to
    # the unitig and then the loop breaks, as it still belongs to
    # the unitig that is created. before forming the unitig by
    # reducing the kmers to a single sequence the func checks if the
    # abundance of each node has the same value, if they do not have
    # the same value it would imply that they do not belong to the same
    # sequence and that some mistake/bug happened with the
    # sequence/kmer/abundance.
    def create_unitig(self, node, out_edge, graph):
        non_branching_path = []
        if graph.in_degree(node) == 0:
            non_branching_path.append((out_edge, node.abundance, node.otu))
        next_node = out_edge[1]
        while list(graph.out_edges(next_node)):
            if graph.in_degree(next_node) != 1:
                break
            non_branching_path.append((list(graph.out_edges(next_node))[0],
                next_node.abundance, next_node.otu))
            if graph.out_degree(next_node) != 1:
                break
            next_node = list(graph.out_edges(next_node))[0][1]
        # this is for the case that the current node and the following
        # node have an indegree greater than 1
        if not non_branching_path:
            non_branching_path.append((out_edge, node.abundance, node.otu))
        edge, abundance, otu = list(zip(*non_branching_path))
        #might have to increase k if one of these errors happen
        if len(set(abundance)) != 1:
            raise ValueError('unitig has edges with different abundances')
        if len(set(otu)) !=1:
            raise ValueError('unitig belongs to different otus')
        return reduce(lambda x,y: x+y[-1],[l.km1mer + r.km1mer[-1] for l,
            r in edge]), next_node, abundance, otu

    # Creates a graph containing all unitigs and their corresponding
    # abundances for further processing/analysis. loops through
    # all nodes, creates a new unitig if it is the or one of the first
    # nodes (no indegree) or if the node has > 1 out_degree
    # it creates new unitigs for each braching path.
    # If the next_node has more than 1 in_degree,
    # it will first create the unitig
    # for this path before continuing,
    # to allow multiple incoming edges to this unitig
    # (otherwise it will only create one edge).
    # Keep in mind, and this is important, that the edge abundance
    # is not actually the edge abundance,
    # but the abundance of the target node.
    def maximal_non_branching_paths(self, graph):
        digraph = nx.MultiDiGraph()
        unitig = None
        unitig_src1 = None
        for node in graph:
            if graph.out_degree(node) != 1 or graph.in_degree(node) == 0:
                if graph.out_degree(node) > 0:
                    if unitig:
                            unitig_src1 = unitig
                    for out_edge in list(graph.out_edges(node)):
                        unitig, next_node, abundance, otu = self.create_unitig(
                                node, out_edge, graph)
                        if graph.out_degree(next_node) != 1 or graph.in_degree(
                                next_node) > 1:
                            if unitig_src1:
                                digraph.add_edge(unitig_src1, unitig,
                                        abundance=abundance[0], otu=otu[0])
                        if graph.in_degree(next_node) > 1:
                            unitig_src = unitig
                            if graph.out_degree(next_node) > 0:
                                for out_edge in list(graph.out_edges(next_node)):
                                    unitig, next_node, abundance, otu = self.create_unitig(
                                            next_node, out_edge, graph)
                                    digraph.add_edge(unitig_src, unitig,
                                            abundance=abundance[0], otu=otu[0])
                            else:
                                digraph.add_edge(unitig_src, next_node.km1mer,
                                        abundance=next_node.abundance, otu=next_node.otu)
        return digraph

class chim_checker:
    def __init__(self, potential_chimera_dict, otus):
        self.potential_chimera_dict = potential_chimera_dict
        self.otus = otus
        
    def levenshtein(self, source, target):
        if len(source) < len(target):
            return self.levenshtein(target, source)
        if len(target) == 0:
            return len(source)
        source = np.array(tuple(source))
        target = np.array(tuple(target))
        previous_row = np.arange(target.size + 1)
        for s in source:
            current_row = previous_row + 1
            current_row[1:] = np.minimum(
                    current_row[1:],
                    np.add(previous_row[:-1], target != s))
            current_row[1:] = np.minimum(
                    current_row[1:],
                    current_row[0:-1] + 1)
            previous_row = current_row
        return previous_row[-1]
    
    def find_pot_chimeras(self, unitig_graph, node, chim_dict, already_searched):
        if unitig_graph.out_degree(node) == 1:
            out_node = list(unitig_graph.out_edges(node))[0][1]
            if (node, out_node) in already_searched:
                    pass
            else:
                already_searched.add((node, out_node))
                self.find_pot_chimeras(unitig_graph, out_node, chim_dict,
                        already_searched)
        elif unitig_graph.out_degree(node) > 1:
            target_list = [(target, data) for source, target, data
                    in unitig_graph.out_edges(node, data=True)]
            highest_abu = max(target_list, key=lambda item:item[1]['abundance'])
            for out_node in target_list:
                if (node, out_node[0]) in already_searched:
                    pass
                else:
                #print((levenshtein(highest_abu[0], out_node[0])/out_node[1]['abundance'], len(out_node[0])))
                    if self.levenshtein(highest_abu[0],
                            out_node[0])/out_node[1]['abundance'] > 5:  #
                        if out_node[1]['otu'] in potential_chimera_dict.keys():
                            chim_dict[out_node[1]['otu']] += 1
                        else:
                            chim_dict[out_node[1]['otu']] = 1
                    already_searched.add((node, out_node[0]))
                    self.find_pot_chimeras(unitig_graph, out_node[0],
                            chim_dict, already_searched)
        return chim_dict
    
    def chim_search(self, unitig_graph, potential_chimera_dict):
        source_nodes = [node for node, degree in unitig_graph.in_degree()
                if degree == 0]
        already_searched = set()
        for node in source_nodes:
            self.find_pot_chimeras(unitig_graph, node, potential_chimera_dict,
                    already_searched)
        return potential_chimera_dict
    
    def run(self, kmer_length=30):
        counter = 1
        high_divergence = kmer_filter(otus, kmer_length).high_divergence
        for node in high_divergence.keys():
            read_list = kmer_filter(otus, kmer_length).reads_contianing_query(node)
            d = DeBruijnGraph(read_list, kmer_length)
            unitig_graph = unitig_creator(d.G).unitig_graph
            self.chim_search(unitig_graph, potential_chimera_dict)
            print(str(counter/len(high_divergence.keys())*100) + '% done')
            counter += 1
        return potential_chimera_dict

potential_chimera_dict = dict()
#otus = '../../workflow_test/amplicon_snakemake_pipeline/comparison/new/finalData/merged_representatives.fasta'
h = chim_checker(potential_chimera_dict, otus).run()
with open('potential_chimeras2.yaml', 'w') as f_:
    yaml.dump(h, f_)
