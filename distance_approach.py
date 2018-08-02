import dinopy
import collections as col
import networkx as nx
from functools import reduce
import os
import yaml

# dict with seqs -> kmers
class kmer_filter:
    def __init__(self, otu_file, k, cutoff=10):
        self.reads = dinopy.FastaReader(otu_file)
        self.de_bruijn_dict = self._make_de_bruijn_file(self.reads, k)
        self.kmer_abundance_sorting = self.kmer_abundance_sorting(self.reads, self.de_bruijn_dict, cutoff)
        
    def _make_de_bruijn_file(self, reads, k):
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
    
    def kmer_abundance_sorting(self, reads, bru_dict, cutoff):
        num_reads = sum(1 for read in reads.entries())
        high_abu = []
        for kmer in bru_dict:
            if bru_dict[kmer]['abu'] > num_reads/cutoff:
                high_abu.append(kmer)
        return high_abu
    
    def softmask(self, reads, kmer_abundance_sorting):
        with dinopy.FastaWriter("masked_reads_chims8.fasta", line_width=1000) as faw:
            for entry in reads.entries():
                seq = entry.sequence.decode()
                for kmer in kmer_abundance_sorting:
                    if kmer in seq:
                        seq = seq.replace(kmer, kmer.lower())
                faw.write_entry((seq.encode(), entry.name))

class chimera_search:
    def __init__(self, masked_reads, k):
        self._seq_dict = self._seq_dict(masked_reads, k)
        self._abu_kmer_zip = self._abu_kmer_zip(self._seq_dict)
        self._intersection_list = self._intersection_list(self._abu_kmer_zip[0], self._abu_kmer_zip[1], 
                                                          self._abu_kmer_zip[2], self._abu_kmer_zip[3])
        self._overlap_graph = self._overlap_graph(self._intersection_list)
        self._high_indegree_graph = self._remove_low_indegree_edges(self._overlap_graph)
        self.chimeric_subgraphs = self.subgraphs(self._high_indegree_graph)
        self.potential_chimeras = self.potential_chimeras(self.chimeric_subgraphs)
     
    @staticmethod
    def _lcs(S,T):
        m = len(S)
        n = len(T)
        counter = [[0]*(n+1) for x in range(m+1)]
        longest = 0
        lcs_set = set()
        for i in range(m):
            for j in range(n):
                if S[i] == T[j] and not S[i].islower() and not T[j].islower():
                    c = counter[i][j] + 1
                    counter[i+1][j+1] = c
                    if c > longest:
                        lcs_set = set()
                        longest = c
                        lcs_set.add(S[i-c+1:i+1])
                    elif c == longest:
                        lcs_set.add(S[i-c+1:i+1])

        return list(lcs_set)

    def _seq_dict(self, masked_reads, k):
        seq_dict = col.defaultdict()
        for f in dinopy.FastaReader(masked_reads).entries():
            seq = f.sequence.decode()
            seq_dict[seq] = {}
            seq_dict[seq]['kmers'] = []
            seq_dict[seq]['abu'] = f.name.decode().split('=')[1][:-1]
            seq_dict[seq]['name'] = f.name.decode().split(';')[0]
            for i in range(len(seq)-k+1):
                node1 = seq[i:i+k-1]
                node2 = seq[i+1:i+1+k-1]
                if node1.isupper() and node2.isupper():
                    seq_dict[seq]['kmers'].append((node1, node2))
            #for i in range(len(seq)-k+1):
             #   node = seq[i:i+k-1]
              #  if node.isupper():
               #     seq_dict[seq]['kmers'].append(node)
        seq_dict2 = {key:value for key, value in seq_dict.items() if seq_dict[key]['kmers']}
        
        return seq_dict2
        
    def _abu_kmer_zip(self, seq_dict):
        keys, values = zip(*seq_dict.items())
        abus = []
        kmers = []
        names = []
        for i in range(len(values)):
            abus.append(values[i]['abu'])
            kmers.append(values[i]['kmers'])
            names.append(values[i]['name'])
        
        return (abus, kmers, keys, names)
    
    def _shortest_common_superstring(self, intersection):
        gu = nx.DiGraph(list(intersection))
        longest_subgraph = [subgraph for subgraph in nx.weakly_connected_component_subgraphs(gu)][0]
        
        return reduce(lambda x,y: x+y[-1], list(nx.topological_sort(longest_subgraph)))
    
    def _intersection_list(self, abus, kmers, keys, names):
        intersec_list = []
        for i in range(len(kmers)):
            for j in range(len(kmers)):
                intersec = set(kmers[i]).intersection(kmers[j])
                if i != j and intersec:
                    intersec_list.append(((keys[i], abus[i], names[i]),
                                          (keys[j], abus[j], names[j]),
                                          self._shortest_common_superstring(intersec)))
        
        return intersec_list

    def _overlap_graph(self, intersec_list):
        t_g = nx.DiGraph()
        for i in range(len(intersec_list)):
            #long_com_subs = self._lcs(intersec_list[i][0][0], intersec_list[i][1][0])[0]
            scs = intersec_list[i][2]
            if int(intersec_list[i][0][1]) > int(intersec_list[i][1][1]):
                if t_g.in_edges(intersec_list[i][1][0]):
                    remove_edge = []
                    add_edge = []
                    flag = 0
                    for edge in t_g.in_edges(intersec_list[i][1][0], data=True):
                        if scs in edge[2]['seq'] and not scs == edge[2]['seq']:
                            flag = 1
                            break
                        elif edge[2]['seq'] in scs and not scs == edge[2]['seq']:
                            remove_edge.append(edge)
                        elif scs == edge[2]['seq']:
                            if int(intersec_list[i][0][1]) > int(edge[2]['abu1']):
                                remove_edge.append(edge)
                            #elif int(intersec_list[i][0][1]) == int(edge[2]['abu1']):
                                #add_edge.append((intersec_list[i][0][0], intersec_list[i][1][0], len(long_com_subs), long_com_subs, intersec_list[i][0][1], intersec_list[i][1][1]))
                            else:
                                flag = 1
                                break
                    if remove_edge:
                        for edge in remove_edge:
                            t_g.remove_edge(edge[0], edge[1])
                        t_g.add_edge(intersec_list[i][0][0], intersec_list[i][1][0], 
                                     length=len(scs), seq = scs, 
                                     abu1 = intersec_list[i][0][1], abu2 = intersec_list[i][1][1],
                                    name1 = intersec_list[i][0][2], name2 = intersec_list[i][1][2])
                    if add_edge:
                        for edge in add_edge:
                            t_g.add_edge(edge[0], edge[1], length=edge[2], seq=edge[3],abu1=edge[4],abu2=edge[5],
                                        name1 = edge[6], name2 = edge[7])
                    if flag == 0:
                        t_g.add_edge(intersec_list[i][0][0], intersec_list[i][1][0], 
                                     length=len(scs), seq = scs, 
                                     abu1 = intersec_list[i][0][1], abu2 = intersec_list[i][1][1],
                                    name1 = intersec_list[i][0][2], name2 = intersec_list[i][1][2])
                else:
                    t_g.add_edge(intersec_list[i][0][0], intersec_list[i][1][0], 
                                 length=len(scs), seq = scs, 
                                 abu1 = intersec_list[i][0][1], abu2 = intersec_list[i][1][1],
                                name1 = intersec_list[i][0][2], name2 = intersec_list[i][1][2])
            else:
                if t_g.in_edges(intersec_list[i][0][0]):
                    remove_edge = []
                    add_edge = []
                    flag = 0
                    for edge in t_g.in_edges(intersec_list[i][0][0], data=True):
                        if scs in edge[2]['seq'] and not scs == edge[2]['seq']:
                            flag = 1
                            break
                        elif edge[2]['seq'] in scs and not scs == edge[2]['seq']:
                            remove_edge.append(edge)
                        elif scs == edge[2]['seq']:
                            if int(intersec_list[i][0][1]) > int(edge[2]['abu1']):
                                remove_edge.append(edge)
                            #elif int(intersec_list[i][0][1]) == int(edge[2]['abu1']):
                                #add_edge.append((intersec_list[i][1][0], intersec_list[i][0][0], len(long_com_subs), long_com_subs, intersec_list[i][1][1], intersec_list[i][0][1]))
                            else:
                                flag = 1
                                break
                    if remove_edge:
                        for edge in remove_edge:
                            t_g.remove_edge(edge[0], edge[1])
                        t_g.add_edge(intersec_list[i][1][0], intersec_list[i][0][0], 
                                     length=len(scs), seq = scs, 
                                     abu1 = intersec_list[i][1][1], abu2 = intersec_list[i][0][1],
                                    name1 = intersec_list[i][1][2], name2 = intersec_list[i][0][2])
                    if add_edge:
                        for edge in add_edge:
                            t_g.add_edge(edge[0], edge[1], length=edge[2], seq=edge[3],abu1=edge[4],abu2=edge[5],
                                        name1=edge[6], name2=edge[7])
                    if flag == 0:
                        t_g.add_edge(intersec_list[i][1][0], intersec_list[i][0][0], 
                                     length=len(scs), seq = scs, 
                                     abu1 = intersec_list[i][1][1], abu2 = intersec_list[i][0][1],
                                    name1 = intersec_list[i][1][2], name2 = intersec_list[i][0][2])
                else:
                    t_g.add_edge(intersec_list[i][1][0], intersec_list[i][0][0], 
                                 length=len(scs), seq = scs, 
                                 abu1 = intersec_list[i][1][1], abu2 = intersec_list[i][0][1],
                                name1 = intersec_list[i][1][2], name2 = intersec_list[i][0][2])
    
        return t_g
    
    def _remove_low_indegree_edges(self, overlap_graph):
        remove_edges = []
        for edge in overlap_graph.edges():
            if overlap_graph.in_degree(edge[1]) < 2:
                remove_edges.append(edge)
        for redge in remove_edges:
            overlap_graph.remove_edge(redge[0], redge[1])
        
        return overlap_graph
    
    def subgraphs(self,overlap_graph):
        subgraph_list = [overlap_graph.subgraph(subg) for subg in nx.weakly_connected_component_subgraphs(overlap_graph)]
        high_indegree_subgraphs = [subgraph for subgraph in subgraph_list if any(x > 1 for x in dict(subgraph.in_degree()).values())]
        
        return high_indegree_subgraphs
    
    def potential_chimeras(self, subgraphs):
        node_list = []
        for subgraph in subgraphs:
            for node in subgraph.nodes():
                if subgraph.in_degree(node) > 1:
                    node_list.append((node, list(subgraph.in_edges(node, data='name2'))[0][2]))
        
        return node_list
    
    def draw_subgraphs(self):
        for i in range(len(self.chimeric_subgraphs)):
            pos = nx.spring_layout(self.chimeric_subgraphs[i], weight='length', k=5/math.sqrt(self.chimeric_subgraphs[i].order()))
            nx.draw(self.chimeric_subgraphs[i], pos)
            #node_labels = nx.get_node_attributes(subgraph_list[i], 'abundance')
            #nx.draw_networkx_labels(subgraph_list[i], pos,labels=node_labels)
            edge_labels1 = nx.get_edge_attributes(self.chimeric_subgraphs[i],'name2')
            #edge_labels2 = nx.get_edge_attributes(self.chimeric_subgraphs[i],'name1')
            nx.draw_networkx_edge_labels(self.chimeric_subgraphs[i], pos, edge_labels=edge_labels1)
            #nx.draw_networkx_edge_labels(self.chimeric_subgraphs[i], pos, edge_labels=edge_labels2)
            plt.show()
            
    def write_fasta(self, filename):
        with dinopy.FastaWriter(str(filename), 'w') as faw:
            for node in enumerate(self.potential_chimeras):
                faw.write_entry((node[1].encode(), str(node[0]).encode()))
