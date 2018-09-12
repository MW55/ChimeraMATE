cimport dinopy
import dinopy
from dinopy.definitions cimport FastaEntryC
import matplotlib.pyplot as plt
import networkx as nx
from functools import reduce
import math
import itertools as it
import sys

cdef class kmer_filter:
    cdef dict __dict__
    def __init__(self, str otu_file, int k, int cutoff):
        self.reads = dinopy.FastaReader(otu_file)
        self.de_bruijn_dict = self._make_de_bruijn_file(self.reads, k)
        self.kmer_abundance_sorting = self.kmer_abundance_sorting(self.reads, self.de_bruijn_dict, cutoff)
    
    # The prior abundance of the otu should also go in the kmer abundance
    def _make_de_bruijn_file(self, reads, int k):
        cdef FastaEntryC f
        cdef str name
        cdef int i
        cdef str node
        cdef dict bru_dict
        cdef list nodes
        bru_dict = {}
        nodes = []
        for f in reads.entries():
            seq = f.sequence.decode()
            name = f.name.decode()
            for i in range(len(seq)-k+1):
                node = seq[i:i+k]
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
    
    def kmer_abundance_sorting(self, reads, dict bru_dict, float cutoff):
        cdef float num_reads
        cdef list high_abu
        cdef str kmer
        cdef float threshold
        num_reads = sum(1 for read in reads.entries())
        high_abu = []
        threshold = ((cutoff/100) * num_reads)
        for kmer in bru_dict:
             if bru_dict[kmer]['abu'] > ((cutoff/100) * num_reads): #num_reads/cutoff:
                high_abu.append(kmer)
        return high_abu
    
    def softmask(self, reads, list kmer_abundance_sorting, str output_file):
        cdef FastaEntryC entry
        cdef str seq
        cdef str kmer
        with dinopy.FastaWriter(str(output_file), line_width=1000, force_overwrite=True) as faw:
            for entry in reads.entries():
                seq = entry.sequence.decode()
                for kmer in kmer_abundance_sorting:
                    if kmer in seq:
                        seq = seq.replace(kmer, kmer.lower())
                faw.write_entry((seq.encode(), entry.name))


cdef class chimera_search:
    cdef dict __dict__
    def __init__(self, str masked_reads, int k, float abskew):
        self._seq_dict = self._seq_dict(masked_reads, k)
        self._abu_kmer_zip = self._abu_kmer_zip(self._seq_dict)
        self._intersection_list = self._intersection_list(
                self._abu_kmer_zip[1])
        self._overlap_graph = self._overlap_graph(self.create_scs_lists(self._abu_kmer_zip[0],
                                      self._abu_kmer_zip[2],
                                      self._abu_kmer_zip[3],
                                      self._intersection_list[0],
                                      self._intersection_list[1]))
        self._high_indegree_graph = self._remove_low_indegree_edges(
                self._overlap_graph, abskew)
        self.chimeric_subgraphs = self.subgraphs(self._high_indegree_graph)#_overlap_graph
        self.potential_chimeras = self.potential_chimeras(
                self.chimeric_subgraphs)
        self.write_nonchims = self.write_nonchims(masked_reads, filename = 'nonchims.fasta')

    def _seq_dict(self, str masked_reads, int k):
            cdef dict seq_dict = {}
            cdef FastaEntryC f
            cdef str seq
            cdef str node1
            cdef str node2
            cdef int i
            for f in dinopy.FastaReader(masked_reads).entries():
                seq = f.sequence.decode()
                for i in range(len(seq)-k+1):
                    node = seq[i:i+k]
                    if node.isupper():
                        if seq not in seq_dict.keys():
                            seq_dict[seq] = {}
                            seq_dict[seq]['kmers'] = []
                            seq_dict[seq]['abu'] = f.name.decode().split('=')[1][:-1]
                            seq_dict[seq]['name'] = f.name.decode().split(';')[0]
                        seq_dict[seq]['kmers'].append(node)
            return seq_dict
    
    def _abu_kmer_zip(self, dict seq_dict):
        cdef tuple keys, values
        cdef list abus
        cdef list names
        cdef int i
        keys, values = zip(*seq_dict.items())
        abus = []
        kmers = []
        names = []
        for i in range(len(values)):
            abus.append(values[i]['abu'])
            kmers.append(values[i]['kmers'])
            names.append(values[i]['name'])

        return (abus, kmers, keys, names)
    
        
    ## this works, but keep in mind that the resulting longest common substring does not neccessarily exists 
    # in the reduced form in any of the target sequences, as there could be a situation in which
    # both sequences share an overlap and then a softmasked region, followed by another overlap
    # the reduce function ignores this "gap" and concatenates the kmers anyway
    def _longest_common_substring(self, tuple pair, list kmers, list kmer_sets):
        cdef list tranverse
        cdef list search
        cdef list s_search
        cdef list longest_entry
        cdef list buffer
        cdef int i
        cdef int j
        if len(kmers[pair[0]]) > len(kmers[pair[1]]):
            tranverse = self._abu_kmer_zip[1][pair[0]]
            search = self._abu_kmer_zip[1][pair[1]]
            s_search = sorted(self._abu_kmer_zip[1][pair[1]])
            k_set = kmer_sets[pair[1]]
        else:
            tranverse = self._abu_kmer_zip[1][pair[1]]
            search = self._abu_kmer_zip[1][pair[0]]
            s_search = sorted(self._abu_kmer_zip[1][pair[0]])
            k_set = kmer_sets[pair[0]]
        longest_entry = []
        buffer = []
        for i in range(len(tranverse)):
            if not buffer:
                if tranverse[i] in k_set:
                    buffer.append(tranverse[i])
                    j = search.index(tranverse[i])
            elif len(search) >= j+2 and tranverse[i] == search[j+1]:
                buffer.append(tranverse[i])
                j += 1
            else:
                if len(buffer) > len(longest_entry):
                    longest_entry = buffer.copy()
                    buffer = []
        if len(buffer) > len(longest_entry):
            longest_entry = buffer.copy()
        return reduce(lambda x,y: x+y[-1], longest_entry)
    
    #add static typisation
    def _intersection_list(self, list kmers):
        kmer_sets = [set(kmer_l) for kmer_l in kmers]
        pairs = it.combinations(range(len(kmers)), 2)
        res = []
        for tup in pairs:
            if not kmer_sets[tup[0]].isdisjoint(kmer_sets[tup[1]]):
                res.append(tup)
        return res, kmer_sets

    def create_scs_lists(self, list abus, tuple keys, list names, list nums, list kmer_sets): #tuple seqs
        cdef list intersec_list
        intersec_list = [((keys[nums[i][0]], abus[nums[i][0]],
            names[nums[i][0]]), (keys[nums[i][1]], abus[nums[i][1]],
                names[nums[i][1]]), self._longest_common_substring(nums[i], self._abu_kmer_zip[1], kmer_sets))
            for i in range(len(nums))]
    
        return intersec_list

    def _compare_sequences(self, list intersec_list, graph, int i,
            tuple direction, str scs):
            cdef list remove_edge
            cdef int flag
            cdef tuple edge
            if graph.in_edges(intersec_list[i][direction[1]][0]):
                remove_edge = []
                flag = 0
                for edge in graph.in_edges(intersec_list[i][direction[1]][0],
                        data=True):
                    if scs in edge[2]['seq'] and not scs == edge[2]['seq']:
                        flag = 1
                        break
                    elif edge[2]['seq'] in scs and not scs == edge[2]['seq']:
                        remove_edge.append(edge)
                    elif scs == edge[2]['seq']:
                        if int(intersec_list[i][direction[0]][1]) > int(
                                edge[2]['abu1']):
                            remove_edge.append(edge)
                        else:
                            flag = 1
                            break
                if remove_edge:
                    for edge in remove_edge:
                        graph.remove_edge(edge[0], edge[1])
                    self._add_edge(intersec_list, graph, i, direction, scs)
                if flag == 0:
                    self._add_edge(intersec_list, graph, i, direction, scs)
            else:
                self._add_edge(intersec_list, graph, i, direction, scs)

    def _add_edge(self, list intersec_list, graph, int i,tuple direction,
            str scs):
        graph.add_edge(intersec_list[i][direction[0]][0],
                intersec_list[i][direction[1]][0], length=len(scs), seq = scs,
                abu1 = intersec_list[i][direction[0]][1],
                abu2 = intersec_list[i][direction[1]][1],
                name1 = intersec_list[i][direction[0]][2],
                name2 = intersec_list[i][direction[1]][2])

    def _overlap_graph(self, list intersec_list):
        cdef int i
        cdef str scs
        t_g = nx.DiGraph()
        for i in range(len(intersec_list)):
            scs = intersec_list[i][2]
            if int(intersec_list[i][0][1]) > int(intersec_list[i][1][1]):
                self._compare_sequences(intersec_list, t_g, i, (0, 1), scs)
            else:
                self._compare_sequences(intersec_list, t_g, i, (1, 0), scs)

        return t_g

    def _remove_low_indegree_edges(self, overlap_graph, float abskew):
        cdef list remove_edges
        remove_edges = list()
        for edge in overlap_graph.edges(data=True):
            if overlap_graph.in_degree(edge[1]) < 2:
                remove_edges.append(edge)
            if int(edge[2]['abu2'])/int(edge[2]['abu1']) >= abskew:
                remove_edges.append(edge)
        for redge in list(remove_edges):
            try:
                overlap_graph.remove_edge(redge[0], redge[1])
            except:
                pass
        
        return overlap_graph

    def subgraphs(self,overlap_graph):
        cdef list subgraph_list
        cdef list high_indegree_subgraphs
        subgraph_list = [overlap_graph.subgraph(subg) for subg
                in nx.weakly_connected_component_subgraphs(overlap_graph)]
        high_indegree_subgraphs = [subgraph for subgraph in subgraph_list
                if any(x > 1 for x in dict(subgraph.in_degree()).values())]

        return high_indegree_subgraphs

    def potential_chimeras(self, list subgraphs):
        cdef list node_list
        node_list = []
        for subgraph in subgraphs:
            for node in subgraph.nodes():
                if subgraph.in_degree(node) > 1:
                    node_list.append((node,
                        list(subgraph.in_edges(node, data='name2'))[0][2]))

        return node_list

    def draw_subgraphs(self, list chimeric_subgraphs):
        cdef int i
        cdef str weight
        cdef int k
        for i in range(len(chimeric_subgraphs)):
            pos = nx.spring_layout(chimeric_subgraphs[i], weight='length',
                    k=5/math.sqrt(chimeric_subgraphs[i].order()))
            nx.draw(chimeric_subgraphs[i], pos)
            edge_labels1 = nx.get_edge_attributes(chimeric_subgraphs[i],'name2')
            nx.draw_networkx_edge_labels(chimeric_subgraphs[i], pos,
                    edge_labels=edge_labels1)
            plt.show()
    
    def write_chims(self, str filename):
        with dinopy.FastaWriter(filename, 'w') as faw:
            for node in enumerate(self.potential_chimeras):
                faw.write_entry((node[1][0].encode(),
                    '{};{}'.format(str(node[0]), node[1][1]).encode()))
    
    def write_nonchims(self, masked_reads, str filename):
        with dinopy.FastaWriter(filename, 'w', force_overwrite=True) as faw: 
            msk = dinopy.FastaReader(masked_reads, 'r')
            for read in msk.entries():
                seq = read.sequence.decode()
                if not read.name.decode().split(';')[0] in set(list(zip(*self.potential_chimeras))[1]):
                    faw.write_entry((seq.upper().encode(), read.name))
