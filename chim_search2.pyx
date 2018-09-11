cimport dinopy
import dinopy
import collections as col
from dinopy.definitions cimport FastaEntryC
import matplotlib.pyplot as plt
import networkx as nx
from functools import reduce
import math
import itertools as it
import multiprocessing as mp
import networkx as nx
import os
import yaml

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
    def __init__(self, str masked_reads, int k, int threads, float abskew):
        self._seq_dict = self._seq_dict(masked_reads, k)
        self._abu_kmer_zip = self._abu_kmer_zip(self._seq_dict)
        self._intersection_list = self._intersection_list(
                self._abu_kmer_zip[1])
        self._overlap_graph = self._overlap_graph(self.run_mp(threads))
        self._high_indegree_graph = self._remove_low_indegree_edges(
                self._overlap_graph, abskew)
        self.chimeric_subgraphs = self.subgraphs(self._overlap_graph)
        self.potential_chimeras = self.potential_chimeras(
                self.chimeric_subgraphs)

    #switch it up so that it only creates the seq_dict entry if both nodes a uppercase
    def _seq_dict(self, str masked_reads, int k):
            cdef dict seq_dict = {}
            cdef FastaEntryC f
            cdef str seq
            cdef str node1
            cdef str node2
            cdef dict seq_dict2
            cdef int i
            for f in dinopy.FastaReader(masked_reads).entries():
                seq = f.sequence.decode()
                seq_dict[seq] = {}
                seq_dict[seq]['kmers'] = []
                seq_dict[seq]['abu'] = f.name.decode().split('=')[1][:-1]
                seq_dict[seq]['name'] = f.name.decode().split(';')[0]
                for i in range(len(seq)-k+1):
                    node1 = seq[i:i+k-1]
                    node2 = seq[i+1:i+k]
                    if node1.isupper() and node2.isupper():
                        seq_dict[seq]['kmers'].append((node1, node2)) #here should come the declarations of seq, kmers,abu,name
            seq_dict2 = {key:value for key, value in seq_dict.items()
                    if seq_dict[key]['kmers']}

            return seq_dict2
    
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
    

    def _shortest_common_superstring(self, intersection):
        cdef str x
        cdef str y
        gu = nx.DiGraph(list(intersection))
        longest_subgraph = max(nx.weakly_connected_component_subgraphs(gu),
                key=len)

        return reduce(lambda x,y: x+y[-1],
                list(nx.topological_sort(longest_subgraph)))
    
    #add static typisation
    def _intersection_list(self, list kmers):
        kmer_sets = [set(kmer_l) for kmer_l in kmers]
        pairs = it.combinations(range(len(kmer_sets)), 2)
        #intersec = lambda a, b: kmer_sets[a].intersection(kmer_sets[b])
        res = []
        for tup in pairs:
            set_intersec = kmer_sets[tup[0]].intersection(kmer_sets[tup[1]])
            if set_intersec:
                res.append((tup, set_intersec))
        nums, seqs = zip(*res)
        return nums, seqs

    def create_scs_lists(self, abus, keys, names, nums, seqs, output, seqrange):
        intersec_list = [((keys[nums[i][0]], abus[nums[i][0]],
            names[nums[i][0]]), (keys[nums[i][1]], abus[nums[i][1]],
                names[nums[i][1]]), self._shortest_common_superstring(seqs[i]))
            for i in seqrange]
    
        output.put(intersec_list)
    
    def _calc_ranges(self, seqs, threads):
        seq_count = len(seqs)
        #threads = 4
        spread = seq_count / threads
        length_of_range = int(spread)
        list_of_ranges = []
        for i in range(0, seq_count+1, length_of_range):
            if i > 0:
                list_of_ranges.append(range(last_i, i))
                #if type(spread) == float:
                 #   i += 2
                  #  list_of_ranges.append(range(last_i, i))
                #else:
                 #   list_of_ranges.append(range(last_i, i))
            last_i = i
        return list_of_ranges
    
    def run_mp(self, threads):
        nums, seqs = self._intersection_list
        mananger = mp.Manager()
        output = mananger.Queue()
        #output = mp.Queue()

        # Setup a list of processes that we want to run
        #processes = [mp.Process(target=rand_string, args=(5, output)) for x in range(4)]
        processes = [mp.Process(target=self.create_scs_lists, 
                                args=(self._abu_kmer_zip[0],
                                      self._abu_kmer_zip[2],
                                      self._abu_kmer_zip[3], 
                                      nums, seqs, output, 
                                      self._calc_ranges(seqs, threads)[x]))
                     for x in range(4)]


        # Run processes
        for p in processes:
            p.start()

        # Exit the completed processes
        for p in processes:
            p.join()

        # Get process results from the output queue
        return [item for sublist in [output.get() for p in processes] for item in sublist]

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

    def _remove_low_indegree_edges(self, overlap_graph, abskew):
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
    
    # still have to write out the sequences that are NOT chimeras, thats the whole point you imbecile
    def write_fasta(self, str filename):
        with dinopy.FastaWriter(filename, 'w') as faw:
            for node in enumerate(self.potential_chimeras):
                faw.write_entry((node[1][0].encode(),
                    '{};{}'.format(str(node[0]), node[1][1]).encode()))


def try_config(mask_kmers_range, otus, softmask_file, output_file):
    results_dict = {}
    for mask_kmers in mask_kmers_range: #range(12, 36)
            print('masked_kmer number {}'.format(mask_kmers))
            results_dict['mask_kmers_{}'.format(mask_kmers)] = {}
            for cutoff in range(1):
                print('cutoff {} of 50'.format(cutoff))
                results_dict['mask_kmers_{}'.format(mask_kmers)]['cutoff_{}'.format(cutoff)] = {}
                g = kmer_filter(otus, mask_kmers, cutoff)
                g.softmask(g.reads, g.kmer_abundance_sorting, softmask_file)
                for chim_kmers in [38]:
                    print('chim_kmer number {} of 36'.format(chim_kmers))
                    results_dict['mask_kmers_{}'.format(mask_kmers)][
                        'cutoff_{}'.format(cutoff)]['chim_kmers_{}'.format(chim_kmers)] = {}
                    for abskew in [0.1]:
                        print('abskew {}'.format(abskew))
                        results_dict['mask_kmers_{}'.format(mask_kmers)][
                        'cutoff_{}'.format(cutoff)]['chim_kmers_{}'.format(chim_kmers)]['abskew_{}'.format(abskew)] = {}
                        try:
                            c = chimera_search(softmask_file, chim_kmers, 4, abskew)
                        except:
                            print('Error')
                            #os.remove(softmask_file)
                            continue
                        try:
                            seq, names = zip(*c.potential_chimeras)
                            co = 0
                            for name in names:
                                if 'ch' in name:
                                    co += 1   
                            results_dict['mask_kmers_{}'.format(mask_kmers)][
                                'cutoff_{}'.format(cutoff)]['chim_kmers_{}'.format(chim_kmers)][
                                'abskew_{}'.format(abskew)][
                                'num_pot_chims'] = len(c.potential_chimeras)
                            results_dict['mask_kmers_{}'.format(mask_kmers)][
                                'cutoff_{}'.format(cutoff)]['chim_kmers_{}'.format(chim_kmers)][
                                'abskew_{}'.format(abskew)][
                                'num_real_chims'] = co
                            results_dict['mask_kmers_{}'.format(mask_kmers)][
                                'cutoff_{}'.format(cutoff)]['chim_kmers_{}'.format(chim_kmers)][
                                'abskew_{}'.format(abskew)][
                                'false_pos'] = len(c.potential_chimeras) - co
                            with open(str(output_file), 'w') as yaml_file:
                                yaml.dump(results_dict, yaml_file, default_flow_style=False)
                        except:
                            print('Error writing output')
                            continue
                os.remove(softmask_file)

def run():
    otus = 'shuff.fasta'
    try_config(38, otus, 'masked.fasta', 'test.yaml')