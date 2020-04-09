import collections, sys, string, re
from collections import OrderedDict
import string


# Initialization of setups for building the graph
class de_bruijn_vertex:
    def __init__(self, unique_kmer_seq, full_seq):
        self.outedges = []
        self.inedges = []
        self.unique_kmer = unique_kmer_seq
        self.mark = 0
        self.sequence = full_seq
        self.unique_kmer = unique_kmer_seq
        self.branching = 0
        self.num_repeats = 0
        self.begin = 0
        self.end = 0

# Build de Bruijin graph to filter inputs to get rid of nodes that occur only once
class debruijin:
    def __init__(self):
        self.debruijin = OrderedDict()


# Add in unique k-mers we defined into our nodes
    def add(self, previous_node, current_node, next_node, reads_node):
        new_node = de_bruijn_vertex(current_node, reads_node)
        
        # If the next node is available, we will add it and edges to graph
        if next_node == 0:
            if current_node not in self.debruijin:
                self.debruijin[current_node] = new_node
                
                # count number of times that the node's string
                # appears in the reads you are processing via num_repeats.
                self.debruijin[current_node].num_repeats = 1
            else:
                # Check whether we re-visited the node or not
                self.debruijin[current_node].num_repeats += 1

            # Check if the previous node is in-degree
            if previous_node != 0:
                if not previous_node in self.debruijin[current_node].inedges:
                    self.debruijin[current_node].inedges.append(previous_node)
                    
        # Similarly, if it's available, add it in the graph
        elif previous_node == 0:
            if current_node not in self.debruijin:
                self.debruijin[current_node] = new_node
                self.debruijin[current_node].num_repeats = 1
            else:
                self.debruijin[current_node].num_repeats += 1

            if not next_node in self.debruijin[current_node].outedges:
                # Check if the edge is out-degree, if it is available, add it in the graph
                self.debruijin[current_node].outedges.append(next_node)


        # Similar process if above condition is false, keep adding in nodes
        # into the graph if conditions are met
        else:
            if current_node not in self.debruijin:
                self.debruijin[current_node] = new_node
                self.debruijin[current_node].num_repeats = 1
            else:
                self.debruijin[current_node].num_repeats += 1
            if not previous_node in self.debruijin[current_node].inedges:
                self.debruijin[current_node].inedges.append(previous_node)
            if not next_node in self.debruijin[current_node].outedges:
                self.debruijin[current_node].outedges.append(next_node)


# Add a string of arbitrary length to the graph to identify good sequences
class good_debruijin:
    def __init__(self):
        self.debruijin = OrderedDict()
        
    # Create an add function to add sequence with the length into the graph
    def add(self, sequence):
        if not sequence in self.debruijin:
            new_node = seq_reads_node(sequence)
            self.debruijin[sequence] = new_node

    # Remove sequences that are in the graph
    def remove(self, sequence):
        if sequence in self.debruijin:
            del self.debruijin[sequence]

   

# Start getting to read in sequences as nodes in graph
class seq_reads_node:
    def __init__(self, sequence):
        self.reads_node = sequence

# Read in the files and make individual lines into unique k-mers' we want
def files_to_reads_node(file_for_read):
   
    file = open(file_for_read, 'r')
    for line in file:
        line = str.replace(line, '\n', '')
        length = len(line)
        previous_node = 0;
        good_debruijin.add(line)

        # Range check to make sure that kmers are not out of bound
        for i in range(0, length):
            if i + kmer_size <= length:  
                if i + kmer_size + 1 <= length: 
                    next_node = line[i + 1: i + kmer_size + 1]
                else:  
                    next_node = 0
                current_node = line[i: (kmer_size + i)]
                debruijin.add(previous_node, current_node, next_node, line)
                previous_node = current_node
    file.close()

    # Define branching node as either indegree or outdegree > 1
    for key, num_occur in debruijin.debruijin.items():
        if len(num_occur.inedges) > 1 or len(num_occur.outedges) > 1:
            num_occur.branching = 1
        if len(num_occur.inedges) == 0:
            num_occur.begin = 1
        if len(num_occur.outedges) == 0:
            num_occur.end = 1


# Vertices that appear only once to be erroneous so we will remove those
# Iterate reads and build de Bruijn graph to keep good reads only
def goodreads_nodes(debruijin, good_debruijin):
    for key, num_occur in debruijin.debruijin.items():

        if num_occur.num_repeats == 1:
            good_debruijin.remove(num_occur.sequence)

        # We only keep good reads that occur more than once
        if len(num_occur.outedges) > 1 or len(num_occur.inedges) > 1:
            num_occur.branching = 1


    # Write to file those good reads only
    f = open('good_reads', 'w')
    for key, num_occur in good_debruijin.debruijin.items():
        f.write(num_occur.reads_node + '\n')
    f.close
    debruijin.debruijin.clear()


# rebuild the de Bruijn graph, this time using only the good reads
def rebuild_debruijin():
    rebuild_file = 'good_reads'
    files_to_reads_node(rebuild_file)


# We will assemble possible sequence by outputs of contigs corresponding to unambiguous
# regions of the graph.
def assemble(debruijin):
    #initialization and start from source node
    contig = []

    for key, num_occur in debruijin.debruijin.items():
        if num_occur.begin == 1:
            result = result_contigs(num_occur, debruijin)
            if len(result) > 1:
                contig.append(result)

    # Walk through the graph, marking node as we visit it,
    # until we reach either a branching node or marked node by counting
    # how many numbers we have visited that node
    for key, num_occur in debruijin.debruijin.items():
        if num_occur.mark == 0:

            # check indegree 1 whose incoming edge comes from a branching node.
            if num_occur.branching == 1:
                continue
            if debruijin.debruijin[num_occur.inedges[0]].branching == 1:
                result = result_contigs(num_occur, debruijin)
                if len(result) > 1:
                    contig.append(result)
    return contig

# We want unambiguous contigs by connecting prefix from previous node to the
# last character of the suffix of the next node except for branching node.
def result_contigs(num_occur, debruijin):
    result = num_occur.unique_kmer
    num_occur.mark = 1
    if num_occur.end == 1:
        return result
    num_occur = debruijin.debruijin[num_occur.outedges[0]]

    # Looking for branching node to remove and then connect prefix and suffix
    while (1):
        if num_occur.branching == 1:
            num_occur.mark = 1
            result += num_occur.unique_kmer[-1:]
            break
        if num_occur.mark == 1:
            break
        result += num_occur.unique_kmer[-1:] #suffix

        # We keep marking visited nodes and move to next edge when we are done
        num_occur.mark = 1
        if num_occur.end == 0:
            num_occur = debruijin.debruijin[num_occur.outedges[0]]
        else:
            break
    return result


# Write results to files
def output_contig(contigs):
    contig_file_full = open('output_contigs_full', 'w')
    contig_length_full = open('contig_lengths_full', 'w')
    contig_file = open('output_contigs', 'w')
    contig_length = open('contig_lengths', 'w')
    for num_occur in contigs:
        contig_file_full.write(num_occur + '\n')
        contig_length_full.write(str(len(num_occur)) + '\n')

        # Find contigs with the size over 100bp and output the results
        if len(num_occur) >100:
            contig_length.write(str(len(num_occur)) + '\n')
        contig_file.write(num_occur + '\n')

    contig_file_full.close
    contig_length_full.close
    contig_file.close
    contig_length.close


debruijin = debruijin()
good_debruijin = good_debruijin()


if __name__ == "__main__":
    # Change the file directory here
    all_sequences = '/Users/user/Downloads/hw3/sequence_reads.txt'

    # get input from user. The default kmer is 31
    kmer_size = int(input("Please enter the k-mers you want to analyze: "))
    
    files_to_reads_node(all_sequences)
    goodreads_nodes(debruijin, good_debruijin)
    rebuild_debruijin()
    output_contig(assemble(debruijin))
