import warnings
warnings.filterwarnings("ignore")
import Tkinter, tkFileDialog
from difflib import SequenceMatcher
import random
import os
import csv
import datetime

# represents a probe sequence
class Probe():

	def __init__(self, seq, percent_gc):

		self.seq = seq
		self.percent_gc = percent_gc
		self.probe_sim = None

# class to represent a specific gene for which a probe will be created
class Gene():

	def __init__(self, name, seq, num_of_probes):

		if seq < 18:

			raise Exception("sequence too small to build probe for")

		self.name = name
		self.seq = seq
		self.num_of_probes = num_of_probes
		self.all_probes = []
		self.accepted_probes = []
		self.current_probe = ""
		self.current_probe_similarities = []
		self.final_probe = ""
		self.log = ""

	# creates all possible seqs from 18-50 bp with the given seq
	def get_all_probes(self):

		self.log += "Getting all probes\n"
			
		for probe_size in range(18, self.get_max_probe_len() + 1):

			for i in range(0, len(self.seq) - probe_size + 1):

				if len(self.all_probes) < self.num_of_probes:

					current_probe = self.seq[i:i + probe_size]

					self.log += "	probe ({0}): {1}\n".format(probe_size, current_probe)

					self.all_probes.append(current_probe)

				else:

					break

	# gets the maximum possible probe length for the given seq
	def get_max_probe_len(self):

		if len(self.seq) > 50:

			return 50

		else:

			return len(self.seq)

	# checks all probes for the following criteria and places only the good probes
	# into the self.accepted_probes
	def qualify_probes(self):

		self.log += "Getting probes that follow criteria\n"

		for probe in self.all_probes:

			if self.get_gc_content(probe) >= 40.0 and self.get_gc_content(probe) <= 60.0 and self.passes_complement(probe) and self.passes_repeats(probe):

				accepted_probe = Probe(probe, self.get_gc_content(probe))

				self.log += "	probe ({0}): {1} | {2} percent gc\n".format(len(accepted_probe.seq), accepted_probe.seq, accepted_probe.percent_gc)

				self.accepted_probes.append(accepted_probe)

	# checks if probe is 40-60% GC
	def get_gc_content(self, probe):

		num_gc = 0

		for base in probe:


			if base == "G" or base == "C":

				num_gc += 1

		percent_gc = ((1.0 * num_gc) / (1.0 * len(probe)) * 100)

		return percent_gc

	# checks if probe does not contain complementary sequences
	def passes_complement(self, probe):

		return True

	# check is probe does not contain repeats of a base of more than 4
	def passes_repeats(self, probe):

		if "AAAAA" in probe:
			return False
		elif "GGGGG" in probe:
			return False
		elif "CCCCC" in probe:
			return False
		elif "TTTTT" in probe:
			return False
		else:
			return True

	# generates all and the qualified probes
	def get_probes(self):

		self.get_all_probes()
		self.qualify_probes()


"""
gathers the seqs as Genes and creates probes iteratively for the whole set with the
following criteria:
	18-50 bp in length
	40-60 percent GC
	no complementary region within the probe
	no stretches of more than 4 of the same base
	no similarity to any region in the genome of more than 70% (not handled here)
"""
class ProbeFinder():

	def __init__(self, probe_num):

		self.genes = []
		self.gene_probe_similarity = [ [] ] # [gene [probe sim % ] ]
		self.probe_num = probe_num
		self.log = ""

	# creates genes from input from fasta file
	def genes_from_fasta(self):

		root = Tkinter.Tk()
		root.withdraw()
		file_path = tkFileDialog.askopenfilename()

		with open(file_path, 'r') as seq:

			data = seq.read()
			# read the string and split it by '>' character
			indiv_seqs = data.split('>')

			for seq in indiv_seqs:

				if len(seq) > 1:

					data = seq.split('\n')
					name = data[0]
					sequence = data[1]

					new_gene = Gene(name, sequence, self.probe_num)

					self.genes.append(new_gene)

					self.log += "Gene: {0}\n".format(name)

	# generates probes for all genes
	def generate_probes_for_all_genes(self):

		for gene in self.genes:

			self.log += "generating probes for {0}".format(gene.name)
			print "generating probes for {0}".format(gene.name)

			gene.get_probes()

			print "{0} probes created".format(len(gene.accepted_probes))

		# remove genes with no probes created
		for gene in self.genes:

			if len(gene.accepted_probes) == 0:

				self.genes.remove(gene)

	# finds unique probes for all of the genes that < 70% similar from any
	# of the other probes
	def find_unique_probes(self):

		probes_pass = False
		i = 0

		while not probes_pass:

			print "Unique iteration: {0}".format(i)

			# holds index of current probe for each gene in the program
			active_probes = []

			for gene in self.genes:

				try:

					# random pick
					len_of_probes = len(gene.accepted_probes)
					probe_idx = random.randint(0, len_of_probes - 1)
					active_probes.append(probe_idx)

				except:

					print gene.name
					print gene.seq

			if self.probes_pass(active_probes):
			
				# return the probes at the given index for each gene
				self.return_probes(active_probes)

				probes_pass = True

			i += 1

	# returns the probes for each gene at the given index
	def return_probes(self, active_probes):

		for index, gene in enumerate(self.genes):

			gene.final_probe = gene.accepted_probes[active_probes[index]]

			print "Probe for gene {0} is {1} | {2} percent sim".format(gene.name, gene.final_probe.seq, gene.final_probe.probe_sim)

	# tests if all probe relations in the active probe array are less than
	# 70%
	def probes_pass(self, active_probes):
		
		pass_percentages = []
		probe_seqs = []

		# get the seqs
		for index, probe_idx in enumerate(active_probes):

			probe_seqs.append(self.genes[index].accepted_probes[probe_idx])

		# calculate relatedness
		for probe in probe_seqs:

			probe_sims = []

			other_probes = list(probe_seqs)
			other_probes.remove(probe)

			for other_probe in other_probes:

				relatedness = SequenceMatcher(None, probe.seq, other_probe.seq).ratio()

				probe_sims.append(relatedness)

			pass_percentages.append(probe_sims)

			# add the average probe_sim % to the probe
			probe_sim_total = 0.0
			avg_probe_sim = 0.0
			for sim in probe_sims:
				probe_sim_total += sim
			avg_probe_sim = ((1.0 * probe_sim_total) / len(probe_sims)) * 100
			probe.probe_sim = avg_probe_sim


		all_less_than_70 = True

		for probe_results in pass_percentages:

			for sim_percentage in probe_results:

				if sim_percentage > 70.0:

					all_less_than_70 = False

		return all_less_than_70

	# writes probe sequences to file
	def final_data_to_file(self):

		print "Writing probe results to file ..."

		save_file = one_dir_back(os.getcwd()) + "/ProbeSummary.csv"

		with open(save_file, 'wb') as f:

			c = csv.writer(f, delimiter = ',', quotechar = '|')

			for index, gene in enumerate(self.genes):

				organism = find_between(gene.name, "[", "]")
				gene_name = find_between(gene.name, "|", "[")
				gene_name = gene_name[::-1][:10][::-1]
				gene_name = find_between(gene_name, "|", " ")
				gene_name = gene_name.replace("|", "")

				if index == 0:
					c.writerow(["Name", "Gene Name", "Organism", "Probe Seq", "Avg Percent Sim", "Probe Len (bp)"])
					c.writerow([gene.name, gene_name, organism, gene.final_probe.seq, gene.final_probe.probe_sim, len(gene.final_probe.seq)])
				
				else:
					c.writerow([gene.name, gene_name, organism, gene.final_probe.seq, gene.final_probe.probe_sim, len(gene.final_probe.seq)])


# gets the name of the organism the gene belongs to
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

def find_between_r( s, first, last ):
    try:
        start = s.rindex( first ) + len( first )
        end = s.rindex( last, start )
        return s[start:end]
    except ValueError:
        return ""

# returns the pwd, minus one level of depth
def one_dir_back(current_directory):

	rev_dir = current_directory[::-1]

	rev_result = ''

	result = ''

	for c in rev_dir:

		if c == '/':

			rev_result += rev_dir[rev_dir.index(c):]
			
			result = rev_result[::-1]

			return result

def main():

	probe_finder = ProbeFinder(1000)

	probe_finder.genes_from_fasta()
	probe_finder.log += "creating probe sets for all genes"
	probe_finder.generate_probes_for_all_genes()
	probe_finder.find_unique_probes()
	probe_finder.final_data_to_file()

if __name__ == "__main__":

	main()






















