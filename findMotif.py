#!/usr/bin/python

import random, numpy, itertools, os, time
from pyfasta import Fasta

BASECOUNT = 4
BASE_INDEX = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }
EPSILON = [ .5, .5, .5, .5 ]
EPSILON_SUM = sum(EPSILON)

def gibbs_alogrithm(data_set):
	#get random index for random sequence
	z = random.randint(0, data_set.sequence_count - 1)

	#select this sequence
	sequence = data_set.sequence(z)

	#create PWM with every sequence but this one: use probability not frequency as per lecture slides
	counts = numpy.zeros((data_set.motif_length, BASECOUNT))
	background_counts = numpy.zeros((BASECOUNT))
	background_total = 0

	#count in frequencies for positions
	for s in range(0, data_set.sequence_count):
		if s != z:
			start = data_set.alignments[s]
			sequence = data_set.sequence(s)
			end = data_set.motif_length + start - 1
			for i in range(0, len(sequence)):
				base = sequence[i]
				j = BASE_INDEX[base]
				if i < start or i > end:
					background_counts[j] += 1
					background_total += 1
				else:
					counts[i - start, j] += 1

	# update the pwm (probabilities)
	for i in range(0, data_set.motif_length):
		for j in range(0, BASECOUNT):
			data_set.pwm[i,j] = (counts[i,j] + EPSILON[j]) / (data_set.sequence_count - 1 + EPSILON_SUM)

	background_length = 0

	for (i, (k, seq)) in enumerate(data_set.sequences.iteritems()):
		if i != z:
			background_length += len(seq) - data_set.motif_length

	# update the background frequencies
	for i in range(0, BASECOUNT):
		data_set.background_probability_matrix[i] = (background_counts[i] + EPSILON[i]) / (background_length + EPSILON_SUM)

	#evaluate alignments using obtained pwm
	max_start = len(sequence) - data_set.motif_length
	qx_over_px_list = numpy.zeros((max_start))
	qx_over_px_list_total = 0

	#aggregate total qx_over_px_list to compare against 
	for start in range(0, max_start):
		qx_over_px_list[start] = qx_over_px(data_set,sequence, start)
		qx_over_px_list_total += qx_over_px_list[start]

	# normalize the qx_over_px_list
	multiplier = 1 / qx_over_px_list_total

	for i in range(0, len(qx_over_px_list)):
		qx_over_px_list[i] *= multiplier

	#like shooting darts at a dartboard randomly choose best position: porportional qx_over_px
	totals = []
	running_total = 0

	for w in qx_over_px_list:
		running_total += w
		totals.append(running_total)

	rnd = random.random() * running_total
	for i, total in enumerate(totals):
		if rnd < total:
			data_set.alignments[z] = i

def qx_over_px(data_set, sequence, start):
	#check qx_over_px of alignment using the entire sequence
	motif_probability = 1
	background_probability = 1

	for (i, base) in enumerate(itertools.islice(sequence, start, start + data_set.motif_length)):
		base_indx = BASE_INDEX[base]
		motif_probability *= data_set.pwm[i, base_indx]
		background_probability *= data_set.background_probability_matrix[base_indx]

	return motif_probability / background_probability


#Data Container
class DataSet:
	def __init__(self, sequences, motif_length):
		self.motif_length = motif_length
		self.sequences = sequences
		self.sequence_count = len(sequences.keys())
		self.pwm = numpy.zeros((self.motif_length, BASECOUNT))
		self.background_probability_matrix = numpy.zeros(BASECOUNT)
		self.alignments = numpy.zeros(self.sequence_count)

	#access helper for pyfasta
	def sequence(self, i):
		#get keys in sequences and return ith sequence
		key = self.sequences.keys()[i]
		return self.sequences[key]


#gibbs takes dictionary with directory and iterations
def gibbs(directory, iterations):
	#pre-process arguments
	start_time = time.time()
	length_path = os.path.join(directory, 'motiflength.txt')
	sequences_path = os.path.join(directory, 'sequences.fa')

	length_file = open(length_path, 'r')
	motif_length = int(length_file.readline())

	sequences = Fasta(sequences_path)

	#create dataset
	data = DataSet(sequences, motif_length)

	#shuffle alignments prior to running gibbs
	for (i, a) in enumerate(data.alignments):
		data.alignments[i] = random.randint(0, len(data.sequence(i)) - data.motif_length - 1)

	#run gibbs for iterations iterations
	for i in range(0, iterations):
		if i%10 == 0:
			print "running iteration %i of %i" %(i,iterations)
		gibbs_alogrithm(data)

	#pring result to command line
	print "Motif:\n************\n %s" % data.pwm

	print "Sites:\n************\n %s" % data.alignments

	#
	_baseArray = ['a','c','g','t']
	for base in data.pwm:
		print _baseArray[numpy.argmax(base)]

	#write files out
	f = os.path.join(directory, "predictedsites.txt")
	f = open(f,'w')
	for site in data.alignments:
		f.write(str(int(site)+1) +"\n")

	f = os.path.join(directory, "predictedmotif.txt")
	f = open(f,'w')
	f.write(">PMOTIF\t"+str(len(data.pwm))+"\n")
	for row in data.pwm:
		for col in row:
			f.write(str(col)+"\t")
		f.write("\n")
	f.write(">")
	exec_time = time.time()-start_time
	print "took time(seconds):" +str(exec_time)

