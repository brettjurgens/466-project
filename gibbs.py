import sys
from pyfasta import Fasta
import itertools

# def produceShifts(sequence_lengths):
# 	shift_possibilities= []
# 	#create individual arrays of shifts
# 	for sequence_length in sequence_lengths:
# 		position = []
# 		for i in sequence_lengths:
# 			position.append(i)
# 			shift_positions.append(position)
# 		for position in positions:

# 	return

	#combine

def product(*args):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

def brutalPossibilities(sequences, motiflength):

	#generate possible slices
	possibile_alignments = []
	for sequence in sequences:
		# print "new sequence"
		alignments = []
		for i in range(len(sequence) - motiflength):
			alignments.append(sequence[i:i+motiflength])
			 # print alignments
			possibile_alignments.append(alignments)
		# print possibile_alignments
	# print possibile_alignments

	print "starting itertools"

	possibile_alignments = list(product(*possibile_alignments))
	# print possibile_alignments
	#score motifs
	motifPossibilities = brutalMotifPossibilities(motiflength, "")

	max_score = 0
	max_motif = ""
	motif_num = 1
	for motif in motifPossibilities:
		print "running motif " + str(motif_num) + " out of " + str(len(motifPossibilities)) #TIME NIGGA!

		# print motif
		score = 0		
		for alignment in possibile_alignments:
			# print "running" + str(alignment) + "on" + motif
			for i in range(motiflength):
				for align in alignment:
					if align[i] == motif[i]:
						score += 1
			# print score
		
		motif_num += 1

		if score > max_score:	
			max_score = score
			max_motif = motif

	return max_motif + "with score" + str(score)


def brutalMotifPossibilities(motiflength, current):
	#create all possibilities for motifs
	# print "starting motif"
	if motiflength == 0:
		return [current]
	# print current
	return brutalMotifPossibilities(motiflength-1, current+"A") + brutalMotifPossibilities(motiflength-1, current+"C") + brutalMotifPossibilities(motiflength-1, current+"G") + brutalMotifPossibilities(motiflength-1, current+"T")


# def brutal(sequences, motiflength):
# 	possibilities = brutalPossibilities(motiflength, "")
# 	for motif in possibilities:
# 		#create and score PWM

def main():
	if len(sys.argv) != 3:
		print "usage: python gibbs.py sequences.fa motiflength.txt"

		#debug code:
		# print brutalPossibilities(["ttgttccggagtggtttactattgtctggt", "tccaaactatttagcattttacgatcgcat", "cttagttctgtacaatgatcgaaaattcat"], 3)

	#take in input
	sequence_fa = sys.argv[1]
	motiflength_txt = sys.argv[2]

	# sequences = Fasta(sequence_fa)
	sequences = []
	sequences_file = open(sequence_fa, 'r')
	line = "yoyoyo"
	while line != "":
		line = sequences_file.readline() 
		sequences.append(sequences_file.readline()[:-1])

	print sequences

	motiflength_file = open(motiflength_txt, 'r')
	motif_length = int(motiflength_file.readline())

	#process input
	print brutalPossibilities(sequences, motif_length)
	return 

main()

