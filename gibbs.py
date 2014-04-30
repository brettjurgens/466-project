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
def brutalPossibilities(sequences, motiflength):

	#generate possible slices
	possibile_alignments = []
	for sequence in sequences:
		alignments = []
		for i in range(len(sequence) - motiflength):
			alignments.append(sequence[i:i+motiflength])
		possibile_alignments.append(alignments)
	possibile_alignments = list(itertools.product(*possibile_alignments))
	# print possibile_alignments
	#score motifs
	motifPossibilities = brutalMotifPossibilities(motiflength, "")

	max_score = 0
	max_motif = ""
	motif_num = 1
	for motif in motifPossibilities:
		# print motif
		score = 0		
		for alignment in possibile_alignments:
			# print "running" + str(alignment) + "on" + motif
			for i in range(motiflength):
				for align in alignment:
					if align[i] == motif[i]:
						score += 1
			# print score
		print "ran motif " + str(motif_num) + " out of " + str(len(motifPossibilities))
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

	sequences = Fasta(sequence_fa)

	motiflength_file = open(motiflength_txt, 'r')
	motif_length = int(motiflength_file.readline())

	#process input
	print brutalPossibilities(sequences, motif_length)
	return 

main()

