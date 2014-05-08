#!/usr/bin/python

import time, numpy, os, csv, math, sys
from findMotif import gibbs

def main():
    # takes in path and iterations
    if(len(sys.argv) < 3):
        print "Usage ./evaluate.py data_set iterations"
    directory = sys.argv[1]
    iterations = int(sys.argv[2])
    evaluate_dir(directory, iterations)

def evaluate_all():
    """
    evaluates the motif-finder on the benchmarks
    """
    count = 0
    for root, dirs, files in os.walk('data'):
        for directory in dirs:
            evaluate_dir(directory, 10)

def evaluate_dir(directory, iterations):
    path = "data/{0}".format(directory)

    print "Processing folder with {0} iterations".format(iterations)

    start_time = time.time()
    # run the gibbs
    gibbs(path, iterations)
    run_time = time.time() - start_time

    print "That took {0} sec!".format(run_time)

    motif_file = open(os.path.join(path, 'motif.txt'), 'r')
    motif = motif_file.readline().split('\t')[2]
    motif_file.close()

    motif_length = len(motif)

    motif_pwm = motif_to_pwm(motif)

    # get the predicted motif
    predicted_motif = []
    p_motif_file = open(os.path.join(path, 'predictedmotif.txt'), 'r')

    # skip the >PMOTIF line
    p_motif_file.readline()

    # read it into memory...
    reader = csv.reader(p_motif_file, delimiter='\t')
    predicted_motif = list(reader)

    # close the file...
    p_motif_file.close()

    # remove the > line
    predicted_motif.pop()

    # remove blank elements from the 2d array
    for arr in predicted_motif:
        arr.pop()

    relative_entropy = compute_relative_entropy(motif_pwm, predicted_motif)

    print relative_entropy

    # open sites
    sites_file = open(os.path.join(path, 'sites.txt'), 'r')
    reader = csv.reader(sites_file, delimiter='\n')
    sites = list(reader)
    sites = numpy.array(sites)
    sites = sites.flatten().tolist()
    sites_file.close()

    # open predicted sites
    p_sites_file = open(os.path.join(path, 'predictedsites.txt'), 'r')
    reader = csv.reader(p_sites_file, delimiter='\n')
    predicted_sites = list(reader)
    predicted_sites = numpy.array(predicted_sites)
    predicted_sites = predicted_sites.flatten().tolist()
    p_sites_file.close()

    # calculate the overlap
    overlap = calculate_overlap(sites, predicted_sites)

    write_stats(path, run_time, relative_entropy, overlap)

def motif_to_pwm(motif):
    """
    converts the motif from motif.txt to a PWM
    """
    nucleotides = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }

    motif = motif.upper()

    pwm = numpy.zeros((len(motif), 4))

    for i, char in enumerate(motif):
        if char != '*':
            pwm[i][nucleotides[char]] = 1
            temp_dict = nucleotides.copy()
            del temp_dict[char]
            for nucleotide, index in temp_dict.iteritems():
                pwm[i][index] = 0.01
        else:
            for j in xrange(4):
                pwm[i][j] = 0.25

    return pwm

def relate_entropy(m1i, m2i):
    """
    computer the relative entropy between
    the i-th row of m1 and m2
    """
    the_sum = 0
    for j in xrange(4):
        m2ij = float(m2i[j])
        the_sum += m1i[j] * math.log(m1i[j]/m2ij)
        print the_sum
    return the_sum

def compute_relative_entropy(motif, predicted_motif):
    """
    computes the relative entropy between
    motif.txt and predictedmotif.txt
    """
    the_sum = 0
    for i in xrange(len(motif)):
        the_sum += relate_entropy(motif[i], predicted_motif[i])
    return the_sum

def calculate_overlap(sites, predicted_sites):
    """
    calculate the number of overlapping sites
    between sites.txt and predictedsites.txt
    """
    overlap = set(sites).intersection(predicted_sites)
    return len(overlap)

def write_stats(directory, run_time, relative_entropy, overlap):
    f = open(os.path.join(directory, 'stats.txt'), 'w')
    write_str = "{0}\t{1}\t{2}".format(run_time, relative_entropy, overlap)
    f.write(write_str)
    f.close()

if __name__ == '__main__':
    main()
