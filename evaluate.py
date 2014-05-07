#!/usr/bin/python

import time, numpy, os

def main():
    """
    evaluates the motif-finder on the benchmarks
    """
    # run the benchmarks
    for root, dirs, files in os.walk('data'):
        for directory in dirs:
            start_time = time.time()
            # run the gibbs
            # like here or something...
            run_time = time.time() - start_time

            path = "data/{}".format(directory)
            f = open(os.path.join(path, 'motif.txt'), 'r')
            motif = f.readline().split('\t')[2]

            motif_pwm = motif_to_pwm(motif)

            # @TODO open the predicted motif...
            predicted_motif = ''

            relative_entropy = compute_relative_entropy(motif_pwm, predicted_motif)

            # @TODO open the sites and predicted sites...
            overlap = 0
            # overlap = calculate_overlap(sites, predicted_sites)

            write_stats(path, run_time, relative_entropy, overlap)


def motif_to_pwm(motif):
    """
    converts the motif from motif.txt to a PWM
    """
    nucleotides = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }

    motif = motif.upper()

    pwm = numpy.zeros((len(motif), 4))

    for i, char in enumerate(motif):
        pwm[i][nucleotides[char]] = 1

    return pwm

def relate_entropy(m1i, m2i):
    """
    computer the relative entropy between
    the i-th row of m1 and m2
    """
    the_sum = 0
    for j in xrange(4):
        the_sum += m1i[j] * math.log(m1i[j]/m2i[j])
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
    overlap = set(sites).intersection(predictedsites)
    return len(overlap)

def write_stats(directory, run_time, relative_entropy, overlap):
    f = open(os.path.join(directory, 'stats.txt'), 'w')
    write_str = '%f\t%f\t%i'.format(run_time, relative_entropy, overlap)
    f.write(write_str)
    f.close

if __name__ == '__main__':
    main()
