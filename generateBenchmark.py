#!/usr/bin/python

import random, os, datetime

def main():
    # constants... for now...
    motif_length = 8
    num_variable = 1
    seq_length = 500
    seq_count = 10
    directory = generate_directory()

    random_sequences = generate_random_sequences(seq_count, seq_length)
    motif = generate_random_motif(motif_length, num_variable)
    binding_sites = generate_binding_sites(motif, seq_count)
    plant_sites = plant_binding_sites(random_sequences, binding_sites, seq_length, motif_length)
    write_sequences(directory, random_sequences)
    write_sites(directory, plant_sites)
    write_motif(directory, motif_length, motif)
    write_motif_length(directory, motif_length)
    
def generate_directory():
    """
    creates a diretory for the data...
    """
    base_dir = './data'
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    now = datetime.datetime.now()
    date_dir = now.strftime('%Y-%m-%d-%H%M%S')

    formatted_dir = base_dir + '/' + date_dir
    if not os.path.exists(formatted_dir):
        os.makedirs(formatted_dir)

    return formatted_dir

def generate_binding_site(motif):
    """
    generates a binding site...
    """
    nucleotides = { 0: 'A', 1: 'C', 2: 'G', 3: 'T' }

    new_motif = motif

    while '*' in new_motif:
        index = new_motif.index('*')
        motif_arr = list(new_motif)
        motif_arr[index] = nucleotides[random.randint(0, 3)]
        new_motif = ''.join(motif_arr)

    return new_motif

def generate_binding_sites(motif, seq_count):
    """
    generates `seq_count` bonding sites...
    """
    binding_sites = []
    for i in xrange(seq_count):
        binding_sites.append(generate_binding_site(motif))
    return binding_sites

def generate_random_motif(motif_length, num_variable):
    """
    generates a random motif of `motif_length` with a random subset of `num_variable` variable spots
    """
    motif = generate_random_sequence(motif_length)
    variable_spots = []

    for i in xrange(num_variable):
        while True:
            # "mark a random subset of the NM positions in this motif (string) as `variable`..."
            if random.randint(0,1) == 1:
                random_index = random.randint(0, motif_length - 1)
                if random_index not in variable_spots:
                    variable_spots.append(random_index)
                    break
            else:
                break

    motif_arr = list(motif)

    for i in xrange(len(variable_spots)):
        variable_index = variable_spots[i]
        motif_arr[variable_index] = '*'

    motif = ''.join(motif_arr)

    return motif

def generate_random_sequence(size):
    """
    generates a random sequence of size `size`
    """
    return ''.join(random.choice('ACGT') for x in range(size))

def generate_random_sequences(count, size):
    """
    generates `count` random sequences of size `size`
    """
    sequences = []
    for i in xrange(count):
        sequences.append(generate_random_sequence(size))
    return sequences

def plant_binding_sites(sequences, binding_sites, seq_length, motif_length):
    """
    plants the binding sites in each of the sequences and returns the sites
    """
    max_start_location = seq_length - motif_length
    sites = []
    for i in xrange(len(sequences)):
        start_site = random.randint(0, max_start_location)
        sequence = sequences[i]
        sequences[i] = sequence[:start_site] + binding_sites[i] + sequence[start_site + motif_length:]
        sites.append(start_site)
    return sites

def write_sequences(directory, sequences):
    f = open(os.path.join(directory, 'sequences.fa'), 'w')
    for i in xrange(len(sequences)):
        f.write('>sequence%d\n' % i)
        f.write(sequences[i] + '\n')
    f.close()

def write_sites(directory, plant_sites):
    f = open(os.path.join(directory, 'sites.txt'), 'w')
    for plant_site in plant_sites:
        f.write(str(plant_site) + '\n')
    f.close()

def write_motif(directory, motif_length, motif):
    f = open(os.path.join(directory, 'motif.txt'), 'w')
    str = 'MOTIF\t{0}\t{1}'.format(motif_length, motif)
    f.write(str)
    f.close()

def write_motif_length(directory, motif_length):
    f = open(os.path.join(directory, 'motiflength.txt'), 'w')
    f.write(str(motif_length))
    f.close()


if __name__ == '__main__':
    main()