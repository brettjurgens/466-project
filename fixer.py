#!/usr/bin/python

import os, csv, numpy
from evaluate import calculate_overlap

def main():
    for root, dirs, files in os.walk('data'):
        for directory in dirs:
            path = "data/{}".format(directory)

            stats_file = open(os.path.join(path, 'stats.txt'), 'r')
            overlap = int(stats_file.readline().split('\t')[2])
            stats_file.close()

            if overlap < 2:            
                print directory

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

                predicted_sites = [int(elem) + 1 for elem in predicted_sites]

                sites = map(int, sites)

                print sites
                print predicted_sites

                # calculate the overlap
                overlap = calculate_overlap(sites, predicted_sites)

                print overlap

                if overlap > 0:
                    predicted_sites = map(str, predicted_sites)
                    f = open(os.path.join(path, 'predictedsites.txt'), 'w')
                    for site in predicted_sites:
                        f.write(str(site) + '\n')
                    f.close()

                    stats_file = open(os.path.join(path, 'stats.txt'), 'r+')
                    old = stats_file.readline().split('\t')
                    write_str = "{0}\t{1}\t{2}".format(old[0], old[1], overlap)
                    stats_file.seek(0)
                    stats_file.write(write_str)
                    stats_file.close() 

if __name__ == '__main__':
    main()
