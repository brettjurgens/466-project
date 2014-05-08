#!/usr/bin/python

import os, csv, numpy, json

def main():
    counter = 0
    data_array = numpy.zeros((7, 3))
    for root, dirs, files in os.walk('data'):
        for directory in dirs:
            index = counter % 7
            path = "data/{}".format(directory)

            stats_file = open(os.path.join(path, 'stats.txt'), 'r')
            stats = stats_file.readline().split('\t')
            stats_file.close()

            for i, stat in enumerate(stats):
                data_array[index][i] += float(stat)

            counter += 1

    print data_array

    # format it...

    d = {}

    for index, row in enumerate(data_array.tolist()):
        d[index] = row

    formatted_json = json.dumps({'data': d}, indent=4)

    with open('compiled_data.json', 'w') as f:
        f.write(formatted_json)

if __name__ == "__main__":
    main()