import subprocess, os, sys

def main():
	iterations = sys.argv[1]
	for root, dirs, files in os.walk('data'):
	    for directory in dirs:
	    	#get path
	        path = "data/{}".format(directory)
	        print path
	        args = ['nohup', './evaluate.py', path, iterations, '&']
	        print args
	        p = subprocess.Popen(args)

main()