import subprocess, os, sys

def main():
	iterations = sys.argv[1]
	for root, dirs, files in os.walk('data'):
	    for directory in dirs:
	        args = ['nohup', './evaluate.py', directory, iterations, '&']
	        print args
	        p = subprocess.Popen(args)

main()