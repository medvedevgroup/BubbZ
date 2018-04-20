import sys
from os import listdir
from os.path import isfile, join

print "##maf version=1"

missing = 0
mypath = './alignment/'
for f in listdir(mypath):
        p = join(mypath, f)
        if isfile(p):
		h = open(p)
		line = [line.strip() for line in h]
		if len(line) > 1:
			print ''
			print '\n'.join(line)
		else:
			missing += 1
		h.close()

print >> sys.stderr, missing

