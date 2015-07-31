#!/bin/env python

def convert(input_file,output_file):
	fh2 = open(output_file,'w')
	
	with open(input_file) as fh:
		for line in fh:
			lines = line.strip()
			if len(lines) > 0 and lines[0] != '#':
				fh2.write('chr'+line)
			else:
				fh2.write(line)
	
	fh2.close()

