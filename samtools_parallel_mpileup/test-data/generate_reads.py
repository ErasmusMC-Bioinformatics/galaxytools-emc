#!/usr/bin/env python


import random
import math


__version_info__ = ('1', '0', '0')
__version__ = '.'.join(__version_info__)


class Region:
	def __init__(self,start,stop,sequence):
		self.start = start
		self.stop = stop
		self.sequence = sequence.strip().replace("\n","").replace(" ","")
		if(len(self.sequence) != self.getSpanningLength()):
			print "ERROR: sequence length: "+str(len(self.sequence))+", while spanning region is: "+str(self.getSpanningLength())
			import sys
			sys.exit()
	
	def getSpanningLength(self):
		return abs(self.stop-self.start+1)

class ReadSynthesizer:
	def __init__(self,chromosome):
		self.regions = []
		self.chromosome = chromosome
	
	def addRegion(self,region):
		self.regions.append(region)
	
	def produceReads(self,readDensity = 1,read_length = 50):
		"""
		Produces uniform reads by walking iteratively over self.regions
		"""
		
		mRNA = self.getTotalmRNA()
		spanning_length = self.getRegionSpanningLength()
		n = spanning_length['total'] - read_length + 1
		
		j = 0
		k = 0
		
		for i in range(n):
			#  "alpha is playing the role of k and beta is playing the role of theta"
			dd = max(0,int(round(random.lognormvariate(math.log(readDensity),0.5))))# Notice this is NOT a binomial distribution!!
			
			for d in range(dd):
				sequence = mRNA[i:i+read_length]
				
				if(random.randint(0,1) == 0):
					strand = 0
				else:
					strand = 16
				flag = strand + 0
				
				print "read_"+str(j)+"."+str(i)+"."+str(d)+"\t"+str(flag)+"\t"+self.chromosome+"\t"+str(self.regions[j].start + k)+"\t60\t"+self.getMappingString(read_length,j,k)+"\t*\t0\t0\t"+str(sequence.upper())+"\t*"
			
			spanning_length['iter'][j] -= 1
			if(k >= self.regions[j].getSpanningLength()-1):
				j += 1
				k = 0
			else:
				k += 1
	
	def getMappingString(self,length,j,offset):
		m = 0
		
		out = ""
		
		for i in range(length):
			k = i + offset
			
			if(k >= self.regions[j].getSpanningLength()):
				j += 1
				
				out += str(m)+"M"
				out += (str(self.regions[j].start - self.regions[j-1].stop-1))+"N"
				m = 1
				
				offset = -k
			else:
				m += 1
		
		out += str(m) + "M"
		
		
		return out
	
	def getRegionSpanningLength(self):
		length = {'total':0,'iter':[]}
		for r in self.regions:
			l = r.getSpanningLength()
			length['iter'].append(l)
			length['total'] += l
		return length
	
	def getTotalmRNA(self):
		mRNA = ""
		for r in self.regions:
			mRNA += r.sequence
		return mRNA



if __name__ == "__main__":
	# Artificial SNP
	rs = ReadSynthesizer('chr1')
	rs.addRegion(Region(  0+1, 59+1,'aaataggtcccaaacgttacgca'+'G'+'tctatgcctgacaaagttgcgaccacttcctctgcc'))#c -> G
	rs.addRegion(Region( 60+1,119+1,'ttgtgtgacacgccggagatagg'+'A'+'catcagcaagtacgttaagtacactgaacgaactgg'))#g -> A
	rs.addRegion(Region(120+1,179+1,'aggtttctacatcgtgcgtgatggc'+'C'+'ctaggagaagtgggtgtatctgcacagcataagt'))#t -> C
	rs.addRegion(Region(180+1,239+1,'tataagacggaagtaaagcgtcttc'+'G'+'ccgttcagcaccccacgctcatagtcaatgctgg'))#a -> G
	#rs.addRegion(Region(240+1,299+1,'ttcagcatagtcaagcgccggtggcctccaaaaagacgcactgagtagcttagctacttt'))
	#rs.addRegion(Region(300+1,359+1,'gctccgcttgcggaagcactaagaggagattgaatttccaaatcccccccgatacctgtg'))
	#rs.addRegion(Region(360+1,419+1,'cggtcgctacgtaagtgcgaagttctgttagatacgctccccttagtatatgggcgttaa'))
	#rs.addRegion(Region(420+1,479+1,'tcggaccgtcggtactcactgcattccaggtctcatatagttcgccctagaagcctggga'))
	rs.addRegion(Region(480+1,539+1,'tgaacgttgaacta'+'GCC'+'ctgatgtaaaccccgcgtgccaattccaggcgtcatgggggca'))#tag -> gcc
	#rs.addRegion(Region(540+1,599+1,'acccctcgcagcctccctcttgctgttggtgcctagtatttcatgatttcgagccgacat'))
	rs.produceReads(2,35)
