#!/usr/bin/env python

# first, fetch all known dbkeys:
# http://genome-test.cse.ucsc.edu/cgi-bin/das/dsn

import urllib2
import xml.etree.ElementTree as ET
import re

def fetch_dbkeys():
	print "Fetching available dbkeys"
	response = urllib2.urlopen('http://genome-test.cse.ucsc.edu/cgi-bin/das/dsn')
	xml = response.read()
	
	root = ET.fromstring(xml)
	
	dbkeys = {}
	for dsn in root.iter('DSN'):
		dbkey = dsn.find('SOURCE').attrib['id']
		full_name = dsn.find('DESCRIPTION').text
		
		dbkeys[dbkey.lower()] = (dbkey,full_name)
	
	return dbkeys

def fetch_liftovers(dbkey_uid, dbkey):
	url_dbkey = "http://hgdownload.cse.ucsc.edu/goldenPath/"+dbkey[0]+"/"
	url_liftover = "http://hgdownload.cse.ucsc.edu/goldenPath/"+dbkey[0]+"/liftOver/md5sum.txt"
	print "Fetching list of liftOver files for dbkey: "+dbkey[0]
	try:
		response = urllib2.urlopen(url_dbkey)
		try:
			response = urllib2.urlopen(url_liftover)
			liftOvers = response.read().strip().split("\n")
			print "   Listed liftOver files for dbkey: "+dbkey[0]
		except:
			print " - Couldn't find liftOver files for dbkey: "+dbkey[0]
			liftOvers = None
	except:
		print " - Couldn't find info at UCSC for dbkey: "+dbkey[0]
		liftOvers = None
	
	return liftOvers

def capitalize_dbkey(dbkey):
	dbkey = list(dbkey)
	for i in range(len(dbkey)):
		if(i % 3 == 0):
			dbkey[i] = dbkey[i].upper()
	
	return ''.join(dbkey)

def __main__():
	dbkeys = fetch_dbkeys()
	
	fh = open("liftOver.loc.available","w")
	fh.write("#<FromSpecies>	<ToSpecies>	<PathToChainFile>	<optional: FromSpeciesDescription>	<optional: ToSpeciesDescription>	<optional: ChainFile md5sum>\n\n")
	
	for dbkey_uid, dbkey in dbkeys.items():
		liftOvers = fetch_liftovers(dbkey_uid, dbkey)
		
		# If the dbkey has corresponding liftOver files, parse them and write them to file
		if(liftOvers):
			for line in liftOvers:
				md5sum, liftOverFile = re.split("[\s]+",line)
				dbkey_to = liftOverFile.split("To",1)[1].lower().replace('.over.chain.gz','')
				if(dbkeys.has_key(dbkey_to)):
					dbkey_to = dbkeys[dbkey_to]
				else:
					dbkey_to = capitalize_dbkey(dbkey_to)
					print "***Warning: unknown dbkey used in liftOver: "+dbkey_to
					dbkey_to = (dbkey_to,"Unknown species")
				liftOverFile = "http://hgdownload.cse.ucsc.edu/goldenPath/"+dbkey[0]+"/liftOver/"+liftOverFile
				fh.write(dbkey[0]+"\t"+dbkey_to[0]+"\t"+liftOverFile+"\t"+dbkey[1]+"\t"+dbkey_to[1]+"\t"+md5sum+"\n")
	
	fh.close()

if __name__ == "__main__":
	__main__()
