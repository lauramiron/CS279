from rosetta import *
from predictor import *
from dateutil.parser import parse
from bs4 import BeautifulSoup
import sys,os,glob,urllib2,xmltodict

chi_dict = {
	1 : ["ARG","ASN","ASP","CYS","GLN","GLU","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"],
	# 1 : ["ASN","CYS","GLN","GLU","HIS","ILE","LEU","LYS","MET","PHE","PRO","THR","TRP","TYR","VAL"],
	2 : ["ARG","ASN","ASP","GLN","GLU","HIS","ILE","LEU","LYS","MET","PHE","PRO","TRP","THY"],
	# 2 : ["ASN","GLN","GLU","HIS","LEU","MET","PHE","PRO","TRP","THY"],
	3 : ["ARG","GLN","GLU","LYS","MET"],
	# 3 : ["GLN"],
	4 : ["ARG","LYS"],
	# 4 : ["ARG"],
	5 : ["ARG"]
}

acid_codes = {
	"ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "ASX":"B", "CYS":"C", "GLU":"E", "GLN":"Q", "GLX":"Z", "GLY":"G", "HIS":"H",
	"ILE":"I", "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P", "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V"
}

def template_predict(pose,chi=False):
	unsolved_residues = []
	hits = blast_search_sequence(pose.sequence())
	newPose = None
	i = 0
	while newPose == None:
		besthit = hits[i]
		newPose, unsolved_residues = apply_blast_fragment(pose,besthit,chi)
		i+=1
	return newPose, unsolved_residues


def template_predict_old(pose,date,chi=False):
	unsolved_residues = []
	hits = blast_search_sequence(pose.sequence())
	newPose = None
	i=0
	while newPose == None:
		besthit,i = next_best_hit_before_date(hits,date,i)
		newPose, unsolved_residues = apply_blast_fragment(pose,besthit,chi)
	return newPose, unsolved_residues


def evaluate_result(predictedPose,referencePose):
	calpha_superimpose_pose(predictedPose,referencePose)
	rmsd = all_atom_rmsd(predictedPose,referencePose)
	return rmsd


def unfoldPdb(pose):
	sequence = pose.sequence()
	flatpose = pose_from_sequence(sequence)
	return flatpose


def blast_search_sequence_HTML(sequence):
    url = 'http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence='+sequence+'&eCutOff=10.0&matrix=BLOSUM62&outputFormat=HTML'
    request = urllib2.Request(url)
    response = urllib2.urlopen(request)
    HTML = response.read()
    soup = BeautifulSoup(HTML,'html.parser')
    hits = soup.pre.find_all('pre')[1:-1]
    return hits


def blast_search_sequence(sequence):
	url = 'http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence='+sequence+'&eCutOff=10.0&matrix=BLOSUM62&outputFormat=XML'
	request = urllib2.Request(url)
	response = urllib2.urlopen(request)
	XML = response.read()
	XML = XML.decode('unicode_escape')
	results = xmltodict.parse(XML, process_namespaces=True)
	hits = results['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']
	return hits


def release_date(pdb_id):
	try:
		url = 'http://www.rcsb.org/pdb/rest/describePDB?structureId='+pdb_id
		request = urllib2.Request(url)
		response = urllib2.urlopen(request)
		XML = response.read()
		results = xmltodict.parse(XML, process_namespaces=True)
		date = results['PDBdescription']['PDB']['@release_date']
		return parse(date)
	except Exception as e:
		print e
		print url
		sys.exit(1)

def next_best_hit_before_date(hits,date,start=0):
	for i in range(start,len(hits)):
		hit = hits[i]
		if (release_date(pdbID_from_hit(hit)) < date):
			return hit,i


def hits_before_date(hits,date):
	for hit in hits:
		pdb_id = hit['Hit_def'].split(':')[0].encode("utf-8")
		if release_date(pdb_id) > date:
			hits.remove(hit)
	return hits


def apply_exact_match(pose,templatePose,chi=False):
	apply_template(pose,1,pose.sequence(),templatePose,1,templatePose.sequence(),chi)


def set_chis(pose,templatePose,p,t):
	aacid = pose.residue(p).name().split(":")[0]
	templateaacid = templatePose.residue(t).name().split(":")[0]
	for i in range(1,5):
		if aacid in chi_dict[i] and templateaacid == aacid:
			try:
				pose.set_chi(i,p,templatePose.chi(i,t))
			except:
				print "aacid: " + aacid + ", i: " + str(i)
				sys.exit(1)
	return pose


def apply_template(pose,start,sequence,templatePose,templateStart,templateSequence,chi=False):
	unsolved_residues = []

	# apply all phi, psi, optionally chi angles of templatePose to matched portion of our pose
	q = start
	h = templateStart
	for i in range(0,len(templateSequence)):
		if templateSequence[i] == '-':
			unsolved_residues.append(q-1)
			unsolved_residues.append(q)
			unsolved_residues.append(q+1)
			q+=1
		elif sequence[i] == '-':
			h+=1
			unsolved_residues.append(q-1)
			unsolved_residues.append(q)
			unsolved_residues.append(q+1)
		else:
			try:
				if templatePose.sequence()[h-1] != templateSequence[i]:
					print "BLAST alignment still incorrect after adjustment, abandoning fragment"
					return None, unsolved_residues
				pose.set_phi(q,templatePose.phi(h))
				pose.set_psi(q,templatePose.psi(h))
				pose.set_omega(q,templatePose.omega(h))
				if chi:
					pose = set_chis(pose,templatePose,q,h)
			except Exception as e:
				print "i:" + str(i) + " q:" + str(q) + " h:" + str(h) + "template_total:" + str(templatePose.total_residue())
				print e
				sys.exit(1)
			q+=1
			h+=1
	
	# add residues before beginning or after end of matched portion to unsolved_residues
	for i in (range(0,start)+range(len(templateSequence),len(sequence))):
		unsolved_residues.append(i+1)
	unsolved_residues = list(set(unsolved_residues))

	print "applied fragment "+pdbID(templatePose)+" ("+str(templateStart)+","+str(len(templateSequence))+") to query ("+str(start)+","+str(len(templateSequence))+")"
	return pose, unsolved_residues


def pdbID_from_hit(hit):
	return hit['Hit_def'].split(':')[0].encode("utf-8")


def pdbID(pose):
	return pose.pdb_info().name().split('.')[0].split('-')[0]


# make sure multiple chains are not messing up pose.sequence() and pose.residue(i) indexing
# (this should not happen, I thought it might)
def sanity_check_pose(pose):
	failed = False
	if len(pose.sequence()) != pose.total_residue():
		# print "length mismatch:"
		# print "pose.sequence():" + str(len(pose.sequence())) +", total_residue():" + str(pose.total_residue())
		failed = True
	for i in range(0,len(pose.sequence())):
		try:
			if pose.sequence()[i] != acid_codes[pose.residue(i+1).name().split(":")[0].split('_')[0]]:
				# print "res mismatch"
				# print str(i) + " sequence:" + pose.sequence()[i] +", residue(i):" + pose.residue(i+1).name().split(":")[0]
				failed = True
		except Exception as e:
			print e
			print "CODE: " + pose.sequence()[i]
			sys.exit(1)
	if failed:
		print pdbID(pose) + " residue(i) and sequence[i+1] mismatch, check chain ID. exiting..."
		sys.exit(1)


# BLAST is returning incorrect alignment data where the sequences are correct, but the start index
# of the hit does not correspond to the correct pose.residue(i) number
# it is usually 5-7 indices off and at the beginning of the hit sequence
def check_hit_alignment(hitPose,hitSequence,hitStart):
	sanity_check_pose(hitPose)
	realStart = hitPose.sequence().find(hitSequence[0:6]) + 1 # residue indexed from 1, sequence from 0
	if realStart != hitStart:
		print "BLAST incorrect hsp_hit-from value"
	return realStart


def apply_blast_fragment(pose,hit,chi=False):
	print "attempting blast fragment..."
	if isinstance(hit['Hit_hsps']['Hsp'], list): matches = hit['Hit_hsps']['Hsp'] # bad api return format
	else: matches = [hit['Hit_hsps']['Hsp']]
	for match in matches:		
		queryStart = int(match['Hsp_query-from'])
		hitStart = int(match['Hsp_hit-from'])
		querySequence = match['Hsp_qseq'].encode("utf-8")
		hitSequence = match['Hsp_hseq'].encode("utf-8")
		hitPdbId = pdbID_from_hit(hit)
		hitPose = pose_from_rcsb(hitPdbId)

		# correct for errors in BLAST results
		hitStart = check_hit_alignment(hitPose,hitSequence,hitStart)
		stop = len(hitPose.sequence()) - hitStart

		templateSequence = hitSequence[:stop]
		newPose, unsolved_residues = apply_template(pose,queryStart,querySequence,hitPose,hitStart,templateSequence,chi)
	return newPose, unsolved_residues
	






