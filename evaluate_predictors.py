from rosetta import *
from predictor import *
from pypdb import *
import sys,os,glob

def evaluateResult(predictedPdb,reference):
	predictedPose = pose_from_pdb(predictedPdb)
	if ".pdb" in reference: rcsbPose = pose_from_pdb(reference)
	else: rcsbPose = pose_from_rcsb(reference)
	rmsd = CA_rmsd(predictedPose,rcsbPose)
	print "rms for " + predictedPdb + ": " + str(rmsd)
	return rmsd

def unfoldPdb(infile,outfile):
	pose = pose_from_pdb(infile)
	sequence = pose.sequence()
	flatpose = pose_from_sequence(sequence)
	flatpose.dump_pdb(outfile)
	return flatpose

# # adapted from pypdb
# def get_blast_from_sequence(pdb_id, chain_id='A'):
#     raw_results = get_raw_blast(pdb_id, output_form='XML', chain_id=chain_id)

#     out = xmltodict.parse(raw_results, process_namespaces=True)
#     out = to_dict(out)
#     out = out['BlastOutput']
#     return out

def get_blast_hits(pdbID):
	result = get_blast(pdbID)
	just_hits = result['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']
	return just_hits

def apply_template_fragment(pose,hit):
	queryStart = hit['Hit_hsps']['Hsp']['Hsp_query-from']
	queryEnd = hit['Hit_hsps']['Hsp']['Hsp_query-to']
	hitStart = hit['Hit_hsps']['Hsp']['Hsp_hit-from']
	hitEnd = hit['Hit_hsps']['Hsp']['Hsp_hit-to']
	querySequence = hit['Hit_hsps']['Hsp']['Hsp_qseq']
	hitSequence = hit['Hit_hsps']['Hsp']['Hsp_hseq']
	hitIdentity = hit['Hit_def']
	hitPdbId = hitIdentity.split(':')[0]
	hitPdbId = hitPdbId.encode("utf-8")
	print "HITPDBID: " + hitPdbId
	hitPose = pose_from_rcsb(hitPdbId)

	if len(querySequence) != len(hitSequence):
		print "ERROR: I have misunderstood these blast results"
		sys.exit(1)

	print "len(querySequence): " + str(len(querySequence))
	querySequence = querySequence.encode("utf-8")
	print "len(querySequence): " + str(len(querySequence))
	print querySequence
	q = int(queryStart)
	h = int(hitStart)
	print "Q: " + str(q)
	print "H: " + str(h)
	for i in range(0,len(querySequence)):
		print "start residue, i:" + str(i) + " q:" + str(q) + " h:" + str(h)
		if hitSequence[i] == '-':
			q+=1
		elif querySequence[i] == '-':
			h+=1
		else:
			pose.set_phi(q,hitPose.phi(h))
			pose.set_psi(q,hitPose.psi(h))
			q+=1
			h+=1
	pose.dump_pdb("one_fragment_applied.pdb")
	


# def apply_template_fragment(pose, )

# if os.path.isfile(sys.argv[1]):
# 	# run on argv[1]
# elif os.path.isdir(sys.argv[1]):
# 	for f in glob.glob('*.pdb'):
# 		# run on f
# else:
# 	print "%s is not a valid path"%sys.argv[1]
# 	sys.exit(1)

# predictors = set(("-dihedral", "-frag9", "-frag3", "-all"))
# options = set(sys.argv[3:])

# try:
#   predictorName = predictors.intersection(options).pop()
#   options.remove(predictorName)
# except KeyError:
#   print "Make sure you've specified a predictor out of:"
#   print [p for p in predictors]
#   sys.exit(1)

# pymol = "-pymol" in options
# if pymol:
#   options.remove("-pymol")

# fa = "-fa" in options
# if fa:
#   options.remove("-fa")

# nIters = 1000
# for o in options:
#   try:
#     nIters = int(o[1:])
#     break
#   except ValueError:
#     continue

# if predictorName == "-all":
# 	for name in [n for n in predictors if n != "-all"]:
# 		outfile = argv[1]+name+"-"+nIters
# 		runPredictor(name,argv[1],outfile,fa)
# else:
# 	outfile = argv[1],argv[1]+predictorName+"-"+nIters
# 	runPredictor(predictorName,outfile,fa)


# def runPredictor(predictorName, infile, outfile, fa):
# 	pose = pose_from_pdb(infile)
# 	pose.pdb_info().name("struct")

# 	if not fa:
# 	  toCentroid = SwitchResidueTypeSetMover('centroid')
# 	  toCentroid.apply(pose)

# 	if predictorName == "-dihedral":
# 	  P = DihedralPredictor(applyPyMOL = pymol, nIters = nIters, full_atom = fa)
# 	elif predictorName == "-frag9":
# 	  P = FragmentPredictor(applyPyMOL = pymol, nIters = nIters, fragSize = 9, fragFile = "frags/frag9", full_atom = fa)
# 	elif predictorName == "-frag3":
# 	  P = FragmentPredictor(applyPyMOL = pymol, nIters = nIters, fragSize = 3, fragFile = "frags/frag3", full_atom = fa)

# 	pose = P.predict(pose)
# 	toFA = SwitchResidueTypeSetMover('fa_standard')
# 	toFA.apply(pose)
# 	pose.dump_pdb("out/%s"%outfile)





