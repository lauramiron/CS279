from putils import *
from logger import *
from rosetta import *
import datetime

sys.stdout = Logger()

def run_dihedral_predictor(templatedPose, restrict_residues=False):
	pose = Pose()
	pose.assign(templatedPose)
	P = DihedralPredictor(applyPyMOL = False, nIters = 1000, full_atom = True, restrict_residues=restrict_residues)
	pose = P.predict(pose)
	toFA = SwitchResidueTypeSetMover('fa_standard')
	toFA.apply(pose)
	return pose

# inputs = ['1IAQ']
inputs = ['1AA9','1IAQ','2rh1']
# inputs = ['2rh1']
print '\n\n'
print datetime.datetime.now()

for pdb_id in inputs:
	rcsbPose = pose_from_rcsb(pdb_id)
	for db in ['fullDB','oldDB']:
		flatPose = unfoldPdb(rcsbPose)
		templatedPose = Pose()
		templatedPose.assign(flatPose)
		templatedPose2 = Pose()
		templatedPose2.assign(flatPose)
		if db == 'fullDB':
			templatedPose, unsolved_residues = template_predict(templatedPose)
			templatedPose2, unsolved_residues = template_predict(templatedPose2,chi=True)
		else:
			date = release_date(pdbID(rcsbPose))
			templatedPose, unsolved_residues = template_predict_old(templatedPose,date)			
			templatedPose2, unsolved_residues = template_predict_old(templatedPose2,date,chi=True)			
		templatedPose.dump_pdb(pdb_id+"-"+db+".pdb")
		rmsd = evaluate_result(templatedPose,rcsbPose)
		print pdb_id+" "+db+" "+str(rmsd)
		templatedPose2.dump_pdb(pdb_id+"-"+db+"-chi.pdb")
		rmsd2 = evaluate_result(templatedPose2,rcsbPose)
		print pdb_id+" "+db+" chi "+str(rmsd2)

		### after initial template prediction, compare refinement methods ###
		
		# dihedral predictor only on unsolved residues 
		if unsolved_residues:
			print "unsolved residues: " + str(unsolved_residues)+", calling dihedral predictor on these residues"
			pose = Pose()
			pose.assign(templatedPose)
			pose = run_dihedral_predictor(templatedPose,unsolved_residues)
			rmsd = evaluate_result(pose,rcsbPose)
			print pdb_id+" "+db+" dihedral restricted "+str(rmsd)
			pose.dump_pdb(pdb_id+"-"+db+"-direstricted.pdb")
		else:
			print "no unsolved residues, continuing..."

		# dihedral predictor on all residues
		pose = Pose()
		pose.assign(templatedPose)
		pose = run_dihedral_predictor(templatedPose)
		rmsd = evaluate_result(pose,rcsbPose)
		print pdb_id+" "+db+" dihedral all "+str(rmsd)
		pose.dump_pdb(pdb_id+"-"+db+"-diall.pdb")

		# small mover
		kT = 1.0
		n_moves = 1000
		movemap = MoveMap()
		movemap.set_bb(True)

		pose = Pose()
		pose.assign(templatedPose)
		small_mover = SmallMover(movemap, kT, n_moves)
		small_mover.apply(pose)
		pose.dump_pdb(pdb_id+"-"+db+"-small.pdb")
		rmsd = evaluate_result(pose,rcsbPose)
		print pdb_id+" "+db+" small "+str(rmsd)

		# shear mover
		pose = Pose()
		pose.assign(templatedPose)
		shear_mover = ShearMover(movemap, kT, n_moves)
		shear_mover.apply(pose)
		pose.dump_pdb(pdb_id+"-"+db+"-shear.pdb")
		rmsd = evaluate_result(pose,rcsbPose)
		print pdb_id+" "+db+" shear "+str(rmsd)
  