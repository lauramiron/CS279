import sys,os

if len(sys.argv) < 3:
  print "Usage: python pack.py <input>.pdb <output>.pdb"
  sys.exit(1)

if not os.path.isfile(sys.argv[1]):
  print "%s is not a valid path to a pdb file"%sys.argv[1]
  sys.exit(1)

from rosetta import *
from toolbox import *

init()

pose = pose_from_pdb(sys.argv[1])

toFA = SwitchResidueTypeSetMover('fa_standard')
toFA.apply(pose)

packer = standard_packer_task(pose)
packer.restrict_to_repacking()
toPacked = PackRotamersMover(get_fa_scorefxn(),packer)
toPacked.apply(pose)

pose.dump_pdb(sys.argv[2])



