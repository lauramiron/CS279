from rosetta import *
from toolbox import *
from math import exp
from random import random,randint,gauss

#############################################################
######          Overview of code structure             ######
# The Predictor class implements a Monte Carlo Marcov Chain
# optimization procedure for protein structure prediction.

# The classes DihedralPredictor and FragmentPredictor
# subclass* Predictor and redefine its sampleMove method.

# All of the changes you make will be to the Predictor class,
# but we encourage you to look through the other two classes!

# * Subclassing means that they inherit all of the methods
# defined in Predictor (for example predict and __init__) and
# override any methods that are redefined.
#############################################################

init()

class Predictor(object):
  def __init__(self,
               applyPyMOL = False,
               nIters = 1000,
               minFreq = 10,
               name = "",
               full_atom = True,
               ):
    self.forcefield = get_fa_scorefxn() if full_atom else get_cen_scorefxn()
    if applyPyMOL:
      self.pmm = PyMOL_Mover()
      self.pmm.keep_history(True)
    else:
      self.pmm = None

    self.nIters = nIters
    self.dumpFreq = self.nIters/10
    self.name = "{}-{}".format("fa" if full_atom else "centroid",name)

  # Method:  Accept Move
  # --------------------
  # Return True to accept the move and False to reject
  # according to the Metropolis Criterion
  #
  # Assume kT = 1.0
  #
  def acceptMove(self, currE, newE):
    ### BEGIN YOUR CODE HERE ###
    ## around 2 lines of code ##
    if newE < currE: return True
    else: return random() < exp(-(newE-currE))

    ###  END YOUR CODE HERE  ###
    

  # Method:  Standard Monte Carlo Method
  # ------------------------------------
  # Accepts a starting pose, start, and use a Monte Carlo
  # algorithm to predict the most likely pose, optPose.
  # 
  # Useful methods:
  # self.forcefield(pose) -- evaluates "energy" of pose
  # pose.assign(newPose) -- assigns newPose's coordinates on pose
  # self.sampleMove(pose) -- returns new pose after sampling a move
  # self.acceptMove(currE, newE) -- returns True if and only if move should be accepted
  #
  def predict(self, start):
    pose = Pose()
    pose.assign(start)
    optPose = Pose()
    optPose.assign(pose)

    currE = self.forcefield(pose)
    optE = currE

    # main Monte Carlo Loop
    for i in xrange(self.nIters):
      ### BEGIN YOUR CODE HERE ###
      # use sampleMove(), acceptMove()
      # update pose, currE, optPose, and optE as appropriate
      ## around 8 lines of code ##
      newPose = self.sampleMove(pose)
      newE = self.forcefield(newPose)
      if self.acceptMove(currE, newE):
        pose.assign(newPose)
        currE = newE
        if currE < optE:
          optE = currE
          optPose.assign(pose)

      ###  END YOUR CODE HERE  ###
      if self.pmm:
        self.pmm.apply(pose)
      if i % self.dumpFreq == 0:
        print "{} dumping pdb at iteration {}".format(self.name, i)
        optPose.dump_pdb("out/{}-{}-{}.pdb".format(self.name, self.nIters, i))
    optPose.dump_pdb("out/{}-{}-final.pdb".format(self.name, self.nIters))
    return optPose

#########################################################################
##### Below classes subclass Predictor and define sampling strategy #####
#########################################################################

class DihedralPredictor(Predictor):

  def __init__(self,
               applyPyMOL = False,
               nIters = 1000,
               name = "dihedral",
               full_atom = True,
               variance = 25,
               restrict_residues = None,
               ):
    self.variance = variance
    self.restrict_residues = restrict_residues
    super(DihedralPredictor,self).__init__(applyPyMOL = applyPyMOL,
                                           nIters = nIters,
                                           name = name,
                                           full_atom = full_atom,
                                           )

  # Method:  sample dihedral move
  # -----------------------------
  # Samples one of three styles of dihedral moves uniformly
  # at random.  Either increments phi by delta, psi by delta,
  # or phi by delta and psi by -delta.
  # Returns new pose after move.
  # 
  # Useful methods:
  # pose.phi(i) -- returns phi of residue i
  # pose.psi(i) -- returns psi of residue i
  # pose.set_phi(i,newPhi) -- sets phi_i to be newPhi
  # pose.set_psi(i,newPsi) -- sets psi_i to be newPsi
  #
  def sampleMove(self, pose):
    newPose = Pose()
    newPose.assign(pose)
    # select a random residue
    if not self.restrict_residues:
      res = randint(1,len(newPose.sequence()))
    else:
      i = randint(0,len(self.restrict_residues)-1)
      res = self.restrict_residues[i]
    # select the size of the deviation
    delta = gauss(0,self.variance)
    # randomly select the style of dihedral move
    r = random()
    if r < 0.33333:
      # dihedral move in phi
      muPhi = newPose.phi(res)
      newPose.set_phi(res,muPhi + delta)

    elif r < 0.66667:
      # dihedral move in psi 
      muPsi = newPose.psi(res)
      newPose.set_phi(res,muPsi + delta)

    else:
      # shear move
      # update phi <- phi + delta and psi <- psi - delta
      muPhi = newPose.phi(res)
      muPsi = newPose.psi(res)

      newPose.set_phi(res, muPhi + delta)
      newPose.set_psi(res, muPsi - delta)

    return newPose


class FragmentPredictor(Predictor):

  def __init__(self,
               applyPyMOL = False,
               nIters = 100,
               name = "frag",
               full_atom = True,
               fragSize = 9,
               fragFile = "",
               ):
    super(FragmentPredictor,self).__init__(applyPyMOL = applyPyMOL,
                                           nIters = nIters,
                                           name = "%s%d"%(name,fragSize),
                                           full_atom = full_atom,
                                           )
    fragset = ConstantLengthFragSet(fragSize)
    fragset.read_fragment_file(fragFile)
    movemap = MoveMap()
    movemap.set_bb(True)
    self.fragMover = ClassicFragmentMover(fragset,movemap)

  # Method: sample fragment move
  # ----------------------------
  # Samples a move by replacing the coordinates for a randomly
  # selected fragment of the pose with those in the fragment
  # set.
  #
  # Useful Methods:
  # fragMover.apply(pose)
  def sampleMove(self,pose):
    newPose = Pose()
    newPose.assign(pose) 
    self.fragMover.apply(newPose)
    return newPose
