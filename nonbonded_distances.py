import sys
from math import sqrt
import matplotlib.pyplot as plt

def read_pdb(fname):
    """
    Return a list of (Resname, Atom Name, (x, y, z))
       for all 'ATOM' entries in FNAME.
    """
    coords = []
    with open(fname) as fp:
        for line in fp:
            # PDB files have fixed width columns!
            if len(line) < 4 or line[:4] != 'ATOM': continue
            fullname = line[12:16].strip()
            resname = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords += [(resname, fullname, (x, y, z))]
    return coords

def extract_atoms(atoms, entry):
    """
    Return a list of (x, y, z) coordinates
    for the given atoms.
    """
    coords = []
    for resname, atomname in atoms:
        for res, atom, coord in entry:
            if res == resname and atom == atomname:
                coords += [coord]
    return coords

def get_minimum_distances(coords1, coords2):
    """
    For each (x, y, z) tuple in coords1,
    return the distance to the closest
    (x, y, z) coordinate in coords2.
    """
    dists = []
    ################################################
    # TODO: Fill in code so that this function
    #       meets the specification given above.
    # Hints:
    #       This should take 1 to 10 lines of code!
    #       Use the sqrt function which we imported from the math module above.

    ################################################
    for c1 in coords1:
      minDistance = sys.maxint;
      for c2 in coords2:
        d = sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2);
        if d < minDistance:
          minDistance = d;
      dists.append(d);
    return dists

####################################################
# You shouldn't have to change anything below here.

# Load atomic coordinates
pdb = []
for struct in sys.argv[1:]:
    print struct
    pdb += [read_pdb(struct)]

# Define positive, negative, and neutrally charged atoms
positive_atoms = [('LYS', 'NZ'),
                  ('ARG', 'NH1'),
                  ('ARG', 'NH2')]

negative_atoms = [('ASP', 'OD1'),
                  ('ASP', 'OD2'),
                  ('GLU', 'OE1'), 
                  ('GLU', 'OE2')]

neutral_atoms = [('VAL', 'CG1'),
                 ('VAL', 'CG2'),
                 ('LEU', 'CD1'),
                 ('LEU', 'CD2')]

# For each entry in pdb, extract
# ... atoms with each charge.
negative_coords, positive_coords, neutral_coords = [], [], []
for entry in pdb:
    positive_coords += [extract_atoms(positive_atoms, entry)]
    negative_coords += [extract_atoms(negative_atoms, entry)]
    neutral_coords  += [extract_atoms(neutral_atoms, entry)]

# For each entry, find the distance from each positive atom
# ... to the nearest negative and neutral atom
pos_neg, pos_neu = [], []
for pos, neg, neu in zip(positive_coords, negative_coords, neutral_coords):
    pos_neg += get_minimum_distances(pos, neg)
    pos_neu += get_minimum_distances(pos, neu)

# Plot minimum distances
plt.hist([filter(lambda x: x < 10, pos_neg),
          filter(lambda x: x < 10, pos_neu)],
         label = ['positive to negative', 'positive to neutral'],
         bins = 20)
plt.ylabel('Count')
plt.xlabel('Distance (Angstroms)')
plt.legend(loc = 'best')
plt.savefig('nonbonded_dists.png')
