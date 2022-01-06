import argparse
from yaml import load, dump
from yaml import CLoader as Loader, CDumper as Dumper
import math
import phonopy
from numpy import pi
import pandas as pd
import argparse

def generate_yaml(supercell_size=[1,1,1],mesh_grid=[1,1,1]):
    print(supercell_size)
    print(mesh_grid)
    #INPUT READ
    data=phonopy.load(supercell_matrix=supercell_size,
                    supercell_filename='SPOSCAR',
                    primitive_matrix='auto',
                    symmetrize_fc=True,
                    force_constants_filename='FORCE_CONSTANTS')

    '''Frequencies and eigenvectors generation'''
    data.run_mesh(mesh=mesh_grid, with_eigenvectors=True,
                    is_mesh_symmetry=True, is_gamma_center=True)
    print('Printing mesh.yaml')
    data.write_yaml_mesh()



parser = argparse.ArgumentParser()
parser.add_argument('--supercell_size',     nargs='+',      help='Supercell size, default size is [1, 1, 1  ]',
                    type=int    ,required=False)
parser.add_argument('--mesh_grid',          nargs='+',      help='Phonopy mesh grid size, default mesh is [1, 1, 1  ]',
                    type=int    ,required=False)
args=parser.parse_args()
generate_yaml(supercell_size=args.supercell_size, mesh_grid=args.mesh_grid)
