from yaml import load, dump
from yaml import CLoader as Loader, CDumper as Dumper
import math
import phonopy
from numpy import pi
import pandas as pd
#PHONOPY MESH GENERATOR
def generate_yaml(supercell_size=[1,1,1],mesh_grid=[1,1,1]):
    supercell_size=supercell_size
    mesh_grid=mesh_grid
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

#PARTITIPATION RATES
def get_rates(i,j,mesh_file, natom,mass):
    data=mesh_file
    natom=natom
    #ATOMIC_PARTITIPATION_RATIO VARIABLES
    weight=float(data['phonon'][i]['weight'])
    NUM_1=float()
    DEM_1=float()
    #RAW_ATOMIC_PARTITIPATION_RATIO VARIABLES
    NUM_2=[]
    DEM_2=float()
    RAW_ATOMIC_PARTITIPATION_RATIO=[]
    ATOMIC_PARTITIPATION_RATIO=[]
    for vec in range(natom):
        a=( data['phonon'][i]['band'][j]['eigenvector'][vec][0][0]**2+
            data['phonon'][i]['band'][j]['eigenvector'][vec][0][1]**2)
        b=( data['phonon'][i]['band'][j]['eigenvector'][vec][1][0]**2+
            data['phonon'][i]['band'][j]['eigenvector'][vec][1][1]**2)
        c=( data['phonon'][i]['band'][j]['eigenvector'][vec][2][0]**2+
            data['phonon'][i]['band'][j]['eigenvector'][vec][2][1]**2)
    #PARTITIPATIO_RATIO
        NUM_1+=(a+b+c)/mass[vec]
        DEM_1+=(a+b+c)**2/mass[vec]**2
    #RAW_ATOMIC_PARTITIPATION_RATIO
        NUM_2.append((a+b+c))
        DEM_2+=(a+b+c)
    #PARTITIPATIO_RATIO
    PARTITIPATION_RATIO=(NUM_1)**2/(natom*DEM_1)
    suma=0.0
    for raw_atom in range(natom):
        RAW_ATOMIC_PARTITIPATION_RATIO.append(NUM_2[raw_atom]/DEM_2)
        suma+=NUM_2[raw_atom]/DEM_2

    index=int()
    for clean_atom in range(natom-1):
        index+=1
        if mass[clean_atom]==mass[index]:
            ATOMIC_PARTITIPATION_RATIO.append((RAW_ATOMIC_PARTITIPATION_RATIO[clean_atom]
                                              +RAW_ATOMIC_PARTITIPATION_RATIO[index]))
    freq=data['phonon'][i]['band'][j]['frequency']
    return PARTITIPATION_RATIO, ATOMIC_PARTITIPATION_RATIO

#MATCH BETWEEN PHONOPY AND BTE QPOINTS
def find_qpoints(phonopy_mesh):
    results=phonopy_mesh
    with open('BTE.qpoints', 'r') as input:
      BTE_qpoints=input.readlines()
    match_qpoints=[]
    no_match_qpoints=[]
    for phonopy_index in range(len(results['phonon'])):
        phonopy_vector=results['phonon'][phonopy_index]['q-position']
        count=0
        for bte_index in range(len(BTE_qpoints)):
            count+=1
            bte_vector=BTE_qpoints[bte_index].split()[3:6]
            error=0
            error=abs(phonopy_vector[0]-float(bte_vector[0]))+abs(phonopy_vector[1]-float(bte_vector[1]))+abs(phonopy_vector[2]-float(bte_vector[2]))
            if error <=0.001:
                match_qpoints.append([phonopy_vector,bte_index])
                break
            if count>=len(BTE_qpoints):
                match_qpoints.append([phonopy_vector,'NO_MATCH'])
                no_match_qpoints.append(phonopy_index)
    qpoinds_df=pd.DataFrame(match_qpoints,columns=['phonopy','bte'])    
    if not no_match_qpoints:
        print('All qpoints were found')
    for i in no_match_qpoints:
        print('No match found for phonopy_index %i' % i)
    return qpoinds_df

#GET SCATERING RATES
def get_scattering_rate(number_qp,qpoint,vmode,scattering_rates):
    n_qpoints,qpoint,vmode=number_qp, qpoint, vmode
    bte_freq=scattering_rates[qpoint+vmode*n_qpoints].split()[0]
    scattering_rate=scattering_rates[qpoint+vmode*n_qpoints].split()[1]
    return float(bte_freq)/(2*pi), float(scattering_rate)
#WRITE RESULTS

def write_results(atoms,data,natom,vmodes,mass,match_qp):
    atom_types, natom, vmodes=atoms, natom, vmodes
    match_qp=match_qp
    data=data
    number_qp=len(match_qp)
    with open('BTE.w_final','r') as input:
        scattering_rates=input.readlines()
    print('Writing results')
    lines=[]
    lines.append("elements:")
    for atom in range(len(atom_types)):
        lines.append('- %s' % atom_types[atom])
    lines.append("phonon:")
    for phonopy_index in range(len(match_qp)):
        lines.append("- q-position: [ %12.7f, %12.7f, %12.7f ]" %tuple(data['phonon'][phonopy_index]['q-position']))
        lines.append("  weight: %d" % data['phonon'][phonopy_index]['weight'])
        lines.append("  band:")
        for vmode in range(vmodes):
            res=get_rates(phonopy_index,vmode,mesh_file=data, natom=natom, mass=mass)
            lines.append('  - #%d' %(vmode+1))
            lines.append('    frequency: %15.10f' % data['phonon'][phonopy_index]['band'][vmode]['frequency'])
            if match_qp['bte'][phonopy_index]!='NO_MATCH':
                bte_freq,scattering_rate=get_scattering_rate(number_qp,phonopy_index,vmode,scattering_rates)
                lines.append('    ShengBTE_frequency: %15.10f' % bte_freq)
                lines.append('    scattering_rate: %15.10f' % scattering_rate)
            lines.append('    pr: %15.10f' % res[0])
            lines.append('    atomic-pr:')
            for k in range(len(atom_types)):
                lines.append('    - # %s:' % str(atom_types[k]))
                lines.append('        %15.10f'% res[1][k])
    with open('res.yaml', 'w') as output:
        output.write("\n".join(lines))
        output.close()

'''
MAIN PROGRAM
'''
def main():
    #mesh.yaml GENERATOR
    generate_yaml(supercell_size=[5,5,3],mesh_grid=[12,12,12])

    #DATA LOADING
    print('Reading mesh.yaml')
    with open('mesh.yaml') as input_data:
        mesh_data=load(input_data,Loader=Loader)

    natom=mesh_data['natom']
    vmodes=natom*3
    mass=[]
    for atom in range(natom):
        mass.append(mesh_data['points'][atom]['mass'])
    atom_types=[]
    index=0
    for n_type in range(natom-1):
        index+=1
        if mass[n_type]==mass[index]:
            atom_types.append(mesh_data['points'][index]['symbol'])

    phonopy_v_bte_index=find_qpoints(phonopy_mesh=mesh_data)
    write_results(atoms=atom_types,data=mesh_data, natom=natom,vmodes=vmodes, mass=mass, match_qp=phonopy_v_bte_index)


main()