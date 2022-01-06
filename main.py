import argparse
from yaml import load, dump
from yaml import CLoader as Loader, CDumper as Dumper
import math
import phonopy
import numpy as np
from numpy import pi
import pandas as pd
import argparse
import warnings

def generate_yaml(supercell_size, mesh_grid):
    #INPUT READ
    data=phonopy.load(supercell_matrix=supercell_size,
                    supercell_filename='SPOSCAR',
                    primitive_matrix='auto',
                    symmetrize_fc=True,
                    force_constants_filename='FORCE_CONSTANTS')

    '''Frequencies and eigenvectors generation'''
    data.run_mesh(mesh=mesh_grid, with_eigenvectors=True, is_mesh_symmetry=True,
                    is_gamma_center=True)
    print('Printing mesh.yaml')
    data.write_yaml_mesh()

def read_raw_data(input_data='mesh.yaml'):
    print('Reading mesh.yaml')
    with open(input_data) as input_data:
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

    return mesh_data, atom_types, mass

def find_qpoints(mesh_data):
    with open('BTE.qpoints', 'r') as input:
        BTE_qpoints=input.readlines()
    if len(mesh_data['phonon'])>len(BTE_qpoints):
        raise ValueError('Phonopy k point mesh is not compatible with ShengBTE mesh')
    if len(mesh_data['phonon'])<len(BTE_qpoints):
        warnings.warn('ShengBTE kpoint mesh is bigger than Phonopy mesh')
    match_qpoints=[]
    no_match_qpoints=[]
    for phonopy_index in range(len(mesh_data['phonon'])):
        phonopy_vector=mesh_data['phonon'][phonopy_index]['q-position']
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

def get_rates(i,j,data, natom,mass):
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

def get_scattering_rate(number_qp,qpoint,vmode,scattering_rates):
    bte_freq=scattering_rates[qpoint+vmode*number_qp].split()[0]
    scattering_rate=scattering_rates[qpoint+vmode*number_qp].split()[1]
    return float(bte_freq)/(2*pi), float(scattering_rate)

def write_results(data,atom_types,mass,match_qp,output_file):
    number_qp=len(match_qp)
    vmodes=len(mass)*3
    natom=len(mass)
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
            res=get_rates(phonopy_index,vmode,data=mesh_data, natom=natom, mass=mass)
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
    with open('{}.yaml'.format(output_file), 'w') as output:
        output.write("\n".join(lines))
        output.close()

def convert_to_dataframe(input_file):
    print('Creating DataFrame')
    with open('{}.yaml'.format(input_file)) as output_data:
        read_data=load(output_data,Loader=Loader)
    nqpoints=len(read_data['phonon'])
    vmodes=len(read_data['phonon'][0]['band'])
    Frequencies=[]
    Scattering_rates=[]
    Pr=[]
    Weight=[]
    element_list=read_data['elements']
    n_elements=len(element_list)
    elements_pr = dict()
    for element in read_data['elements']:
        elements_pr[element] = []
    for qpoint in range(nqpoints):
        for vmode in range(vmodes):
            Frequencies.append(float(read_data['phonon'][qpoint]['band'][vmode]['frequency']))
            Scattering_rates.append(float(read_data['phonon'][qpoint]['band'][vmode]['scattering_rate']))
            #Pr.append(read_data['phonon'][qpoint]['band'][vmode]['pr'])
            Weight.append(float(read_data['phonon'][qpoint]['weight']))
            for element in range(n_elements):
                elements_pr[element_list[element]].append(read_data['phonon'][qpoint]['band'][vmode]['atomic-pr'][element])
    #Scattering_rates=[0 if x=='nan' else x for x in Scattering_rates]

    data=pd.DataFrame({'Freq': Frequencies, 'Scr': Scattering_rates, #'Pr': Pr
                       'Weight': Weight})
    for i in element_list:
        data[i]=elements_pr[i]
    return data, element_list

def clean_data(data):
    data=data.sort_values(by=['Scr'],na_position='last',ascending=False)
    mean=data['Scr'].mean()
    for i in range(len(data)):
        comparation=data.iloc[i]['Scr']-mean
        if comparation>=10000:
            data=data.iloc[1:,:]
        else:
            break
    #El cutoff puede ser una variable
    data['Scr']=data['Scr'].where(data['Scr']>1e-4)
    data.dropna(subset = ["Scr"], inplace=True)
    data.dropna(subset = ["Freq"], inplace=True)
    return data

def create_mesh(data,grid_size):
    gridx,gridy=np.meshgrid(np.linspace(data['Freq'].min(), data['Freq'].max(), grid_size[0]),np.logspace(np.log10(data['Scr'].min()), np.log10(data['Scr'].max()), grid_size[1]))
    return gridx,gridy

def create_polygons(grid_x,grid_y):
    index=0
    verts=[]
    for vecx in gridx:
        for i in range(len(vecx)-1):
            a=[[vecx[i],vecx[i],vecx[i+1],vecx[i+1]],
                                   [gridy[index+1][i],gridy[index][i],
                                    gridy[index][i],gridy[index+1][i]]]
            vert=np.array(a).T
            verts.append(vert)
        index+=1
        if index==len(gridy)-1:
            break
    return verts

def assign_polygons(data,polygons):
    import matplotlib.path as mpath
    print('Assigning polygons')
    polygon_order=[]
    checked=[]
    for point in range(len(data)):
        p=data['Freq'].iloc[point],data['Scr'].iloc[point]
        for polygon in range(len(polygons)):
            path = mpath.Path(polygons[polygon],closed=False)
            if path.contains_point(p,radius=1e-9)==True:
                polygon_order.append(polygon)
                if polygon not in checked:
                    checked.append(polygon)
                if polygon in checked:
                    break
    data['Polygon']=polygon_order
    return data

def convert_to_heatmap(data,column,shape,polygon_number):
    data=data.sort_values(by=['Polygon'])
    ocuppied_polygons=data['Polygon'].unique()
    all_polygons=[i for i in range(polygon_number)]
    heatmap_data=[]
    if column!='Scr':
        for poly in all_polygons:
            scr_weight=0
            if poly in ocuppied_polygons:
                for element in range(len(data.loc[data['Polygon']==poly])):
                    scr_weight+=data.loc[data['Polygon']==poly].iloc[element][column]*data.loc[data['Polygon']==poly].iloc[element]['Weight']
                heatmap_data.append(scr_weight)
            else:
                heatmap_data.append(0)
    if column=='Scr':
        for poly in all_polygons:
            if poly in ocuppied_polygons:
                heatmap_data.append(data['Weight'].loc[data['Polygon']==poly].sum())
            else:
                heatmap_data.append(0)

    return np.reshape(heatmap_data,shape)

def multi_plot(input_data,elements,output,robust,range,cmap,interpolation_method):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rc('text', usetex=True)
    mpl.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    elements.insert(0,'Scr')
    nplots=len(elements)
    extent=[0,data['Freq'].max(), data['Scr'].min(),data['Scr'].max()]
    fig, axs = plt.subplots(nrows=nplots, ncols=1,sharex=True, figsize=[6.4,4.8*nplots], squeeze=True)#,
                            #subplot_kw={'xticks': xticks} )#, 'yticks': [], 'yticklabels':[]                              })
    plt.subplots_adjust(hspace =0.0)
    if robust=='True':
        Scr=convert_to_heatmap(data,column='Scr',shape=[len(gridx)-1,len(gridy[0])-1],polygon_number=(len(gridx)-1)*(len(gridy[0])-1))
        max=Scr.max()
        min=Scr.min()
    if robust=='False':
        min=range[0]
        max=range[1]
    for ax, j in zip(axs.flat, elements):
        print('Creating {} heatmap'.format(j))
        matrix=convert_to_heatmap(data,column=j,shape=[len(gridx)-1,len(gridy[0])-1],polygon_number=(len(gridx)-1)*(len(gridy[0])-1))
        im=ax.imshow(matrix,extent=extent, rasterized=True, interpolation=interpolation_method,
                     cmap=cmap,origin='lower',aspect='auto',vmin=min,vmax=max)
        ax.set_yticks([])
        ax2=ax.twinx()
        ax2.tick_params(axis='y',which='major',labelsize=14,length=6)
        ax2.tick_params(axis='y',which='minor',length=3)
        ax.tick_params(axis='x',labelsize=14,length=6,width=1)
        ax2.tick_params(axis='y',which='both',width=1)
        ax2.plot(data['Freq'],data['Scr'],'rx',alpha=0)
        ax2.set_yscale('log')
        ax2.yaxis.set_ticks_position('left')
        ax.set_ylabel(r'$\Gamma_{tot}$ $\big($ ps$^{-1}\big)$',labelpad=35,fontsize=18)
        if j=='Scr':
            ax2.text((data['Freq'].max()*0.85),data['Scr'].max()*0.5,str('Total'),fontsize=18)
        else:
            ax2.text((data['Freq'].max()*0.9),data['Scr'].max()*0.5,str(j),fontsize=18)
    ax.set_xlabel('Frequency (THz)',fontsize=18)
    left, bottom, width, height = ax.get_position().bounds
    cax = fig.add_axes([left, 0.075, width, height * 0.05])
    #nticks=len(ax.get_xticks())
    #ticks=np.linspace(min, max, nticks, endpoint=True)
    cbar=fig.colorbar(im, orientation='horizontal', cax=cax)#, ticks=ticks)
    cbar.set_label('Weight',fontsize=18)
    plt.savefig('{}.pdf'.format(output),bbox_inches='tight', pad_inches=0.2)

def single_plot(input_data,elements,robust,range,cmap,interpolation_method):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rc('text', usetex=True)
    mpl.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    compount_name=''.join(elements)
    elements.insert(0,'Scr')
    print(compount_name)
    extent=[0,data['Freq'].max(), data['Scr'].min(),data['Scr'].max()]
    if robust=='True':
        Scr=convert_to_heatmap(data,column='Scr',shape=[len(gridx)-1,len(gridy[0])-1],polygon_number=(len(gridx)-1)*(len(gridy[0])-1))
        max=Scr.max()
        min=Scr.min()
    if robust=='False':
        min=range[0]
        max=range[1]
    for i in elements:
        print('Creating {} heatmap'.format(i))
        fig,ax=plt.subplots(figsize=[6.4,4.8])
        matrix=convert_to_heatmap(data,column=i,shape=[len(gridx)-1,len(gridy[0])-1],polygon_number=(len(gridx)-1)*(len(gridy[0])-1))
        im=ax.imshow(matrix,extent=extent, rasterized=True, interpolation=interpolation_method, cmap=cmap,origin='lower',aspect='auto',vmin=min,vmax=max)
        ax.set_yticks([])
        ax2=ax.twinx()
        ax2.tick_params(axis='y',which='major',labelsize=14,length=6)
        ax2.tick_params(axis='y',which='minor',length=3)
        ax.tick_params(axis='x',labelsize=14,length=6,width=1)
        ax2.tick_params(axis='y',which='both',width=1)
        ax2.plot(data['Freq'],data['Scr'],'rx',alpha=0)
        ax2.set_yscale('log')
        ax2.yaxis.set_ticks_position('left')
        ax.set_ylabel(r'$\Gamma_{tot}$ $\big($ ps$^{-1}\big)$',labelpad=35,fontsize=18)
        ax.set_xlabel('Frequency (THz)',fontsize=18)
        #nticks=len(ax.get_xticks())
        #ticks=np.linspace(min, max, nticks, endpoint=True)
        cbar=plt.colorbar(im,ax=ax)
        cbar.set_label('Weight',fontsize=18)
        if i=='Scr':
            i='Total'
        plt.savefig('{}_{}.pdf'.format(compount_name,i),bbox_inches='tight', pad_inches=0.2)


parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mode',         nargs='?', help='Modes: \n \t w: write yaml file with the results \n \t r: read an existing yaml file',
                    choices={"r", "w"},             default='w')
parser.add_argument('-p', '--plot_mode', nargs='?', help='Plotting mode, if m only one .pdf will be generated',
                    choices={'m','s'}, default='m')
parser.add_argument('--supercell_size',     nargs='+',      help='Supercell size, default size is [1, 1, 1  ]',
                    type=int    ,required=False,    default=[1,1,1])
parser.add_argument('--mesh_grid',          nargs='+',      help='Phonopy mesh grid size, default mesh is [1, 1, 1  ]',
                    type=int    ,required=False,    default=[1,1,1])
parser.add_argument('--output_yaml',        nargs='?',        help='Output with the scattering rates, partitipation rates and atomic partitipation rates',
                    type=str    ,required=False,     default='results')
parser.add_argument('--multiplot_name',        nargs='?',        help='Multiplot .pdf file name',
                    type=str    ,required=False,     default='multiplot')
parser.add_argument('--heatmap_xy',         nargs='+',        help='Heatmap x,y  grid, default is [10,10]',
                    type=int    , required=False,    default=[5,5])
parser.add_argument("--robust", nargs='?', type=str,
                        default='True', help='If true heatmap range is fixed to  the extreme values')
parser.add_argument('--heatmap_range', type=int, nargs='+', default=[0,1000], required=False,
                    help='Heatmap colorbar minimum and maximum')
parser.add_argument("--cmap", nargs='?', type=str,
                        default='Blues', help='Heatmap colormap, colormaps from Matplotlib')
parser.add_argument("--interpolation_method", nargs='?', type=str,
                        default='bicubic', help='Interpolation method for imshow, all methods from Matplotlib')


args=parser.parse_args()

if args.mode=='w':
    generate_yaml(supercell_size=args.supercell_size, mesh_grid=args.mesh_grid)
    mesh_data,atom_types,mass=read_raw_data(input_data='mesh.yaml')
    match_qpoints=find_qpoints(mesh_data)
    write_results(data=mesh_data, atom_types=atom_types, mass=mass,match_qp=match_qpoints,output_file=args.output_yaml)

data,element_list=convert_to_dataframe(input_file=args.output_yaml)
data=clean_data(data)
gridx,gridy =create_mesh(data,grid_size=args.heatmap_xy)

polygons=create_polygons(gridx,gridy)
data=assign_polygons(data,polygons)
data.to_csv('culo', sep='\t')
if args.plot_mode=='s':
    print('Single plot mode')
    single_plot(data,elements=element_list,robust=args.robust,range=args.heatmap_range,cmap=args.cmap,
                interpolation_method=args.interpolation_method)
if args.plot_mode=='m':
    print('Multiplot mode')
    multi_plot(data,elements=element_list,output=args.multiplot_name,robust=args.robust,
                range=args.heatmap_range,cmap=args.cmap, interpolation_method=args.interpolation_method)
