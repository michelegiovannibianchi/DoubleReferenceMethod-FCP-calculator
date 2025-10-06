from ase.geometry import get_distances
from ase.io import read,write

import macrodensity as md
import math
import os
import numpy as np


def average_potential(input_file,output_file):
    """Function computing the average potential with MacroDensity"""
    df, fig = md.plot_planar_average(
        input_file=input_file,
        output_file=output_file,
        lattice_vector=1,# variable required to compute macroscopic average. Here this variable is set to 1 since we are not interested in macroscopic value
        axis='z', # averaging over the direction normal to the surface
        )
    return df["Planar"].values

###################################################################################

def add_vacuum(snap,n_molecule):
    """ Function adding a vacuum region in the middle of the water region:
        Input 
            - snap: ase object, geometry without a vaccum region
            - n_molecule: int, number of water molecules
        
        Return: ase object, geometry with a vaccum region    """

    ############1) Associate the different atoms to the different water molecules############
    
    snap.wrap()# To be safe
    n_atoms=len(snap)
    
    # Use Oxygen to identify the water molecules
    O_pos=np.zeros((n_molecule,4))#[O atom index, x, y ,z] of the oxygen in the water molecules
    molecule_label=np.zeros(n_atoms)# label to associate the different atoms to the molecules

    symbol=snap.get_chemical_symbols()
    pos=snap.get_positions()
    w_O=0
    molecule_counter=0
    
    #loop to identify oxygen atoms
    for j in range(0,n_atoms):
        if symbol[j]=='O':
            O_pos[w_O,0]=j
            O_pos[w_O,1:]=pos[j]
            w_O=w_O+1
            molecule_counter=molecule_counter+1
            molecule_label[j]=molecule_counter

    #Associate hydrogen atoms to the closest oxygen atom to associate H atoms to the different molecules
    for j in range(0,n_atoms):
        if symbol[j]=='H':
            distance=get_distances(pos[j,:],O_pos[:,1:],cell=snap.get_cell(),pbc=True)
            point_index=np.argmin(distance[1])#index of O atom in array "O_pos"
            index_O=O_pos[point_index,0]#index of O atom in array "pos"
            molecule_label[j]=molecule_label[int(index_O)]
    
    #Add to the atoms the label of the corresponding molecules    
    snap.set_tags(molecule_label)


    ############2)Add vacuum############
    
    box=snap.get_cell()
    c_old=box[2,2]# dimension of the box without a vacuum
    snap.center(vacuum=7.5, axis=2)# along z perpendicular to the slab
    box=snap.get_cell()
    c_new=box[2,2]# # dimension of the box with a vacuum
    
    
    ############3)Identify and restore the water molecules on the eadge of the box############
    
    # When a vacuum region is added, the O-H bond within water molecules on the edge of the box, are broken
    # Here duplicate these "broken" molecules to avoid unphysical broken O-H bond at the water-vacuum interface
    
    symbol=snap.get_chemical_symbols()
    pos=snap.get_positions()
    molecule_tag=snap.get_tags()# Labels to indentify to which molecule atom belong
    
    #Loop over molecules to check if they are on the edge
    for w in range(1, molecule_counter+1):#Index from 1, atoms with molecule_tag[w]=0 do not belong to any water molecule
        index_O=1000000000000# i.e., infinity
        index_H=[]
        count_H=0
        
        #Loop over atoms
        for j in range(0,len(snap)):
        
            if molecule_tag[j]==w:#Found an atom that belongs to molecule with label w
                if symbol[j]=='O':
                    index_O=j
                
                elif symbol[j]=='H':
                    index_H.append(j)
                    count_H=count_H+1
            
                
    
        O_up=0#Flag to identify if the Oxygen atom of the water molecule is in the water region above (or below) the slab
        H_up=np.zeros(count_H)
        
    
        if pos[index_O,2]> c_new/2:#O atom is in the water region above the slab
            O_up=1
        #Loop over H atoms inside the same molecule (Usually 2 H atom per molecule but here more general formulation to consider cases of hydronium and hydroxide ions)
        for kk in range(0,count_H):#
            
            if pos[index_H[kk],2]> c_new/2: #H is in the water region above the slab
                H_up[kk]=1
        
     
    
        if (O_up==0 and max(H_up)==1) or (O_up==1 and min(H_up)==0):
        # Condition 1)O in water region below slab AND one H abow the slab at least  OR 2) O in region above AND one H is below
        # This identifies the condition when a O-H bond was broken adding the vacuum
                if O_up==0:
                    pos[index_O,2]=pos[index_O,2]+c_old #Duplicate O atom on the edge
                
                #Loop over H atoms inside the same molecule
                for kk in range(0,count_H):
                    if H_up[kk]==0:
                        pos[index_H[kk],2]=pos[index_H[kk],2]+c_old #Duplicate H atom on the edge
                
                
                            
    snap.set_positions(pos)
    
    return snap

###################################################################################

def read_fermi_nelect(path_OUTCAR):
    """ Function to read Fermi level and number of electrons from OUTCAR"""   
    with open('neutral/OUTCAR', 'r') as fp:
        # read all lines using readline()
            lines = fp.readlines()
            for row in lines:
            # check if string present on a current line
                word = 'Fermi energy:'
            
                # find() method returns -1 if the value is not found,
                # if found it returns index of the first occurrence of the substring
                if row.find(word) != -1:
                    word=row.split()
                    #Fermi level
                    phi_prime_0_f=float(word[2])
                if row.find('NELECT') != -1:
                    words=row.split()
                    #Number of electrons without extra charge
                    nelect0=float(words[2])
    return phi_prime_0_f,nelect0

###################################################################################

def compute_Bader(case_dir,bader_folder):
    """ Function to compute Bader charge analysis"""

    #It requires Bader code by Henkelman group (https://theory.cm.utexas.edu/henkelman/code/bader) and utility "VTST" tools to handle the Charge Density generated by VASP (https://theory.cm.utexas.edu/vtsttools/scripts.html)
            
    # For NON-VASP users: modify the following lines to adjust the Bader charge calculation depending on your DFT code
            
    os.system(f'cp chgsum.pl {bader_folder}/')# move the utility script inside the calculation directory
    os.system(f'cp bader {bader_folder}/')# move the Bader code inside the calculation directory
    os.chdir(case_dir+f'/{bader_folder}/')# get inside the calculation directory
    os.system('./chgsum.pl AECCAR0 AECCAR2')# add core and valence electron charge
    os.system('./bader CHGCAR -ref CHGCAR_sum')# Bader charge calculation
    return

###################################################################################

def read_vaccum_level(path_OUTCAR):
    """ Function to read vacuum level from OUTCAR"""   
    with open(path_OUTCAR, 'r') as fp:
        # read all lines using readline()
        lines = fp.readlines()
        for row in lines:
            # check if string present on a current line
            word = ' vacuum level on the upper side and lower side of the slab'
        # find() method returns -1 if the value is not found,
        # if found it returns index of the first occurrence of the substring
            if row.find(word) != -1:
                word=row.split()
                #word[12] is the vaccum level in the side above the slab
                #word[13]))is the vaccum level in the side below the slab
                phi_prime_v_r=(float(word[12])+float(word[13]))/2# average between two-side vacuum level
            
    return phi_prime_v_r

###################################################################################

def read_log_fcp(log_file):
    """ Function to read the log file of the FCP calculation"""
    with open(log_file, 'r') as fp:
        fp.readline()# label
        stringa=fp.readline()
        substringa=stringa.split()
        nelect_last=float(substringa[1]) #number of electrons
        U_last=float(substringa[5])# applied potential
        C_cal=float(substringa[11])# calculated capacitance
    return nelect_last, U_last, C_cal

###################################################################################

def convert_V_to_label(V):

    """ Function to convert the numeric format of the applied potential to a label.
        Convert the sign of the potentials to a single letter prefix:
        V > 0 -> 'p'
        V < 0 -> 'm'

        example V = -0.5 V => m05
        
        Input:
            -V: float, potential value
        Output:
            -V_label: str, associated labels       """

   
    mantissa = str(abs(V)).split('.')[0]
    digits = str(abs(V)).split('.')[1]
    if digits == '0':
        digits = ''
    if V > 0:
        V_labels=str('p'+mantissa+digits)
    else:
        V_labels=str('m'+mantissa+digits)

    return V_labels