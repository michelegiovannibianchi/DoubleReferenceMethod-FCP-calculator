from ase.io import read,write
from ase.calculators.DoubleReferenceMethod.FCPelectrochem2rm_bader import FCP2rm
from ase.calculators.vasp import Vasp
from ase.io.bader import attach_charges
from ase.calculators.DoubleReferenceMethod.utils import average_potential,add_vacuum,read_fermi_nelect,read_vaccum_level,compute_Bader,read_log_fcp,convert_V_to_label

import math
import os
import time
import numpy as np
import sys 


# Define a wrapper class to encapsulate the Double Reference Workflow
class DoubleReferenceEvaluator:
    def __init__(self, workflow_func, **kwargs):
        self.workflow_func = workflow_func
        self.kwargs = kwargs
        self.results = None

    def evaluate(self, atoms):
        # Call the workflow with the atoms and stored kwargs
        self.results = self.workflow_func(snap=atoms, **self.kwargs)
        return self.results



def DoubleReferenceWorkflow(snap,external_bias_vector,calc_neutral_no_vacuum,calc_neutral_vacuum_no_dipole,calc_neutral_vacuum_dipole,calc_charge,guess_extra_electrons,C_guess,restart,V_SHE=None):
    """Function implementing the Double Reference Workflow
    For teoretical foundations, see the original article DOI: https://doi.org/10.1103/PhysRevB.73.165402
    
    Input: 
        - snap: ase atoms, 
            atomic geometry

        - external_bias_vector: list, 
            values of applied potential for which the Double Reference Method will be applied

        - calc_neutral_no_vacuum: ase calculator, 
            calculator for system without a vacuum region and no extra charge

        - calc_neutral_vacuum_no_dipole: ase calculator, 
            calculator for system with a vacuum region, but not dipole corrections

        - calc_neutral_vacuum_dipole: ase calculator, 
            calculator for system with a vacuum region, and with dipole corrections

        - calc_charge: ase calculator, 
            calculator for system without a vacuum region and extra charge

        - guess_extra_electrons: int, 
            initial guess of extra electrons to add to the system to the first point of the external_bias_vector 
            default:  0, i.e. start from the neutral system

        - C_guess: float, 
            initial guess for the capacitance of the interface
            default: 1/80 e/(V A^2) (same default value of the FCP2rm calculator)

        - restart: bool, 
            if True the calculation will start from the calculation of the system with a vacuum region, and with dipole corrections
            (suppose that the calculations of the neutral system without a vacuum region and system with a vaccum region, but not dipole corrections have been already performed and the charge densities are available)
            default: False

        - V_SHE: float, 
            value of the standard hydrogen electrode potential
            default: 4.44 V
        """
    
    if snap is None:
        raise ValueError("Input atomic geometry is required")
    if external_bias_vector is None:
        raise ValueError("Input external_bias_vector is required")
    if calc_neutral_no_vacuum is None:
        raise ValueError("Input calc_neutral_no_vacuum is required")
    if calc_neutral_vacuum_no_dipole is None:
        raise ValueError("Input calc_neutral_vacuum_no_dipole is required")
    if calc_neutral_vacuum_dipole is None:
        raise ValueError("Input calc_neutral_vacuum_dipole is required")
    if calc_charge is None:
        raise ValueError("Input calc_charge is required")
    if guess_extra_electrons is None:
        guess_extra_electrons=0 # Start from the electron number of the neutral system
    if C_guess is None:
        C_guess=1/80# e/(V A^2), default value of the FCP2rm calculator
    if restart is None:
        restart=False
    if V_SHE is None:   
        V_SHE=4.44# V, default value of the standard hydrogen electrode potential

    
    #################   Start calculation of the "Double Reference Method" workflow   #################

    #snap=read('POSCAR',format='vasp')
    element=snap.get_chemical_symbols()
    n_O=element.count('O')
    case_dir=os.getcwd()

    #### 1) No extra charge + no vacuum
    t0 = time.time()
    snap_neutral=snap
    snap_neutral.calc=calc_neutral_no_vacuum

    if restart==False:
        print('Neutral: Etot=',snap_neutral.get_total_energy())

    #Average potential with MacroDensity
    input_file = 'neutral/LOCPOT'
    output_file = 'neutral/Neutral_no_vac_planar.dat'

    pot=average_potential(input_file,output_file)

    # Potential in the middle of the slab
    phi_prime_0_m=min(pot)
    # Potential in the middle of the water region
    phi_prime_0_w=(pot[0]+pot[-1])/2


    #Read Fermi level and electron number from OUTCAR
    path_OUTCAR='neutral/OUTCAR'
    phi_prime_0_f,nelect0=read_fermi_nelect(path_OUTCAR)


    # Write a RECAP file where store relevant info
    with open('RECAP.dat', 'w') as fp:
        fp.writelines('No extra charge + no vacuum: Potential in the slab [V], Potentential in the middle of water region [V], Fermi level [V]: \n')
        fp.writelines(str(phi_prime_0_m))
        fp.writelines('\n')
        fp.writelines(str(phi_prime_0_w))
        fp.writelines('\n')
        fp.writelines(str(phi_prime_0_f))
        fp.writelines('\n')

    #BADER Charge analysis
    bader_folder="neutral/"
    compute_Bader(case_dir,bader_folder)

    #Resume info in xyz

    #Read energy and forces
    snap=read('OUTCAR',format='vasp-out')

    #Attach the electron number
    snap.info["NELECT0"]= nelect0# electrons without extra charge
    snap.info["NELECT"]= nelect0# electrons with extra charge (in this case no extra charge)

    #Attach Bader charge           
    attach_charges(snap, 'ACF.dat')

    write("snap.xyz",snap,format='extxyz')
    os.chdir(case_dir)

    #Move all relevant files outside computing directory
    os.system('cp neutral/OUTCAR OUTCAR_neutral_no_vaccum')
    os.system('cp neutral/snap.xyz OUTCAR_neutral_no_vaccum.xyz')
    os.system('cp neutral/Neutral_no_vac_planar.dat Neutral_no_vac_planar.dat')
    os.system('cp neutral/ACF.dat ACF_neutral_no_vac.dat')
    os.system('cp neutral/BCF.dat BCF_neutral_no_vac.dat')
    os.system('cp neutral/AVF.dat AVF_neutral_no_vac.dat')

    t1 = time.time()

    #### 2) No extra charge + vacuum            

    snap_vacuum=add_vacuum(snap,n_O)

    ####### 2.1) No extra charge + vacuum + no dipole corrections           
    snap_vacuum.calc=calc_neutral_vacuum_no_dipole
    if restart==False:
        print('Neutral vacuum no dipole: Etot=',snap_vacuum.get_total_energy())

    ####### 2.1) No extra charge + vacuum + dipole corrections           
    snap_vacuum.calc=calc_neutral_vacuum_dipole
    print('Neutral vacuum dipole: Etot=',snap_vacuum.get_total_energy())

    #Average potential with MacroDensity
    input_file = 'neutral_vacuum/LOCPOT'
    output_file = 'neutral_vacuum/Neutral_vac_dipole_planar.dat'
    pot=average_potential(input_file,output_file)

    # Potential in the middle of the slab
    phi_prime_m_r=min(pot)# the minimum of the potential is in the fixed layers of the slab

    #Read vacuum level from OUTCAR
    path_OUTCAR='neutral_vacuum/OUTCAR'
    phi_prime_v_r=read_vaccum_level(path_OUTCAR)   

    # Append to the RECAP file other relevant info
    with open('RECAP.dat', 'a') as fp:
        fp.writelines('No extra charge + no vacuum: Potential in the slab [V], Vaccum level [V]:\n')
        fp.writelines(str(phi_prime_m_r))
        fp.writelines('\n')
        fp.writelines(str(phi_prime_v_r))
        fp.writelines('\n')

    os.system('cp neutral_vacuum/OUTCAR OUTCAR_neutral_vacuum_dipole')
    os.system('cp neutral_vacuum/Neutral_vac_dipole_planar.dat Neutral_vacuum_dipole_planar.dat')

    ####### Compute PZC value via double reference method
    #potential in the neutral cell in the middle of water, without a vacuum, referred to the system with a vacuum
    phi_0_w=phi_prime_0_w-phi_prime_0_m+phi_prime_m_r-phi_prime_v_r
    #Fermi level in the neutral cell without a vacuum, referred to the vacuum level
    phi_0_f=phi_prime_0_f-phi_prime_0_w+phi_0_w
    #PZC value referred to the SHE
    U_0_SHE=-V_SHE-phi_0_f


    #Append to the RECAP file other relevant info
    with open('RECAP.dat', 'a') as fp:
        fp.writelines('Recap Neutral cell: PZC vs SHE [V], phi_0_w [V]:\n')
        fp.writelines(str(U_0_SHE))
        fp.writelines('\n')
        fp.writelines(str(phi_0_w))
        fp.writelines('\n')

    temp=read('OUTCAR_neutral_no_vaccum.xyz',format='extxyz')
    temp.info["U_vs_SHE"]=U_0_SHE
    write('OUTCAR_neutral_no_vaccum.xyz',temp,format='extxyz')

    t2 = time.time()
    print((t2-t1)/60)

    #Append timing info for the neutral calculation to RECAP file
    with open('RECAP.dat', 'a') as fp:
        fp.writelines('Timing [min]\n')
        fp.writelines('Neutral no vacuum:\n')
        fp.writelines(str((t1-t0)/60))
        fp.writelines('\n')
        fp.writelines('Neutral vacuum:\n')
        fp.writelines(str((t2-t1)/60))
        fp.writelines('\n')


    #### 3) Extra charge + no vacuum

    #initialize vectors
    nelect_vector=[]
    U_vector=[]

    #Store the point corresponding to no extra charge
    nelect_vector.append(nelect0)
    U_vector.append(U_0_SHE)    
    
    #Use the wavefunction and the charge density of the neutral stystem as a good starting point for the charged calculations
    os.system('mkdir charge')
    os.system('cp -r neutral/CHG charge/')
    os.system('cp -r neutral/CHGCAR charge/')
    os.system('cp -r neutral/WAVECAR charge/') 

    snap_charge=read('POSCAR',format='vasp')

    for V,index_V in zip(external_bias_vector,range(0,len(external_bias_vector))):

        V_label=convert_V_to_label(V)
        
        if index_V==0: #For the first point start from the PZC informations
            tic=time.time()
            cal_FCP=FCP2rm(innercalc=calc_charge,fcptxt=f'log-fcp_U_{V_label}.txt',U=V,NELECT =nelect0+guess_extra_electrons, NELECT0=nelect0,work_ref=V_SHE,C=C_guess, phi_0_w=phi_0_w,max_FCP_iter=20) 
            snap_charge.calc=cal_FCP
            print(snap_charge.get_total_energy())
            toc = time.time()

        else: #For the other points start from the info of the previous points

            #read the last point to calculate a new guess of the capacitance and the number of electrons
            nelect_last, U_last, C_cal =read_log_fcp(f'log-fcp_U_{convert_V_to_label(external_bias_vector[index_V-1])}.txt')
            nelect_vector.append(nelect_last)
            U_vector.append(U_last)

            #Linear fit to estimate the number of electrons to add
            coef=np.polyfit(np.array(U_vector),np.array(nelect_vector), 1)

            #start the FCP calculation
            tic=time.time()
            cal_FCP=FCP2rm(innercalc=calc_charge,fcptxt=f'log-fcp_U_{V_label}.txt',U=V,NELECT =np.polyval(coef,V), C=C_cal, NELECT0=nelect0,work_ref=V_SHE,phi_0_w=phi_0_w,max_FCP_iter=20) 
            snap_charge.calc=cal_FCP
            print(snap_charge.get_total_energy())
            toc = time.time()


        #Move all relevant files outside computing directory
        os.system(f'cp charge/OUTCAR OUTCAR_{V_label}')
        os.system(f'mv charge/snap.xyz OUTCAR_{V_label}.xyz')
        os.system(f'cp charge/Charge_planar.dat Charge_planar_{V_label}.dat')
        os.system(f'cp charge/ACF.dat ACF_{V_label}.dat')
        os.system(f'cp charge/BCF.dat BCF_{V_label}.dat')
        os.system(f'cp charge/AVF.dat AVF_{V_label}.dat')


        #Append timing info for the neutralcharge calculation to RECAP file
        with open('RECAP.dat', 'a') as fp:                
            fp.writelines(f'Charge {V} V:\n')
            fp.writelines(str((toc-tic)/60))# in minutes
            fp.writelines('\n')

    print("Calculation terminated")
    return


def DoubleReferenceWorkflow_PZC(snap,calc_neutral_no_vacuum):
    """Function similar to DoubleReferenceWorkflow but performing only the calculation of the PZC
       Implemented for reference to get the same formats of the output files of the full DoubleReferenceWorkflow
    
    Input: 
        - snap: ase atoms, 
            atomic geometry

        - calc_neutral_no_vacuum: ase calculator, 
            calculator for system without a vacuum region and no extra charge

        """
    
    if snap is None:
        raise ValueError("Input atomic geometry is required")
    
    if calc_neutral_no_vacuum is None:
        raise ValueError("Input calc_neutral_no_vacuum is required")
    

    
    #################   Start calculation of PZC workflow   #################

    #snap=read('POSCAR',format='vasp')
    element=snap.get_chemical_symbols()
    n_O=element.count('O')
    case_dir=os.getcwd()

    #### 1) No extra charge + no vacuum
    t0 = time.time()
    snap_neutral=snap
    snap_neutral.calc=calc_neutral_no_vacuum

    print('Neutral: Etot=',snap_neutral.get_total_energy())


    #Read Fermi level and electron number from OUTCAR
    path_OUTCAR='neutral/OUTCAR'
    phi_prime_0_f,nelect0=read_fermi_nelect(path_OUTCAR)

    #BADER Charge analysis
    bader_folder="neutral/"
    compute_Bader(case_dir,bader_folder)

    #Resume info in xyz

    #Read energy and forces
    snap=read('OUTCAR',format='vasp-out')

    #Attach the electron number
    snap.info["NELECT0"]= nelect0# electrons without extra charge
    snap.info["NELECT"]= nelect0# electrons with extra charge (in this case no extra charge)

    #Attach Bader charge           
    attach_charges(snap, 'ACF.dat')

    write("snap.xyz",snap,format='extxyz')
    os.chdir(case_dir)

    #Move all relevant files outside computing directory
    os.system('cp neutral/OUTCAR OUTCAR_neutral_no_vaccum')
    os.system('cp neutral/snap.xyz OUTCAR_neutral_no_vaccum.xyz')
    os.system('cp neutral/Neutral_no_vac_planar.dat Neutral_no_vac_planar.dat')
    os.system('cp neutral/ACF.dat ACF_neutral_no_vac.dat')
    os.system('cp neutral/BCF.dat BCF_neutral_no_vac.dat')
    os.system('cp neutral/AVF.dat AVF_neutral_no_vac.dat')

    t1 = time.time()
    #print("Timing [min]:",(t1-t0)/60)
   
    return
