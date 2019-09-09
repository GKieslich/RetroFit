#!/usr/bin/env python

"""
RetroFit program:
    Fits guest molecules to open metal sites of MOFs using an interaction potential of a
    model system (MIP)."""


__author__ = ["David Bodesheim", "Christian Schneider", "Gregor Kieslich"]
__version__ = "1.0"
__date__ = "06.09.2019"
__email__ = "gregor.kieslich@tum.de"

import matplotlib.pyplot as plt
import numpy as np
import os
import time
import glob
from scipy.interpolate import griddata
import scipy.io
from scipy.interpolate import RegularGridInterpolator
import copy
from argparse import ArgumentParser

# ASE library must be installed prior to usage
from ase.io import read
from ase import Atoms




def credits(__version__, __author__, __date__, __email__):
    print('\n******************************************************************************')
    print('\t\t\t\t*** RetroFit ***\n')
    auths = ", ".join(__author__)
    print('\tAuthors: {}'.format(auths))
    print('\tVersion: {} (Last modified: {})'.format(__version__, __date__))
    print('\tContact: {}\n'.format(__email__))
    print('******************************************************************************\n\n')

    return None

def calc_angle(v1, v2):
    ''''
    Calculates the angle between two vectors v1, v2 in degrees.
    '''
    return np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1) *\
        np.linalg.norm(v2))) * (180./np.pi)

def read_config(workdir, codedir):
    '''
    Reads the config file which is in the directory workdir
    Returns the set of Parameters c_params'''

    os.chdir(workdir)

    with open('RetroFit.config', 'r') as config:
        lines = config.readlines()
        c_params = ['delta_steps',
                'a_steps',
                'rnn_steps',
                'guest_dir',
                'data_file',
                'MMM_angle',
                'RMM']

        for line in lines:
            if line.startswith('#'):
                continue
            for i, cp in enumerate(c_params):
                if line.startswith(cp):
                    c_params[i] = line.split()[1]

    os.chdir(codedir)
    return c_params

def energy_import(path, codedir):
    '''
    Imports energy data of the form [R,theta,delta,E]
    and converts it from Hartree to kcal/mol.
    '''

    # reading in energy file
    data  = np.array([[float(j) for j in i.split() if j != '']\
            for i in open(path).read().split('\n') \
                if i !=''],dtype='float64')
    data[:,3] -= np.min(data[:,3])  # set minimal value as 0
    data[:,3] /= 0.0015936015752548128  # h2kcalmol

    # sorts R, theta and delta and returns unique values only
    R = sorted(np.unique(data[:,0]))
    theta = sorted(np.unique(data[:,1]))
    delta = sorted(np.unique(data[:,2]))

    # creaitng a meshgrid
    Rmesh,Tmesh,Dmesh = np.meshgrid(R,theta,delta, indexing='ij')

    # creating a grid of Energy with values of 0
    E = np.zeros(Rmesh.shape)  # has dimensions of len(R)xlen(theta)xlen(delta)

    # inserting correct Energy values to the Energy grid
    for i, dat in enumerate(data):
        r, t, d, e = dat
        E[R.index(r), theta.index(t), delta.index(d)] = e

    os.chdir(codedir)
    return np.array(R, dtype = float),\
            np.array(theta, dtype = float),\
            np.array(delta, dtype = float),\
            E,\
            data

def import_guest(guest_dir):
    '''
    Reads out the first 4 atoms in the .xyz files of the guest molecule.
    These should be in order N, C, N, C.

    Return:
        atomlist: positions of N, C, N, C
        names: names of xyz - files'''
    os.chdir(guest_dir)
    atomlist = []
    names = []
    for filename in glob.glob('*.xyz'):
        guest = read(filename, format = 'xyz')
        atomlist_tmp = []
        for i in range(4):
            atomlist_tmp.append(guest.get_positions()[i])
        atomlist.append(atomlist_tmp)
        filename = filename.replace('.xyz', '')
        names.append(filename)
    os.chdir('..')
    return atomlist, names

def interpol_reg(R, theta, delta, E):
    '''
    Creates an Interpolation object from a regular Grid
    Input:
        R: list of all R values
        theta: list of all theta values
        delta: list of all delta values
        E: grid of Energy values
    Output:
        interpolation_object
    '''

    # sorting R, theta and delta in case not already sorted
    R_ind = [int(x) for x in range(len(R))]
    R_sort, R_ind = zip(*sorted(zip(R, R_ind)))
    theta_ind = [int(x) for x in range(len(theta))]
    theta_sort, theta_ind = zip(*sorted(zip(theta, theta_ind)))
    delta_ind = [int(x) for x in range(len(delta))]
    delta_sort, delta_ind = zip(*sorted(zip(delta, delta_ind)))

    # rearranges Energy grid according to sorting
    Ecop = copy.deepcopy(E)  # creates a deepcopy of Energy
    for i, rind in enumerate(R_ind):
        for j, tind in enumerate(theta_ind):
            for k, dind in enumerate(delta_ind):
                E[i,j,k] = Ecop[rind,tind,dind]

    interpolation_object = \
            RegularGridInterpolator((R_sort,theta_sort,delta_sort), E)

    return interpolation_object

def calc_from_Mol(atomlist, gamma, delta, RMM):
    '''
    calculates possible RMN values from delta and calculates theta value
    from alpha.
    Input:
        atomlist: list of imported xyz-files
        gamma
        delta: list of delta values
        RMM
    Output:
        RMN: list of list of RMN values (one list for each molecule)
        theta: list of theta values (one value for each molecule)
        delta: list of delta values
        alpha_list: list of alpha values (one value for each molecule)
        R_NN_list: list of R_NN values (one value for each molecule)
    '''

    q = np.pi/180. # conversion factor from deg -> rad

    def calc_RMN(RMM, R_NN, delta, gamma):
        return 0.5 * (RMM - R_NN) / (np.sin((delta - gamma/2.)*q))
    def calc_theta(alpha, gamma):
        return 180 - alpha/2. + gamma/2.

    RMN = []
    theta = []
    alpha_list = []
    R_NN_list = []
    for i in range(len(atomlist)):

        N1 = atomlist[i][0]; C1 = atomlist[i][1]
        N2 = atomlist[i][2]; C2 = atomlist[i][3]

        alpha = calc_angle(np.array(N1) - np.array(C1), \
                    np.array(N2) - np.array(C2))
        R_NN = np.linalg.norm(N1-N2)

        RMN.append(calc_RMN(RMM, R_NN, delta, gamma))
        theta.append(calc_theta(alpha, gamma))

        alpha_list.append(alpha)
        R_NN_list.append(R_NN)

    return RMN, theta, delta, alpha_list, R_NN_list

def output_Molecule_data(int_func, atomlist, names, R_list, theta_list,\
        delta_list, alpha_list, R_NN_list,\
        codedir, data_file, guest_dir, workdir, save_MOL_LIST = False):
    '''
    Screens all possible R, theta and delta values for a molecule
    to obtain optimal values for them.
    Writes out an output file with the resp. information.
    '''
    E_min_list = []; Ropt_list = []; thetaopt_list = []; deltaopt_list = []
    for i in range(len(atomlist)):

        Epath = []
        error_n = 0
        for j in range(len(delta_list)):
            try:
                Epath.append(int_func([R_list[i][j], theta_list[i],\
                    delta_list[j]]))

            except:
                error_n += 1
                Epath.append(np.nan)
                continue

        try:
            Emin = np.nanmin(Epath)
            if isinstance(Emin, np.ndarray):  # in case Emin is np-array
                Emin=Emin[0]
            index = Epath.index(Emin)
            Ropt = R_list[i][index]
            thetaopt = theta_list[i]
            deltaopt = delta_list[index]
            print('\n')
            print(names[i])
            print('-----------------')
            print('alpha : ', alpha_list[i])
            print('R_NN : ', R_NN_list[i])
            print('Emin : ', Emin)
            print('R_opt : ', Ropt)
            print('theta_opt : ', thetaopt)
            print('delta_opt : ', deltaopt)
            print(str(error_n) + ' of ' + str(len(delta_list)) + ' values '
                    'were out of bounds for interpolation!')
            E_min_list.append(Emin)
            Ropt_list.append(Ropt)
            thetaopt_list.append(thetaopt)
            deltaopt_list.append(deltaopt)
        except:
            print('\n')
            print(names[i])
            print('-----------------')
            print('!!!Did not work!!! ' + str(error_n) + ' of ' +
                str(len(delta_list)) + ' out of bounds.')
            Ropt_list.append(-1)
            thetaopt_list.append(-1)
            deltaopt_list.append(-1)
            E_min_list.append(-1)
    if save_MOL_LIST == True:
        os.chdir(workdir)
		# replace -1 values with masked value "--"
        Ropt_list = np.ma.masked_equal(Ropt_list, -1)
        thetaopt_list = np.ma.masked_equal(thetaopt_list, -1)
        deltaopt_list = np.ma.masked_equal(deltaopt_list, -1)
        E_min_list = np.ma.masked_equal(E_min_list, -1)
        
        # sort molecules according to energy penalty
        E_min_list, Ropt_list, thetaopt_list, deltaopt_list, names, alpha_list, R_NN_list  = \
        zip(*sorted(zip(E_min_list, Ropt_list, thetaopt_list, deltaopt_list, names, alpha_list, R_NN_list)))
		


        with open('MOL_LIST.out', 'w') as out:
            out.write('### with codedir (dir where code was executed): {}\n'\
                    .format(codedir))
            out.write('### with data_file: {}\n'.format(data_file))
            out.write('### with guest_dir: {}\n'.format(guest_dir))
            out.write('### with '
                    'len(delta_list) = {}, deltamin = {}, deltamax = {}\n'\
                            .format(len(delta_list),\
                    np.min(delta_list), np.max(delta_list)))
            out.write('### name | alpha | R_NN | Ropt | thetaopt | deltaopt |'
                    'Emin\n')


            for i in range(len(names)):

                out.write(str(names[i]) + '\t\t' + str(alpha_list[i]) + '\t' +\
                        str(R_NN_list[i]) + '\t'  + str(Ropt_list[i]) + '\t' +\
                        str(thetaopt_list[i]) + '\t' + str(deltaopt_list[i]) +\
                        '\t' + str(E_min_list[i]) + '\n')
        os.chdir(codedir)
    return names, alpha_list, R_NN_list, Ropt_list,\
            thetaopt_list, deltaopt_list, E_min_list

def calc_from_Mol_noatominput(alpha, R_NN, gamma, delta, RMM):
    '''
    Calculates theta and RMN with no guest molecule input needed.
    Input:
        alpha
        R_NN
        gamma
        delta
        RMM
    Output:
        RMN
        theta,
        delta'''

    q = np.pi/180.  # conversion factor for deg -> RAD

    def calc_RMN(RMM, R_NN, delta, gamma):
        return 0.5 * (RMM - R_NN) / (np.sin((delta-gamma/2.)*q))
    def calc_theta(alpha, gamma):
        return 180 - alpha/2. + gamma/2.

    RMN = calc_RMN(RMM, R_NN, delta, gamma)
    theta= calc_theta(alpha, gamma)

    return RMN, theta, delta

def output_energysurface(int_func, delta_list, gamma, RMM,\
        codedir, workdir, data_file,\
        a_min, a_max, a_steps, rnn_min, rnn_max, rnn_steps):
    '''
    Writes out a E_SURF.out file which contains Energy-alpha-R_NN pairs.
    '''

    # creating testlists of alpha and R_NN from which to scan
    alpha_testlist = np.linspace(a_min,a_max,a_steps, dtype = 'float')
    R_NN_testlist = np.linspace(rnn_min,rnn_max,rnn_steps, dtype = 'float')



    # looping trough alpha_testlist and R_NN_testlist
    # at each loop, a set of R(list), theta, delta(list) is created
    E_min_plot_list = []; alpha_plot_list = []; R_NN_plot_list = []
    coord_list = []
    for i, aa in enumerate(alpha_testlist):
        for j, RnnRnn in enumerate(R_NN_testlist):
            R_n, t_n, d_n = calc_from_Mol_noatominput(aa, RnnRnn, gamma, \
                delta_list, RMM)

            # looping trough R-delta pairs
            E_list = []
            coords = []
            for k, dd in enumerate(d_n):
                try:
                    E_list = np.append(E_list, int_func([R_n[k], t_n, dd]))
                    coords.append([R_n[k], t_n, dd])

                except:
                    E_list = np.append(E_list, np.inf)
                    coords.append([np.inf, np.inf, np.inf])
            E_min = np.nanmin(E_list)
            E_min_plot_list.append(E_min)
            ind = np.where(E_list == E_min)[0]
            alpha_plot_list.append(aa)
            R_NN_plot_list.append(RnnRnn)
            if E_min==np.inf:
                coord_list.append(np.inf)
            else:
                coord_list.append(coords[int(ind[0])])


    for i, ee in enumerate(E_min_plot_list):
        if ee == np.inf:
            E_min_plot_list[i] = '--'
            coord_list[i] = ['--', '--', '--']

    os.chdir(workdir)
    with open('E_SURF.out', 'w') as out:
        out.write('### with data_file: {}\n'.format(data_file))
        out.write('### with codedir (dir where code was executed): {}\n'\
                .format(codedir))
        out.write('### with '
                'len(delta_list) = {}, deltamin = {}, deltamax = {}\n'\
                        .format(len(delta_list),\
                np.min(delta_list), np.max(delta_list)))
        out.write('### alpha \t| R_NN \t| E_min \t| '
                'R_CuN \t| theta \t| delta \n')
        for i in range(len(alpha_plot_list)):
            out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(alpha_plot_list[i],\
                R_NN_plot_list[i], E_min_plot_list[i],\
                coord_list[i][0], coord_list[i][1], coord_list[i][2]))
    os.chdir(codedir)
    return None

def eval_alpha_R_NN_boundaries(gamma, RMM, R, theta, delta,\
        R_steps, delta_steps):
    '''Calculates boundaries for alpha and R_NN values for the dataset of R,
    theta and gamma values. For R_NN it is calculated by screening all
    possible values of R and delta and taking max and min.
    Input:
        gamma
        RMM
        R: list of imported R-values
        theta: list of imported theta-values
        delta: list of imported delta-values
        R_steps: number of steps for evaluation. The higher, the more precise.
        delta_steps: number of steps for evaluation. The higher, the more
                    precise.
    Output:
        alpha_min
        alpha_max
        R_NN_min
        R_NN_max

    '''
    q = np.pi/180. # conversion from deg -> rad

    # change theta according to symmetry
    theta_sym = np.append(theta, 360-theta)
    # evaluation of alpha boundaries
    alpha_min = 360. + gamma - 2*np.max(theta_sym)
    alpha_max = 360. + gamma - 2*np.min(theta_sym)

    #  evaluation of R_NN boundaries
    def eval_R_NN(R, delta, RMM, gamma):
        return RMM - 2*np.sin((delta - 0.5*gamma)*q) * R
    R_NN_list = []
    R = np.linspace(np.min(R), np.max(R), R_steps)
    delta = np.linspace(np.min(delta), np.max(delta), delta_steps)

    # double for loop could probably improved
    R_lst = []
    d_lst = []
    for r in R:
        for d in delta:
            R_lst.append(r)
            d_lst.append(d)
            R_NN_list.append(eval_R_NN(r, d, RMM, gamma))


    R_NN_min = np.min(np.array(R_NN_list))
    R_NN_max = np.max(np.array(R_NN_list))

    return alpha_min, alpha_max, R_NN_min, R_NN_max


def main():
    t0_tot = time.clock()  # timing
    t_MOL_LIST = None; t_E_SURF = None  # timing
    
    credits(__version__, __author__, __date__, __email__)

    codedir = os.getcwd()  # codedir is current working directory


    # parsing arguments from command line input
    parser = ArgumentParser()
    parser.add_argument("workdir", nargs = '?')
    parser.add_argument("-sMOL", "--save_MOL_LIST", action = 'store_true',\
            default = False)
    parser.add_argument("-sE", "--save_E_SURF", action = 'store_true',\
            default = False)

    args = parser.parse_args()


    # reading in parser arguments
    workdir = args.workdir
    save_MOL_LIST = args.save_MOL_LIST
    calc_E_SURF = args.save_E_SURF


    # reading in .config file
    delta_steps, a_steps, rnn_steps, guest_dir, data_file, MMM_angle, RMM = \
            read_config(workdir, codedir)


    # making RMM to float and calculating gamma
    RMM = np.array(RMM, dtype = float)
    print('RMM = {}'.format(RMM))
    gamma = 180 - 2*(180 - np.array(MMM_angle, dtype = float))
    print('gamma = {}'.format(gamma))


    # importing energy file
    R, theta, delta, E, data = energy_import(data_file, codedir)


    # importing all .xyz file from guest molecules
    atomlist, names = import_guest(guest_dir)

    # creating interpolation object of regular interpolation
    int_func = interpol_reg(R, theta, delta, E)

    # defining a delta-testlist to screen through
    delta_list = np.linspace(np.min(delta), np.max(delta), delta_steps,\
            dtype = float)
            
    if save_MOL_LIST == True:
        t0_MOL_LIST = time.clock()  # timing
    # calculate alpha and all possible RNN (from delta_list) for guest molecules
    R_list, theta_list, delta_list, alpha_list, R_NN_list = \
            calc_from_Mol(atomlist, gamma, delta_list, RMM)


    # finding optimal RNN and alpha and writing out the data of molecules
    names, alpha_list, R_NN_list, Ropt_list, \
          thetaopt_list, deltaopt_list, E_min_list =\
    output_Molecule_data(int_func, atomlist, names, R_list, theta_list,\
    delta_list, alpha_list, R_NN_list,\
    codedir, data_file, guest_dir, workdir, save_MOL_LIST)
    if save_MOL_LIST == True:
        t_MOL_LIST = time.clock() - t0_MOL_LIST  # timing

    # if set true: writing out a energysurface alpha vs RNN for MOF system
    if calc_E_SURF == True:
        t0_E_SURF = time.clock()  # timing
        
        alpha_min, alpha_max, R_NN_min, R_NN_max =\
             eval_alpha_R_NN_boundaries(gamma, RMM, R, theta, delta_list,\
         1000, 1000)

        print('\n\nCalculated boundaries for alpha: ' + str(alpha_min) + \
                     ', ' + str(alpha_max))
        print('Calculated boundaries for R_NN: ' + str(R_NN_min) + \
                ', ' + str(R_NN_max))

        output_energysurface(int_func, delta_list, gamma, RMM,\
                codedir, workdir, data_file, alpha_min, alpha_max, a_steps,\
                R_NN_min, R_NN_max, rnn_steps)
                
        t_E_SURF = time.clock()-t0_E_SURF  # timing
        
    # printing Timings
    print('\n\n------------------------------------------------------------------')
    if t_MOL_LIST!=None:
        print('*Time for calculation and exporting of molecule list: {:.2f}s'.format(t_MOL_LIST))
    if t_E_SURF!=None:
        print('*Time for calculation and exporting of energy surface: {:.2f}s'.format(t_E_SURF))
    t_fin = time.clock()
    print('*Total Time: {:.2f}s\n'.format(t_fin-t0_tot))

    return None


if __name__ == '__main__':

    main()
