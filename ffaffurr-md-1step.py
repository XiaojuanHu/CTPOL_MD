from openmm.app import *
from openmm.app.internal import *
from openmm import *
from openmm.unit import *
from sys import stdout
import os, sys, re
import traceback
from xml.etree import ElementTree as ET
import time
import argparse
#from IPython.core import ultratb
#sys.excepthook = ultratb.FormattedTB(
    #mode='Verbose',
#    color_scheme='Linux', 
#    call_pdb=False
#)
time_start=time.time()

# arguments
parser = argparse.ArgumentParser(
    description="Runs simulation with CT+POL if applicable",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("parameter", help="parameter file, can be OPLS-AA.xml or ffaffurr-oplsaa.xml")
parser.add_argument("-s", "--structure", help="pdb file", default='input.pdb')
args=parser.parse_args()

# parameter file we choose
force_file = args.parameter
# initio structure
try:
    pdb_path=args.structure
except IndexError:
    pdb_path='input.pdb'
    
# basic parameters of the simulation
n_steps = 20000000
time_step = 0.002*picoseconds
trajectory_output_frequency = 500
potential_output_frequency = 10

################amber vdW > oplsaa vdW#########################################################################
def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            #print (p1,p2,sig,eps)
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)  #eps = fudge_factor * esp14
    #for i in range(nonbonded_force.getNumExceptions()):
    #    (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
    #    print(p1, p2, q, sig, eps)
    return system

################### read paramters from CustomForce.xml ###########################################################
def get_Coulomb_factor(force_file):
        tree = ET.ElementTree(file = force_file)
        
        if tree.getroot().find('NonbondedForce') is not None:
                f14 = float(tree.getroot().find('NonbondedForce').attrib['coulomb14scale'])

        return(f14)


def get_custombondforce_para(file='CustomForce.xml'):
        
        type__ChargeTransfer_parameters = {}
        type__zero_transfer_distance = {}
        type__averageAlpha0eff = {}
        type__averageR0eff = {}
        
        tree = ET.ElementTree(file=file)
                        
        if tree.getroot().find('CustomChargeTransfer') is not None:
                CN_factor = str(tree.getroot().find('CustomChargeTransfer').attrib['CN_factor'])
                f15 = float(tree.getroot().find('CustomChargeTransfer').attrib['f15'])
                for atom in tree.getroot().find('CustomChargeTransfer').findall('Atom'):
                        type__ChargeTransfer_parameters[atom.attrib['type']] = [float(atom.attrib['a']), float(atom.attrib['b'])]
                        type__zero_transfer_distance[atom.attrib['type']] = float(atom.attrib['r'])
                        
        if tree.getroot().find('CustomPoleForce') is not None:
                for atom in tree.getroot().find('CustomPoleForce').findall('Polarize'):
                        type__averageAlpha0eff[ atom.attrib['type'] ] = float( atom.attrib['polarizability'] )
                        type__averageR0eff[ atom.attrib['type'] ] = float( atom.attrib['R0'] )
                        
        return(CN_factor, \
                        f15, \
                        type__ChargeTransfer_parameters, \
                        type__zero_transfer_distance, \
                        type__averageAlpha0eff, \
                        type__averageR0eff)

def get_CN(n_atom__xyz, cation_pairs, type__zero_transfer_distance):
        
        CN = 0
        for pair in cation_pairs:
                #print("DEBUG: ",pair)
                distance = get_distance(pair[0], pair[1], n_atom__xyz)
                if n_atom__type[pair[0]] in type__zero_transfer_distance.keys():
                        if distance <= type__zero_transfer_distance[n_atom__type[pair[0]]]:
                                CN += 1
                elif n_atom__type[pair[1]] in type__zero_transfer_distance.keys():
                        if distance <= type__zero_transfer_distance[n_atom__type[pair[1]]]:
                                CN += 1
        #sys.exit()
        return(CN)

def get_atom__type():
        # build system with OpenMM
        
        pdb = PDBFile(pdb_path)
        modeller = Modeller(pdb.topology, pdb.positions)
        forcefield = ForceField('OPLS-AA.xml')
        
        modeller.addSolvent(forcefield, padding=1.0*nanometers)
        
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None, removeCMMotion=False)
        
        # get some atom information
         ## data.atoms, data.atomType, data.atomType[atom]
        try:
            data = forcefield._SystemData(modeller.topology)
        except TypeError:
            data = forcefield._SystemData()
        
         ## Make a list of all atoms ('name', 'element', 'index', 'residue', 'id')
        data.atoms = list(modeller.topology.atoms())
        
         ## Make a list of all bonds
        bondedToAtom = forcefield._buildBondedToAtomList(modeller.topology)
        
         ## data.atomType, residues__atoms
        for chain in modeller.topology.chains():
                
                for res in chain.residues():
                              
                        [template, matches] = forcefield._getResidueTemplateMatches(res, bondedToAtom, ignoreExternalBonds=False)
                        
                        if matches is None:
                                raise Exception('User-supplied template does not match the residue %d (%s)' % (res.index+1, res.name)) 
                        else:
                                data.recordMatchedAtomParameters(res, template, matches)                        
        
         ## n_atom__type                
        n_atom__type = {}
        for atom in data.atoms:
                n_atom__type[atom.__dict__['index']] = data.atomType[atom]
                
        return(n_atom__type)
        
global n_atom__type
n_atom__type = get_atom__type()
#print(n_atom__type)

def forcegroupify(system):
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force.setForceGroup(i)
        forcegroups[force] = i
    return forcegroups

def getEnergyDecomposition(context, forcegroups):
    energies = {}
    for f, i in forcegroups.items():
        energies[f.__class__.__name__] = context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()
    return energies

##################functions for each conformer###############################################################################################

# get n_atom__xyz
def get_pdb__info(file_path):
    pdb = PDBFile(file_path)
    modeller = Modeller(pdb.topology, pdb.positions)
    forcefield = ForceField('OPLS-AA.xml')
        
    modeller.addSolvent(forcefield, padding=1.0*nanometers)
    
    
    coords = modeller.getPositions()
    arr = coords.value_in_unit(angstrom)
    n_atom__xyz = dict()
    for n in range(len(arr)):
        n_atom__xyz[n] = arr[n]

    return n_atom__xyz

# get distance(nm) between two atoms
def get_distance(atom1, atom2, n_atom__xyz):

        distance = math.sqrt(   ( n_atom__xyz[atom1][0] - n_atom__xyz[atom2][0] ) ** 2. \
                                                + ( n_atom__xyz[atom1][1] - n_atom__xyz[atom2][1] ) ** 2. \
                                                + ( n_atom__xyz[atom1][2] - n_atom__xyz[atom2][2] ) ** 2. )
        distance = distance * 0.1                     # angstrom to nm
        return(distance)

# get E0
def get_E0(n_atom__xyz, n_atom__charge, n_atom_nohyd, ion_index):
        
        #type__averageR0eff = get_type__R0eff()
        
        n_atom__E0 = {}
        n_atom__E0[ion_index] = 0
        
        for i in n_atom_nohyd:
                if i != ion_index:
                        try:
                            #print("DEBUG: ion_index")
                            
                            ind1 = n_atom__xyz[ ion_index ]
                            #print("DEBUG: i ")
                            ind2 = n_atom__xyz[ i ]

                            vec_r_iMe = numpy.array(  n_atom__xyz[ ion_index ]  ) * 0.1 \
                                        - numpy.array( n_atom__xyz[ i ] ) * 0.1    # angstrom > nm
                            vec_r_Mei = numpy.array( n_atom__xyz[ i ] ) * 0.1 - numpy.array( n_atom__xyz[ ion_index ] ) * 0.1
                                    
                            r_iMe = get_distance(i, ion_index, n_atom__xyz)
                            
                            # add r_cutoff
                            r_vdw_radius = type__averageR0eff[ n_atom__type[i] ] * 0.1 + type__averageR0eff[ n_atom__type[ion_index] ] * 0.1   # angstrom > nm
                            r_cutoff = 0.92 * r_vdw_radius
                            
                            if r_iMe <= r_cutoff:
                                    r_iMe = r_cutoff
                            
                            factor_i = n_atom__charge[ ion_index ] / ( r_iMe ** 3 )
                            n_atom__E0[i] = factor_i * vec_r_iMe
                            
                            factor_Me = n_atom__charge[ i ] / ( r_iMe ** 3 )
                            n_atom__E0[ion_index] += factor_Me * vec_r_Mei
                        except KeyError as KE:
                            print("DEBUG: KeyError when i, ion_index = ",i,ion_index)
                            traceback.print_exc(file=sys.stdout)
                            #print(n_atom__xyz) 
                            sys.exit("DEBUG EXIT")
        return(n_atom__E0)
        
# get Tij
def get_Tij(n_atom__xyz, i, j):
        
        #type__averageR0eff = get_type__R0eff()
        
        vec_r_ij = numpy.array( n_atom__xyz[ j ] ) * 0.1 - numpy.array( n_atom__xyz[ i ] ) * 0.1         # angstrom > nm
                
        r_ij = get_distance(i, j, n_atom__xyz)
        
        # add r_cutoff
        r_vdw_radius = type__averageR0eff[ n_atom__type[i] ] * 0.1 + type__averageR0eff[ n_atom__type[j] ] * 0.1   # angstrom > nm
        r_cutoff = 0.92 * r_vdw_radius
        
        if r_ij <= r_cutoff:
                r_ij = r_cutoff
                
        I = numpy.eye(3)
        
        rij_rij = numpy.dot(vec_r_ij.reshape(-1,1), vec_r_ij.reshape(1,-1))
        Tij = ( 3  / r_ij ** 2 * rij_rij  - I ) / r_ij ** 3
        return(Tij)

# get induced dipole
def get_induced_dipole(n_atom__xyz, n_atom__alpha0eff, n_atom__charge, n_atom_nohyd, ion_index):
        
        # get E0
        n_atom_nohyd__E0 = get_E0(n_atom__xyz, n_atom__charge, n_atom_nohyd, ion_index)
        
        # initial guess of the induced dipole of cation
        dipole_Me_initial = numpy.array( [0., 0., 0.] ) 
        
        dipole_Me_1 = dipole_Me_initial
        
        # induced dipole
        n_atom_nohyd__dipole ={}
        dipole_Me_2 = n_atom__alpha0eff[ion_index] *  n_atom_nohyd__E0[ion_index]
        for i in n_atom_nohyd:
                if i != ion_index:
                        TiMe = get_Tij(n_atom__xyz, i, ion_index)
                        TMei = get_Tij(n_atom__xyz, ion_index, i)
                        
                        n_atom_nohyd__dipole[i] = n_atom__alpha0eff[i] *  n_atom_nohyd__E0[i] + n_atom__alpha0eff[i] * numpy.dot(dipole_Me_1 , TiMe)
                        dipole_Me_2 += n_atom__alpha0eff[ion_index] * numpy.dot(n_atom_nohyd__dipole[i] , TMei)
                        
        N = 0
        while numpy.linalg.norm( dipole_Me_2 - dipole_Me_1 ) > 2.0819434e-8:           # 1D = 0.020819434 e*nm   10**(-6) D > e*nm 
                
                dipole_Me_1 = dipole_Me_2
                
                n_atom_nohyd__dipole ={}
                dipole_Me_2 = n_atom__alpha0eff[ion_index] *  n_atom_nohyd__E0[ion_index]
                for i in n_atom_nohyd:
                        if i != ion_index:
                                TiMe = get_Tij(n_atom__xyz, i, ion_index)
                                TMei = get_Tij(n_atom__xyz, ion_index, i)
                                
                                n_atom_nohyd__dipole[i] = n_atom__alpha0eff[i] *  n_atom_nohyd__E0[i] + n_atom__alpha0eff[i] * numpy.dot(dipole_Me_1 , TiMe)
                                dipole_Me_2 += n_atom__alpha0eff[ion_index] * numpy.dot(n_atom_nohyd__dipole[i] , TMei)
                N += 1
        
        n_atom_nohyd__dipole[ion_index] = dipole_Me_2
        
        return(n_atom_nohyd__dipole, n_atom_nohyd__E0)
        
def write_CT(sp, total_transq):
    with open('charge_tranfer', 'w') as out:
        out.write('step: ', sp, 'charge transfer: ',total_transq)
        
        
##########################################################################################################################
global f14
global polar_index
global ion_index
#pdb = PDBFile(pdb_path)

print("Initializing...")
pdb = PDBFile(pdb_path)

print("Loading topology and parameters")
modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField(force_file)

#modeller.addSolvent(forcefield, padding=1.0*nanometers)

#system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=None)
print("Creating System")
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,nonbondedCutoff=1.2*nanometer, constraints=HAngles)

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

print("Applying parameters")
if force_file == 'OPLS-AA.xml':
        system = OPLS_LJ(system)
        
# if we use newFF parameters
elif 'ffaffurr-oplsaa' in force_file:
        print("Applying modified parameters and CT a+ POL")
        CN_factor, \
         f15, \
         type__ChargeTransfer_parameters, \
         type__zero_transfer_distance, \
         type__averageAlpha0eff, \
         type__averageR0eff     = get_custombondforce_para(file='CustomForce.xml')
        
        # get n_atom__xyz
        n_atom__xyz = get_pdb__info(pdb_path)
        
        #print("DEBUG", n_atom__xyz) 
        # ion_index, polar_index
        polar_index = []
        for (atom_index, atom) in enumerate(modeller.topology.atoms()):
            if atom_index < 424:
                if atom.element.symbol == 'Zn': # or atom.element.symbol == 'Na':
                        ion_index = atom_index
                elif atom.element.symbol == 'S' or atom.element.symbol == 'O' or atom.element.symbol == 'N':
                        polar_index.append(atom_index)
        
        if type__ChargeTransfer_parameters:

                f14 = get_Coulomb_factor(force_file)
                
                # get CN
                cation_pairs = []               
                for i in polar_index:
                        if i < ion_index:
                                cation_pairs.append((i, ion_index))
                        else:
                                cation_pairs.append((ion_index, i))
                               
                CN = get_CN(n_atom__xyz, cation_pairs, type__zero_transfer_distance)
                
                # load nonbondedforce
                forces = { force.__class__.__name__ : force for force in system.getForces() }
                nbforce = forces['NonbondedForce']
                
                # charge on L
                total_transq = 0
                index__charge = {}
                index__charge_ini = {}
                for index in range(nbforce.getNumParticles()):
                        
                        # charge
                        [charge, sigma, epsilon] = nbforce.getParticleParameters(index)
                        
                        index__charge_ini[index] = charge
                        
                        #if index < 424:
                        if n_atom__type[index] in type__ChargeTransfer_parameters.keys():
                                r = get_distance(index, ion_index, n_atom__xyz)
                                
                                # only consider r <= zero_transfer_distance 
                                if r < type__zero_transfer_distance[n_atom__type[index]] :
                                        
                                        if (CN == 1) or (CN_factor == 'No_Charge_transfer'):
                                                transq = (type__ChargeTransfer_parameters[n_atom__type[index]][0] * r + type__ChargeTransfer_parameters[n_atom__type[index]][1]) * (math.sqrt(f15))
                                        else:
                                                transq = ( type__ChargeTransfer_parameters[n_atom__type[index]][0] * r + type__ChargeTransfer_parameters[n_atom__type[index]][1] ) * (math.sqrt(f15)) / math.pow( CN, 1.0/float(CN_factor) ) 
                                        
                                        total_transq += transq
                                        charge_new = ( charge/charge.unit + transq)* charge.unit
                                        
                                        nbforce.setParticleParameters(index, charge_new, sigma, epsilon)
                                else:
                                        charge_new = charge
                        else:
                                charge_new = charge
                        index__charge[index]=charge_new
                
                [charge, sigma, epsilon] = nbforce.getParticleParameters(ion_index)
                charge_ion = ( charge/charge.unit - total_transq ) * charge.unit
                nbforce.setParticleParameters(ion_index, charge_ion, sigma, epsilon )
                index__charge[ion_index]=charge_ion
                        
                for i in range(nbforce.getNumExceptions()):
                        (p1, p2, q, sig, eps) = nbforce.getExceptionParameters(i)
                        if q._value != 0.0:
                                q_new = index__charge[p1] * index__charge[p2] * f14
                                nbforce.setExceptionParameters(i, p1, p2, q_new, sig, eps)
        
        # implement polarization energy
        if type__averageAlpha0eff:
                
                # implement polarization energy, ev > kJ/mol
                pol_force = openmm.CustomExternalForce(' -0.5 * 96.485309 * ( dipole0 * E0 + dipole1 * E1 + dipole2 * E2) ')      
                
                pol_force.addPerParticleParameter("dipole0")
                pol_force.addPerParticleParameter("dipole1")
                pol_force.addPerParticleParameter("dipole2")
                pol_force.addPerParticleParameter("E0")
                pol_force.addPerParticleParameter("E1")
                pol_force.addPerParticleParameter("E2")
                
                # get n_atom_nohyd list, ion_index
                n_atom_nohyd = []
                for (atom_index, atom) in enumerate(modeller.topology.atoms()):
                    if atom_index < 424:
                        if atom.element.symbol != 'H':          # TIP3
                                n_atom_nohyd.append(atom_index)
                
                # n_atom_nohyd__alpha0eff
                n_atom_nohyd__alpha0eff = {}
                for index in n_atom_nohyd:
                        ind = n_atom__type[index]
                        #print(index, ind, end=" | ")
                        try:
                            n_atom_nohyd__alpha0eff[index] = type__averageAlpha0eff[ ind ]
                        except KeyError:
                            print("Missing key %s when index = %s "%( ind, index ) )
                            print(f"n_atom_nohyd[{index}] : {n_atom_nohyd[index]}")
                            print(f"n_atom__type[{index}] : {n_atom__type[index]}")
                            print(f"Atom in Top: {list(pdb.topology.atoms())[index]}")
                            raise(KeyError)        
                # get E0
                n_atom__charge = {}
                # load nonbondedforce
                forces = { force.__class__.__name__ : force for force in system.getForces() }
                nonbonded_force = forces['NonbondedForce']
                for index in range(nonbonded_force.getNumParticles()):
                        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(index)
                        n_atom__charge[index] = charge/charge.unit
                
                #n_atom_nohyd__E0 = get_E0(n_atom__xyz, n_atom__charge, n_atom_nohyd, ion_index)
                
                #get induced dipole
                n_atom_nohyd__dipole, n_atom_nohyd__E0 = get_induced_dipole(n_atom__xyz, n_atom_nohyd__alpha0eff, n_atom__charge, n_atom_nohyd, ion_index)
                
                
                for index in n_atom_nohyd:
                        #le+=1
                        pol_force.addParticle( index, [ n_atom_nohyd__dipole[index][0], n_atom_nohyd__dipole[index][1], n_atom_nohyd__dipole[index][2], n_atom_nohyd__E0[index][0], n_atom_nohyd__E0[index][1], n_atom_nohyd__E0[index][2] ] )
                        
                system.addForce(pol_force)      
        

print("Getting ready for simulation")        
#platform = Platform.getPlatformByName('CUDA')
#prop = dict(CudaPrecision='mixed')    

fgrps=forcegroupify(system)
        
#integrator = VerletIntegrator(0.002*picoseconds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, time_step)


simulation = Simulation(
        modeller.topology, 
        system, 
        integrator
#        platform,
#        prop
)

simulation.context.setPositions(modeller.positions)

print("Minimizing")

simulation.minimizeEnergy()
simulation.reporters.append(
        StateDataReporter(sys.stdout, 10, step=True, time=True, speed=True)
)
simulation.reporters.append(DCDReporter('output.dcd', trajectory_output_frequency))
simulation.reporters.append(StateDataReporter('data.csv', potential_output_frequency, time=True,
        kineticEnergy=True, potentialEnergy=True))
#simulation.reporters.append(PDBReporter('output.pdb', 500))
simulation.reporters.append(CheckpointReporter('checkpnt.chk', 50000))

print("Simulating...")

with open('time_out.dat', 'w') as timeout:
        headers="Step SimTime charge_trans".split() 
        head_form = "{:10s} "*len(headers)
        timeout.write(head_form.format(*headers) )
        timeout.write('\n')
        
        for sp in range(n_steps):
                simulation.step(1)
                #print(simulation.currentStep)
                if sp/10 % 1 != 0: continue
                #print('sp',sp)
                timeout.write("%10d "%sp)
                start = time.time()
                
                # get position
                state = simulation.context.getState(getEnergy=True,getPositions=True) # getForces =True)
                position = state.getPositions()
                index__position = {}
                a=-1
                for i in position:
                    
                        a += 1
                        index__position[a] = [i[0]/nanometer * 10, i[1]/nanometer * 10, i[2]/nanometer * 10]   # nm > angstrom
                        
                # get CN
                cation_pairs = []               
                for i in polar_index:
                        if i < ion_index:
                                cation_pairs.append((i, ion_index))
                        else:
                                cation_pairs.append((ion_index, i))
                            
                CN = get_CN(index__position, cation_pairs, type__zero_transfer_distance)
                
                # update new charges
                # load nonbondedforce
                forces = { force.__class__.__name__ : force for force in system.getForces() }
                nbforce = forces['NonbondedForce']
                
                # charge on L
                #start = time.time()
                total_transq = 0
                index__charge = {}
                for index in range(nbforce.getNumParticles()):
                        
                        # charge
                        [charge, sigma, epsilon] = [index__charge_ini[index], 0.0, 0.0]
                        #print(index, charge)
        
                        #if atom_index < 424:
                        if n_atom__type[index] in type__ChargeTransfer_parameters.keys():
                                r = get_distance(index, ion_index, index__position)
                                print('index:  ', float(index), '  position:  ', float(index__position[index][0]), float(index__position[index][1]), float(index__position[index][2]))
                                # only consider r <= zero_transfer_distance 
                                if r < type__zero_transfer_distance[n_atom__type[index]] :
                                        
                                        if (CN == 1) or (CN_factor == 'No_Charge_transfer'):
                                                transq = (type__ChargeTransfer_parameters[n_atom__type[index]][0] * r + type__ChargeTransfer_parameters[n_atom__type[index]][1]) * (math.sqrt(f15))
                                        else:
                                                transq = ( type__ChargeTransfer_parameters[n_atom__type[index]][0] * r + type__ChargeTransfer_parameters[n_atom__type[index]][1] ) * (math.sqrt(f15)) / math.pow( CN, 1.0/float(CN_factor) ) 
                                        
                                        total_transq += transq
                                        charge_new = ( charge/charge.unit + transq)* charge.unit
                                        
                                        nbforce.setParticleParameters(index, charge_new, sigma, epsilon)
                                else:
                                        charge_new = charge
                                
                        else:
                                charge_new = charge
                        index__charge[index]=charge_new/charge.unit
                print('index:  ', float(ion_index), '  position:  ', float(index__position[ion_index][0]), float(index__position[ion_index][1]), float(index__position[ion_index][2]))
                
                # charge on ion
                [charge, sigma, epsilon] = [index__charge_ini[ion_index], 0.0, 0.0]
                charge_ion = ( charge/charge.unit - total_transq ) * charge.unit
                nbforce.setParticleParameters(ion_index, charge_ion, sigma, epsilon )
                index__charge[ion_index]=charge_ion/charge.unit
                
                # write charge transfered to ion
                print('step: ',sp, 'CT: ',total_transq)
                #print(index__charge)
                
                # set exceptions
                for i in range(nbforce.getNumExceptions()):
                        (p1, p2, q, sig, eps) = nbforce.getExceptionParameters(i)
                        if q._value != 0.0:
                                q_new = (index__charge[p1]*charge.unit) * (index__charge[p2]*charge.unit) * f14
                                nbforce.setExceptionParameters(i, p1, p2, q_new, sig, eps)
                
                nbforce.updateParametersInContext(simulation.context)
                
                ## polarization
                #get induced dipole
                n_atom_nohyd__dipole, n_atom_nohyd__E0 = get_induced_dipole(index__position, n_atom_nohyd__alpha0eff, index__charge, n_atom_nohyd, ion_index)
                
                cmforce = forces['CustomExternalForce']
                
                ind=-1
                for index in n_atom_nohyd:
                        ind+=1
                        cmforce.setParticleParameters( ind, index, [n_atom_nohyd__dipole[index][0], n_atom_nohyd__dipole[index][1], n_atom_nohyd__dipole[index][2], n_atom_nohyd__E0[index][0], n_atom_nohyd__E0[index][1], n_atom_nohyd__E0[index][2]] )
                
                #for i in range(system.getNumForces()):
                #        force = system.getForce(i)
                #        force.setForceGroup(i)	
                #for i in range(system.getNumForces()):
                #        force = system.getForce(i)
                #        energy = simulation.context.getState(getEnergy=True, groups=1<<i).getPotentialEnergy()
                        #energy__decomposition[force.__class__.__name__ ] = energy
                        #if force.__class__.__name__ == 'CustomBondForce':
                        #        timeout.write("%10.3f "%(energy/4.184/kilojoules_per_mole))
                        #print(force.__class__.__name__, energy/4.184/kilojoules_per_mole )
                #        print( energy/4.184/kilojoules_per_mole, end='    ' ),
                
                
                cmforce.updateParametersInContext(simulation.context)
                
                tt= getEnergyDecomposition(simulation.context, fgrps)
                for idd in tt.keys():
                    print(idd,tt[idd]) 
                print('\n')
                print('################\n')
                
#simulation.step(20000)
simulation.saveCheckpoint('state.chk')
sys.exit(0)
# Write restart file
if not (args.orst or args.ochk): args.orst = 'output.rst'
if args.orst:
    state = simulation.context.getState( getPositions=True, getVelocities=True )
    with open(args.orst, 'w') as f:
        f.write(XmlSerializer.serialize(state))
if args.ochk:
    with open(args.ochk, 'wb') as f:
        f.write(simulation.context.createCheckpoint())
if args.opdb:
    crd = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(top.topology, crd, open(args.opdb, 'w'))



sys.exit(0)

# time cost
# 100 steps : 2.336803674697876 s
# 1000 steps : 136.5139036178589 s

# updata parameter every 20 step, 33 atoms
# 1000 steps : 6.988223314285278 s
# 2000 steps : 26.175992965698242 s
# 3000 steps : 57.635138511657715 s

#savedStdout = sys.stdout
#sys.stdout = f
#for i in range(system.getNumForces()):
#       force = system.getForce(i)
#       force.setForceGroup(i)
#for i in range(system.getNumForces()):
#       force = system.getForce(i)
#       energy = simulation.context.getState(getEnergy=True, groups=1<<i).getPotentialEnergy()
#       #energy__decomposition[force.__class__.__name__ ] = energy
#       if force.__class__.__name__ == 'CustomBondForce':
#               print(file.rsplit('/',2)[1], energy/4.184/kilojoules_per_mole )
#print(file.split('.')[0], state.getPotentialEnergy()/4.184/kilojoules_per_mole)
#forces = { force.__class__.__name__ : force for force in system.getForces() }
#bondforce = forces['HarmonicBondForce']
#list_12_interacts = []
#for index in range(bondforce.getNumBonds()):
#       [particle1, particle2, r0, Kb]= bondforce.getBondParameters(index)
#       if int(particle1) <= int(particle2):
#               atuple = (int(particle1),int(particle2))
#       else:
#               atuple = (int(particle2),int(particle1))
#       list_12_interacts.append( atuple )
##print(list_12_interacts)
#torforce = forces['PeriodicTorsionForce']
#a = 0
#list_imp1234_interacts = []
#for index in range(torforce.getNumTorsions()):
#       [ particle1, particle2, particle3, particle4, perio, phase, k] = torforce.getTorsionParameters(index)
#       a += 1
#       if particle1 <= particle2:
#               atuple1 = (particle1, particle2)
#       else:
#               atuple1 = (particle2, particle1)
#       if particle2 <= particle3:
#               atuple2 = (particle2, particle3)
#       else:
#               atuple2 = (particle3, particle2)
#       if particle3 <= particle4:
#               atuple3 = (particle3, particle4)
#       else:
#               atuple3 = (particle4, particle3)
#       if (atuple1 or atuple2 or atuple3) not in list_12_interacts:            
#               atuple = (int(particle1), int(particle2), int(particle3), int(particle4))               
#               list_imp1234_interacts.append( atuple )
#print(list_imp1234_interacts)



#print(a)
#print(file.split('.')[0])
#energy__decomposition = {}
#for i in range(system.getNumForces()):
#       force = system.getForce(i)
#       energy = simulation.context.getState(getEnergy=True, groups=1<<i).getPotentialEnergy()
#       #energy__decomposition[force.__class__.__name__ ] = energy
#       print(force.__class__.__name__, energy/4.184/kilojoules_per_mole)
#print('FF: ', state.getPotentialEnergy()/4.184/kilojoules_per_mole)

savedStdout = sys.stdout
sys.stdout = f
print(file.rsplit('/',2)[1], state.getPotentialEnergy()/4.184/kilojoules_per_mole) 
##print(file.split('.')[0], state.getForces(asNumpy=True)/4.184)   Force

#os.system('sort energies_openmm_pre.kj > energies_openmm.kj')
