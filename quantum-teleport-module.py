"""Quantum Teleport Module

This module contains definitions for Quantum Teleport simulation on QuTip.
It allows the user to implement the Standart Teleport Protocol and the
Realistic Teleport Protocol by choosing which noise and decoherence will be
applied to the Quantum Channel and the state to be teleported, respectively.

Also allows plotting of resulting Bloch Sphere Deformation and the resulting
quantum fidelity distribution histrogram.

This script requires that qutip be installed whitin the Python environment you
are running the script on. You need to import numpy library and matplotlib,random,
qutip modules. For 3D plots, import Axes3D from "mpl_toolkits.mplot3d". 

"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import qutip as qt
import random
from mpl_toolkits.mplot3d import Axes3D 

def state_drawing():
    """ gets a random quantum state for testing teleport
        
        Paramaters
        -----------
        None
        
        Returns
        --------
        Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    """
    
    def draw(start, end):
        return random.random() * (end - start) + start
    th = draw(0,np.pi)
    ph = draw(0,2*np.pi)
    return ket2dm(np.cos(th/2)*basis(2,0)+np.exp(1j*ph)*np.sin(th/2)*basis(2,1))
"""
options for noise and decoherence maps

"""
def bit_flip(p):
    """Bit-Flip Channel's list of Kraus operators.
    
    Parameters
    -----------
    p: number in 0 to 1 range
        probability of flipping the state
        
    Returns
    --------
    list: list of quantum operators
        [
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True]
        ]
    """
    return [np.sqrt(p)*qeye(2),np.sqrt(1-p)*sigmax()]

def phase_flip(p):
    """Phase-Flip Channel's list of Kraus operators.
    
    Parameters
    -----------
    p: number in 0 to 1 range
        probability of flipping the state
        
    Returns
    --------
    list: list of quantum operators
        [
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True]
        ]
    """
    return [np.sqrt(p)*qeye(2),np.sqrt(1-p)*sigmaz()]

def bitphase_flip(p):
    """BitPhase-Flip Channel's list of Kraus operators.
    
    Parameters
    -----------
    p: number in 0 to 1 range
        probability of flipping the state
        
    Returns
    --------
    list: list of quantum operators
        [
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True]
        ]
    """
    return [np.sqrt(p)*qeye(2),np.sqrt(1-p)*sigmay()]

def depolarizing_channel(p):
    """Depolarizing Channel's list of Kraus operators.
    
    Parameters
    -----------
    p: number in 0 to 1 range
        probability of flipping the state
        
    Returns
    --------
    list: list of quantum operators
        [
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True]
        ]
    """
    return [np.sqrt(1-3*p/4)*qeye(2),(np.sqrt(p)/2)*sigmax(),(np.sqrt(p)/2)*sigmay(),(np.sqrt(p)/2)*sigmaz()]

def amplitude_damping(gamma):
    """Amplitude Damping Channel's list of Kraus operators.
    
    Parameters
    -----------
    Gamma: number in 0 to 1 range
        Amplitude damping coefficient
        
    Returns
    --------
    list: list of quantum operators
        [
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True]
        ]
    """
    return [Qobj([[1,0],[0,np.sqrt(1-gamma)]]),Qobj([[0,np.sqrt(gamma)],[0,0]])]

def general_amplitude_damping(gamma,p):
    """General Amplitude Damping Channel's list of Kraus operators.
    
    Parameters
    -----------
    p: number in 0 to 1 range
        probability of flipping the state
    Gamma: number in 0 to 1 range
        Amplitude damping coefficient
        
    Returns
    --------
    list: list of quantum operators
        [
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True]
        ]
    """
    return [np.sqrt(p)*Qobj([[1,0],[0,np.sqrt(1-gamma)]]),np.sqrt(p)*Qobj([[0,np.sqrt(gamma)],[0,0]]),np.sqrt(1-p)*Qobj([[np.sqrt(1-gamma),0],[0,1]])],np.sqrt(1-p)*Qobj([[0,0],[np.sqrt(gamma),0]])

def phase_damping(alpha):
    """Phase Damping Channel's list of Kraus operators.
    
    Parameters
    -----------
    Alpha: number in 0 to 1 range
        Phase damping coefficient
        
    Returns
    --------
    list: list of quantum operators
        [
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True],
        [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True]
        ]
    """
    return [np.sqrt(alpha)*qeye(2),np.sqrt(1-alpha)*sigmaz()]
    
def noise(kraus_list):
    """Applies the Choi Matrix (associated to the chosen noise) to quantum channel, which is by
    definition the singlet state.
    
    Parameters
    -----------
    kraus_list: list
        list of kraus operators. You can choose the built-in channels for noise simulation listed above.
        
    Returns
    ---------
    Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
    """
    return sum([tensor(qeye(2),E)*ket2dm(bell_state('00'))*tensor(qeye(2),E.dag()) for E in kraus_list])

def decoherence(state,kraus_list):
    """Applies the chosen decoherence process to the state to be teleported.
    
    Parameters
    -----------
    kraus_list: list
        list of kraus operators. You can choose the built-in channels for noise simulation listed above.
    state: Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
        state to be teleported
        
    Returns
    ---------
    Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
    """    
    return sum([K*state*K.dag() for K in kraus_list])

def STP_teleport(state):
    """Applies the STP teleport with no noise in quantum channel
    
    Parameters
    ------------
    state: Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
        state to be teleported
    
    Returns
    --------
    Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
        output state
    
    """
    state_channel = tensor(state, ket2dm(bell_state('00')))
    A = [bell_state('00').proj(), bell_state('01').proj(), bell_state('10').proj(),bell_state('11').proj()]
    B = [qeye(2),sigmaz(),sigmax(),1j*sigmay()]
    return sum([(tensor(A[i],B[i])*state_channel*tensor(A[i].dag(),B[i].dag())) for i in [0,1,2,3]]).ptrace(2)

def STP_noisy_teleport(state, noise_list):
    """Applies the STP teleport through noisy quantum channel
    
    Parameters
    ------------
    state: Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
        state to be teleported
    noise_list: list
        list of kraus operators. You can choose the built-in channels for noise simulation listed above.
    
    Returns
    --------
    Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
        output state
    
    """
    state_channel = tensor(state, noise(noise_list))
    A = [bell_state('00').proj(), bell_state('01').proj(), bell_state('10').proj(),bell_state('11').proj()]
    B = [qeye(2),sigmaz(),sigmax(),1j*sigmay()]
    return sum([(tensor(A[i],B[i])*state_channel*tensor(A[i].dag(),B[i].dag())) for i in [0,1,2,3]]).ptrace(2)

def realistic_teleport_protocol(state,noise_list,decoherence_list):
    """Applies the realistic teleport protocol
    
    Parameters
    ------------
    state: Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
        state to be teleported
    noise_list: list
        list of kraus operators. You can choose the built-in channels for noise simulation listed above.
    decoherence_list: list
        list of kraus operators. You can choose the built-in channels for noise simulation listed above.
    
    Returns
    --------
    Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
        output state
    
    """
    state_channel = tensor(decoherence(state,decoherence_list), noise(noise_list))
    A = [bell_state('00').proj(), bell_state('01').proj(), bell_state('10').proj(),bell_state('11').proj()]
    B = [qeye(2),sigmaz(),sigmax(),1j*sigmay()]
    return sum([(tensor(A[i],B[i])*state_channel*tensor(A[i].dag(),B[i].dag())) for i in [0,1,2,3]]).ptrace(2)

def fidelity_histogram(noise_list,decoherence_list):
    """plots a histogram of quantum states fidelity resulting of the realistic teleport protocol with chosen
    noise channel and decoherence channel.
    
    Parameters
    -----------
    noise_list: list
        list of kraus operators. You can choose the built-in channels for noise simulation listed above.
    decoherence_list: list
        list of kraus operators. You can choose the built-in channels for noise simulation listed above.
        
    Returns
    --------
    Plot:
        Histogram plot with 100 bins and figsize=(10,4) and 'font.size':20
    
    """
    psi_in = [ket2dm(np.cos(th/2)*basis(2,0)+np.exp(1j*ph)*np.sin(th/2)*basis(2,1)) for th in np.linspace(0,np.pi,30) for ph in np.linspace(0,2*np.pi,2*30)]
    F_list = [fidelity(psi, realistic_teleport_protocol(psi,noise_list,decoherence_list)) for psi in psi_in]
       
    
    plt.rcParams.update({'font.size':20})

    plt.figure(figsize=(10,4))
    plt.xlabel('Fidelity of quantum states',fontsize=14)
    
    plt.title('Fidelity Distribution ', fontsize=18)
    
    n, bins, _ = plt.hist(F_list, bins=100)
    
"""coordinates transformations for Bloch Sphere's Deformation Plot"""

def psi_to_polar(psi):
    """ gets polar coordinates of quantum state in Bloch Sphere representation
    
    Parameters
    -----------
    state: Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
        
    Returns
    --------
    list: [r,theta,phi]
        list of polar coordinates
    
    """
    r = (np.sqrt(((2*psi-qeye(2))*(2*psi-qeye(2))).tr()/2))
    theta = (np.arccos((2*psi-qeye(2))[0,0]/r)).real
    if (np.sin(theta)<=0.00001): 
        phi = 0
    else:
        phi = (1j*np.log(((2*psi-qeye(2))[0,1])/np.sin(theta))).real
    return [r,theta,phi]

def polar_to_cartesian(P):
    """gets the list of polar coordinates and transforms to cartesian coordinates
    
    Parameters
    -----------
    list: [r,theta,phi]
        list of polar coordinates
    Returns
    --------
    list: [x,y,z]
        list of cartesian coordinates        
    """
    return [P[0]*np.sin(P[1])*np.cos(P[2]),P[0]*np.sin(P[1])*np.sin(P[2]),P[0]*np.cos(P[1])]

def Bloch_Sphere_Contraction(noise_list,decoherence_list):
    """plots the Bloch Sphere deformation resulting from the realistic teleport protocol with chosen
    noise channel and decoherence channel.
    
    Parameters
    -----------
    noise_list: list
        list of kraus operators. You can choose the built-in channels for noise simulation listed above.
    decoherence_list: list
        list of kraus operators. You can choose the built-in channels for noise simulation listed above.
        
    Returns
    --------
    Plot:
        Bloch Sphere Shape
    """
    
    psi_in = [ket2dm(np.cos(th/2)*basis(2,0)+np.exp(1j*ph)*np.sin(th/2)*basis(2,1)) for th in np.linspace(0,np.pi,30) for ph in np.linspace(0,2*np.pi,2*30)]
    psi_out = [realistic_teleport_protocol(psi,noise_list,decoherence_list) for psi in psi_in]
    P = [polar_to_cartesian(psi_to_polar(psi)) for psi in psi_out]
    X = [P[i][0].real for i in range(0,len(P))]
    Y = [P[i][1].real for i in range(0,len(P))]
    Z = [P[i][2].real for i in range(0,len(P))]
   
    
    
    fig = plt.figure()
    ax = plt.axes(projection="3d")


    X = [P[i][0].real for i in range(0,len(P))]
    Y = [P[i][1].real for i in range(0,len(P))]
    Z = [P[i][2].real for i in range(0,len(P))]

    ax.scatter(X, Y, Z, c='r', marker='o')

    plt.title('Bloch Sphere', fontsize=18)
    
    plt.show()