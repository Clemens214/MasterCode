# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 08:37:09 2024

@author: cleme
"""
import random
import numpy as np
from mpmath import zeta
import scipy.linalg as lin
import matplotlib.pyplot as plt

def randomNum (eigenenergy, size):
    values = list()
    for i in range(size):
        randVal = random.random()*2 - 1
        values.append(randVal*eigenenergy)
    return values

def makeHamiltonian (eigenenergies, hopping, size=None, alignment='top left'):
    alignment = alignment.lower()
    if size == None:
        size = len(eigenenergies)
    maxSize = max(len(eigenenergies), size)
    
    if 'top' in alignment:
        startRow = 0
    elif 'bottom' in alignment:
        startRow = maxSize - size
    if 'left' in alignment:
        startColumn = 0
    elif 'right' in alignment:
        startColumn = maxSize - size
    
    H = np.zeros((maxSize, maxSize), dtype=complex)
    for row in np.arange(len(eigenenergies)):
        for column in np.arange(len(eigenenergies)):
            if row == column:
                H[row, column] = eigenenergies[row]
    for row in np.arange(size):
        for column in np.arange(size):
            if column == row+1 or column == row-1:
                H[startRow+row, startColumn+column] = hopping
    return np.array(H)

def prepareEM (sample, hopping, lengthSample, *fermiParameters, hoppingLead=None, hoppingInter=None):
    lengthTotal, maxVal, decay, offset = fermiParameters
    if hoppingLead == None:
        hoppingLead = hopping
    if hoppingInter == None:
        hoppingInter = hopping
    
    def fermifuncLeft (sites, maxVal, decay, offset):
        values = [maxVal/(1+np.exp(decay*(i+1-offset))) for i in sites]
        return np.array(values)

    def fermifuncRight (sites, maxVal, decay, offset, lengthTotal):
        values = [maxVal/(1+np.exp(decay*(lengthTotal-i-offset))) for i in sites]
        return np.array(values)
    
    def combineH (center, left, right, hoppingInter, length):
        sizeMax = max(len(center), len(left), len(right))
        arrayNew = np.zeros((sizeMax, sizeMax), dtype=complex)
        for row in np.arange(length-1, sizeMax-length):
            for column in np.arange(length-1, sizeMax-length):
                if row-length >= 0 and column-length >= 0:
                    arrayNew[row, column] = center[row-length, column-length]
        values = [length, sizeMax-length]
        for row in np.arange(sizeMax):
            for column in np.arange(sizeMax):
                newFactor = 0
                if row == values[0]-1 and column == values[0]:
                    newFactor = hoppingInter
                elif row == values[0] and column == values[0]-1:
                    newFactor = hoppingInter
                if row == values[1]-1 and column == values[1]:
                    newFactor = hoppingInter
                if row == values[1] and column == values[1]-1:
                    newFactor = hoppingInter
                arrayNew[row, column] += newFactor + left[row][column] + right[row][column]
        return arrayNew
    
    # compute the fermi functions
    sites = range(lengthTotal)
    fermifuncLeft = fermifuncLeft(sites, maxVal, decay, offset)
    fermifuncRight = fermifuncRight(sites, maxVal, decay, offset, lengthTotal)
    
    # generate the Hamiltonians of the leads
    length = int((lengthTotal-lengthSample)/2)
    sigmaL = makeHamiltonian(1j*fermifuncLeft, hoppingLead, size=length, alignment = 'top left')
    sigmaR = makeHamiltonian(1j*fermifuncRight, hoppingLead, size=length, alignment = 'bottom right')
    
    # compute the coupling strengths
    gammaL = 1j*(sigmaL - sigmaL.conj().T)
    gammaR = 1j*(sigmaR - sigmaR.conj().T)
    
    #generate the Hamiltonian of the total system
    totalSystem = combineH(sample, sigmaL, sigmaR, hoppingInter, length)
    return totalSystem, gammaL, gammaR

def normalize (leftEVs, rightEVs):
    for i in range(len(leftEVs)):
        leftEV = leftEVs[:,i].conj().T
        rightEV = rightEVs[:,i]
        value = np.dot(leftEV, rightEV)
        
        leftEVs[:,i] = leftEVs[:,i]/np.sqrt(value).conj()
        rightEVs[:,i] = rightEVs[:,i]/np.sqrt(value)
        print(np.dot(leftEVs[:,i].conj().T, rightEVs[:,i]))
    return leftEVs, rightEVs

def TestEV (Matrix, Eigenvals, leftEVs, rightEVs, eta = 1e-10):
    TestLeft = list()
    DiffLeft = list()
    TestRight = list()
    DiffRight = list()
    # [ np.allclose(np.dot(vl[:,i].T, T), w[i]*vl[:,i].T) for i in range(3) ]
    for i in range(len(Eigenvals)):
        # check the left Eigenvectors
        leftEV = leftEVs[:,i].conj().T
        MatrixMultLeft = np.dot(leftEV, Matrix)
        EigMultLeft = Eigenvals[i] * leftEV
        TestLeft.append(np.allclose(MatrixMultLeft, EigMultLeft, atol=eta))
        DiffLeft.append(MatrixMultLeft-EigMultLeft)
        # check the right Eigenvectors
        rightEV = rightEVs[:,i]
        MatrixMultRight = np.dot(Matrix, rightEV)
        EigMultRight = Eigenvals[i] * rightEV
        TestRight.append(np.allclose(MatrixMultRight, EigMultRight, atol=eta))
        DiffRight.append(MatrixMultRight-EigMultRight)
    return TestLeft, TestRight, DiffLeft, DiffRight

def eigenvectors (totalSystem):
    Eigenvals, leftEVs, rightEVs = lin.eig(totalSystem, left=True, right=True)
    
    leftEVs, rightEVs = normalize(leftEVs, rightEVs)
    
    TestLeft, TestRight, DiffLeft, DiffRight = TestEV(totalSystem, Eigenvals, leftEVs, rightEVs)
    if False in TestLeft or False in TestRight:
        print("Non matching Eigenvector found.")
    
    ProductAll = list()
    for i in range(len(Eigenvals)):
        leftEigVec = leftEVs[:,i].conj().T
        leftEigVec = leftEigVec.reshape(1, -1)
        rightEigVec = rightEVs[:,i]
        rightEigVec = rightEigVec.reshape(-1, 1)
        ProductAll.append(np.matmul(rightEigVec, leftEigVec))
    Product = sum(ProductAll)
    
    offDiagMax = 0
    Diag = list()
    for i in range(len(Product)):
        Diag.append(abs(Product[i][i]))
        for j in range(len(Product)):
            if i!=j and offDiagMax < abs(Product[i][j]):
                offDiagMax = abs(Product[i][j])
    
    Control = [Product, min(Diag), offDiagMax]
    return Eigenvals, leftEVs, rightEVs, Control

def HurwitzZeta (x, beta):
    factor = np.sign(x.imag)/(2*np.pi*1j*beta)
    value = -1*factor * zeta(2, q=0.5+factor*x)
    #return value
    return complex(value.real, value.imag)

def particleElement (eig1, eig2, T):
    beta = 1/T
    factor = 1/(eig1 -eig2)
    value = factor*HurwitzZeta(eig1, beta) - factor*HurwitzZeta(eig2, beta)
    return value

def energyElement (eig1, eig2, T):
    beta = 1/T
    factor = 1/(eig1 -eig2)
    value = eig1*factor*HurwitzZeta(eig1, beta) - eig2*factor*HurwitzZeta(eig2, beta)
    return value

def currentElement (stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp, chemPot):
    print('starting calculation of the current element.')
    value = 0
    for i in range(len(leftEVs)):
        for j in range(len(leftEVs)):
            # get the normal left and right Eigenvectors
            EigVal = Eigenvals[i]
            leftEV = leftEVs[:,i].conj().T
            leftEV = leftEV.reshape(1, -1)
            rightEV = rightEVs[:,i]
            rightEV = rightEV.reshape(-1, 1)
            
            # get the daggered Eigenvectors
            EigValDagger = Eigenvals[j].conj()
            leftEVdagger = leftEVs[:,j]
            leftEVdagger = leftEVdagger.reshape(-1, 1)
            rightEVdagger = rightEVs[:,j].conj().T
            rightEVdagger = rightEVdagger.reshape(1, -1)
            
            #compute the matrix element for chosen i and j
            ProductLeft = np.matmul(stateLeft, rightEV)
            ProductRight = np.matmul(rightEVdagger, stateRight)
            ProductMid = np.matmul(np.matmul(leftEV, midFactor), leftEVdagger)
            
            Product = np.matmul(np.matmul(ProductLeft, ProductMid), ProductRight)
            
            #compute the additional matrix element
            factor = energyElement(EigVal, EigValDagger, Temp) - chemPot*particleElement(EigVal, EigValDagger, Temp)
            
            value += Product*factor
    return value[0][0]

def current (totalSystem, gammaL, gammaR, Temp, chemPot, value=0):
    print('starting current calculation.')
    # define the states to be used for calculating the current
    stateLeft = totalSystem[:,value].conj().T
    stateLeft = stateLeft.reshape(1, -1)
    stateRight = totalSystem[:,value+1]
    stateRight = stateRight.reshape(-1, 1)
    
    #compute the Eigenvectors and the Eigenvalues of the system
    Eigenvals, leftEVs, rightEVs, ControlEV = eigenvectors(totalSystem)
    
    #compute the entropy current
    midFactor = gammaL - gammaR
    value = currentElement (stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp, chemPot)
    return value, ControlEV

#%% Variables
# variables for the sample
lengthSample = 96
eigenenergy = 0.1
hopping = 1

# variables for the leads
lengthTotal = 256
maxVal = 1
decay = 0.2 #0.3
offset = 32 #32
#offset shoeld be at most halt the length of the leads
hoppingLead = hopping
hoppingInter = hopping

# variables for the calculation of the current
Temp = 1e-10
chemPot = 1

#%% Calculation
#%%% compute the Hamiltonian of the Sample
sample = makeHamiltonian(randomNum(eigenenergy, lengthSample), hopping)

# prepare the Hamiltonian of the extended molecule formalism
length = int((lengthTotal-lengthSample)/2)
fermiParameters = [lengthTotal, maxVal, decay, offset]

#%%% compute the currents
currentsEM = list()

# current for the extended molecule formalism
totalSystemEM, gammaL_EM, gammaR_EM = prepareEM(sample, hopping, lengthSample, *fermiParameters)
currentEM, ControlEV_EM = current (totalSystemEM, gammaL_EM, gammaR_EM, Temp, chemPot, value=length)
currentsEM.append(currentEM)
print('done: Extended Molecule:', currentEM)

if False:
    def prepareSI (omega, Hamiltonian, size, hopping, eta=1e-10, hoppingInter=None, hoppingBath=None):
        sigmaL = np.zeros((size, size), dtype=complex)
        sigmaR = np.zeros((size, size), dtype=complex)
        if hoppingBath == None:
            hoppingBath = hopping
        if hoppingInter == None:
            hoppingInter = hopping
        
        G = (omega+1j*eta)/(2*hoppingBath**2)*(1-np.sqrt(1-(4*hoppingBath**2)/(omega+1j*eta)**2))
        
        # generate the Hamiltonians of the leads
        sigmaL[0,0] = hoppingInter**2 * G
        sigmaR[size-1,size-1] =hoppingInter**2 * G
        
        # compute the coupling strengths
        gammaL = 1j*(sigmaL - sigmaL.conj().T);
        gammaR = 1j*(sigmaR - sigmaR.conj().T);
    
        #generate the Hamiltonian of the total system
        totalSystem = Hamiltonian + sigmaL + sigmaR
        return totalSystem, gammaL, gammaR
    
    # variables for the final calculation
    step = 1
    omegas = np.arange(-2, 2+step, step)
    currentsSI = list()
    for omega in omegas:
        # current for the semi-infinite leads
        totalSystemSI, gammaL_SI, gammaR_SI = prepareSI(omega, sample, lengthSample, hopping)
        currentSI, ControlEV_SI = current (totalSystemSI, gammaL_SI, gammaR_SI, Temp, chemPot, value=0)
        currentsSI.append(currentSI)
        print('done: semi-infinite leads:', currentSI, 'for omega:', omega)