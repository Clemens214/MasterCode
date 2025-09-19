# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 08:37:09 2024

@author: cleme
"""
import random
import numpy as np
from mpmath import mp
from mpmath import zeta, sqrt
from scipy.linalg import eig
from timeit import default_timer as timer
mp.dps = 20
start = timer()

def randomNum (eigenenergy, size):
    values = list()
    for i in range(size):
        randVal = mp.mpf(random.random()*2 - 1)
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
    
    H = mp.matrix(maxSize)
    for row in range(len(eigenenergies)):
        for column in range(len(eigenenergies)):
            if row == column:
                H[row, column] = eigenenergies[row]
    for row in range(size):
        for column in range(size):
            if column == row+1 or column == row-1:
                H[startRow+row, startColumn+column] = hopping
    return H

def prepareEM (sample, hopping, lengthSample, *fermiParameters, hoppingLead=None, hoppingInter=None):
    lengthTotal, maxVal, decay, offset = fermiParameters
    if hoppingLead == None:
        hoppingLead = hopping
    if hoppingInter == None:
        hoppingInter = hopping
    
    def fermifuncLeft (sites, maxVal, decay, offset):
        values = [mp.mpf(maxVal/(1+mp.exp(decay*(i+1-offset)))) for i in sites]
        return mp.matrix(values)

    def fermifuncRight (sites, maxVal, decay, offset, lengthTotal):
        values = [mp.mpf(maxVal/(1+mp.exp(decay*(lengthTotal-i-offset)))) for i in sites]
        return mp.matrix(values)
    
    def combineH (center, left, right, hoppingInter, length):
        sizeMax = max(len(center), len(left), len(right))
        arrayNew = mp.matrix(sizeMax)
        for row in range(length-1, sizeMax-length):
            for column in range(length-1, sizeMax-length):
                if row-length >= 0 and column-length >= 0:
                    arrayNew[row, column] = center[row-length, column-length]
        values = [length, sizeMax-length]
        for row in range(sizeMax):
            for column in range(sizeMax):
                newFactor = 0
                if row == values[0]-1 and column == values[0]:
                    newFactor = hoppingInter
                elif row == values[0] and column == values[0]-1:
                    newFactor = hoppingInter
                if row == values[1]-1 and column == values[1]:
                    newFactor = hoppingInter
                if row == values[1] and column == values[1]-1:
                    newFactor = hoppingInter
                arrayNew[row, column] += newFactor + left[row, column] + right[row, column]
        return arrayNew
    
    # compute the fermi functions
    sites = range(lengthTotal)
    fermifuncLeft = fermifuncLeft(sites, maxVal, decay, offset)
    fermifuncRight = fermifuncRight(sites, maxVal, decay, offset, lengthTotal)
    
    # generate the Hamiltonians of the leads
    length = int((lengthTotal-lengthSample)/2)
    sigmaL = makeHamiltonian(mp.mpc(1j)*fermifuncLeft, hoppingLead, size=length, alignment = 'top left')
    sigmaR = makeHamiltonian(1j*fermifuncRight, hoppingLead, size=length, alignment = 'bottom right')
    
    # compute the coupling strengths
    gammaL = 1j*(sigmaL - sigmaL.conjugate().T)
    gammaR = 1j*(sigmaR - sigmaR.conjugate().T)
    
    #generate the Hamiltonian of the total system
    totalSystem = combineH(sample, sigmaL, sigmaR, hoppingInter, length)
    return totalSystem, gammaL, gammaR

def normalize (leftEVs, rightEVs):
    for i in range(len(leftEVs)):
        leftEV = leftEVs[i,:]
        rightEV = rightEVs[:,i]
        value = leftEV * rightEV
        
        squareRoot = sqrt(value[0,0])
        leftEVs[i,:] = leftEVs[i,:]/squareRoot#.conjugate()
        rightEVs[:,i] = rightEVs[:,i]/squareRoot
        #print(leftEVs[:,i].conjugate().T * rightEVs[:,i])
    return leftEVs, rightEVs

def TestEV (Matrix, Eigenvals, leftEVs, rightEVs, eta=1e-10):
    TestLeft = list()
    DiffLeft = list()
    TestRight = list()
    DiffRight = list()
    # [ np.allclose(np.dot(vl[:,i].T, T), w[i]*vl[:,i].T) for i in range(3) ]
    for i in range(len(Eigenvals)):
        # check the left Eigenvectors
        leftEV = leftEVs[i,:]
        MatrixMultLeft = np.array(leftEV * Matrix)
        MatrixMultLeft = np.array([complex(val.real, val.imag) for val in MatrixMultLeft])
        EigMultLeft = np.array(Eigenvals[i] * leftEV)
        EigMultLeft = np.array([complex(val.real, val.imag) for val in EigMultLeft])
        TestLeft.append(np.allclose(MatrixMultLeft, EigMultLeft, atol=eta))
        DiffLeft.append(MatrixMultLeft-EigMultLeft)
        # check the right Eigenvectors
        rightEV = rightEVs[:,i]
        MatrixMultRight = Matrix * rightEV
        MatrixMultRight = np.array([complex(val.real, val.imag) for val in MatrixMultRight])
        EigMultRight = Eigenvals[i] * rightEV
        EigMultRight = np.array([complex(val.real, val.imag) for val in EigMultRight])
        TestRight.append(np.allclose(MatrixMultRight, EigMultRight, atol=eta))
        DiffRight.append(MatrixMultRight-EigMultRight)
    return TestLeft, TestRight, np.array(DiffLeft), np.array(DiffRight)

def eigenvectors (totalSystem, accurate = True):
    print('Starting calculation of the Eigenvectors.')
    startEV = timer()
    if accurate == True:
        Eigenvals, leftEVs, rightEVs = mp.eig(totalSystem, left=True, right=True)
    else:
        Eigenvals, leftEVs, rightEVs = eig(totalSystem, left=True, right=True)
    endEV = timer()
    print('Time spent calculating Eigenvectors:', (endEV - startEV)/60)
    
    print('Starting normalization of Eigenvectors.')
    leftEVs, rightEVs = normalize(leftEVs, rightEVs)
    print('Finished normalizing Eigenvectors')
    
    # check whether the eigenvectors are correct
    print('Checking Eigenvectors.')
    startCheck = timer()
    MatchLeft, MatchRight, DiffLeft, DiffRight = TestEV(totalSystem, Eigenvals, leftEVs, rightEVs)
    if False in MatchLeft or False in MatchRight:
        print("Non matching Eigenvector found.")
    else:
        print("All Eigenvectors match.")
    TestLeft = [MatchLeft, DiffLeft]
    TestRight = [MatchRight, DiffRight]
    
    ProductAll = list()
    for i in range(len(Eigenvals)):
        leftEV = leftEVs[i,:]
        rightEV = rightEVs[:,i]
        Mult = rightEV * leftEV
        MultNP = [[complex(Mult[a,b].real, Mult[a,b].imag) for b in range(Mult.cols)] for a in range(Mult.rows)]
        ProductAll.append(np.array(MultNP))
    Product = sum(ProductAll)
    
    offDiagMax = 0
    Diag = list()
    DiagImag = list()
    for i in range(len(Product)):
        Diag.append(abs(Product[i,i]))
        DiagImag.append(Product[i,i].imag)
        for j in range(len(Product)):
            if i!=j and offDiagMax < abs(Product[i,j]):
                offDiagMax = abs(Product[i,j])
    
    Control = [np.array(Product), TestLeft, TestRight, [max(Diag), min(Diag)], max(DiagImag), offDiagMax]
    endCheck = timer()
    print('Time spent checking Eigenvectors:', (endCheck - startCheck)/60)
    return Eigenvals, leftEVs, rightEVs, Control

def HurwitzZeta (x, beta):
    factor = mp.sign(x.imag)/(2*mp.pi*1j*beta)
    value = -1*factor * zeta(2, q=0.5+factor*x)
    return value

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
    value = 0
    for i in range(len(leftEVs)):
        for j in range(len(leftEVs)):
            # get the normal left and right Eigenvectors
            EigVal = Eigenvals[i]
            leftEV = leftEVs[i,:]
            rightEV = rightEVs[:,i]
            
            # get the daggered Eigenvectors
            EigValDagger = Eigenvals[j].conjugate()
            leftEVdagger = leftEV.conjugate().T
            rightEVdagger = rightEV.conjugate().T
            
            #compute the matrix element for chosen i and j
            ProductLeft = stateLeft * rightEV
            ProductRight = rightEVdagger * stateRight
            ProductMid = leftEV * midFactor * leftEVdagger
            
            Product = ProductLeft * ProductMid * ProductRight
            
            #compute the additional matrix element
            factor = energyElement(EigVal, EigValDagger, Temp) - chemPot*particleElement(EigVal, EigValDagger, Temp)
            
            value += Product*factor
    return value

def current (totalSystem, gammaL, gammaR, Temp, chemPot, value=0):
    # define the states to be used for calculating the current
    stateLeft = totalSystem[:,value].conjugate().T
    stateRight = totalSystem[:,value+1]
    
    #compute the Eigenvectors and the Eigenvalues of the system
    Eigenvals, leftEVs, rightEVs, ControlEV = eigenvectors(totalSystem)
    
    #compute the entropy current
    midFactor = gammaL - gammaR
    startCurrent = timer()
    print('Starting calculation of the current.')
    value = currentElement (stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp, chemPot)
    endCurrent = timer()
    print('Time spent calculating current:', (endCurrent - startCurrent)/60)
    return value[0,0], ControlEV

#%% Variables
# variables for the sample
lengthSample = 96 #96
#normal lengthSample: 96
eigenenergy = 0.1
hopping = 1

# variables for the leads
lengthTotal = 256 #256
#normal lengthTotal: 256
lengthLead = int((lengthTotal-lengthSample)/2)
print('Sample:', lengthSample, '\n',
      'Leads:', lengthLead, '\n', 
      'Total:', lengthTotal)
maxVal = 1
decay = 0.2 #0.3
offset = 32 #32
#offset should be at most half the length of the leads; normal: 32
hoppingLead = hopping
hoppingInter = hopping

# variables for the calculation of the current
Temp = 1e-10
chemPot = 1

#%% Calculation
#%%% compute the Hamiltonian of the Sample
sample = makeHamiltonian(randomNum(eigenenergy, lengthSample), hopping)

# prepare the Hamiltonian of the extended molecule formalism
fermiParameters = [lengthTotal, maxVal, decay, offset]

#%%% compute the currents
# current for the extended molecule formalism
print('Preparing the Extended Molecule Hamiltonian.')
totalSystemEM, gammaL_EM, gammaR_EM = prepareEM(sample, hopping, lengthSample, *fermiParameters)
print('Starting the calculation.')
currentEM, ControlEV_EM = current(totalSystemEM, gammaL_EM, gammaR_EM, Temp, chemPot, value=lengthLead)
print('done: Extended Molecule:', currentEM)

# print the elapsed time
end = timer()
print('Elapsed time:', (end - start)/60)