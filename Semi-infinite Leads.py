# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 15:12:46 2024

@author: cleme
"""
import numpy as np
import matplotlib.pyplot as plt
import random

def randomNum (eigenenergy, size):
    values = list()
    for i in range(size):
        randVal = random.random()*2 - 1
        values.append(randVal*eigenenergy)
    return values

def makeHamiltonian (eigenenergies, hopping, size=None, alignment = 'top left'):
    alignment = alignment.lower()
    if size==None:
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

def fermifuncLeft (sites, maxVal, decay, offset):
    values = [maxVal/(1+np.exp(decay*(i+1-offset))) for i in sites]
    return np.array(values)

def fermifuncRight (sites, maxVal, decay, offset, sizeTotal):
    values = [maxVal/(1+np.exp(decay*(sizeTotal-i-offset))) for i in sites]
    return np.array(values)

def combineH (sample, leadLeft, leadRight, hoppingInter):
    sizeMax = max(len(sample), len(leadLeft), len(leadRight))
    
    arrayNew = np.zeros((sizeMax, sizeMax), dtype=complex)
    
    length = int((lengthTotal-lengthSample)/2)
    for row in np.arange(length-1, sizeMax-length):
        for column in np.arange(length-1, sizeMax-length):
            if row-length >= 0 and column-length >= 0:
                arrayNew[row, column] = sample[row-length, column-length]
    
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
            
            arrayNew[row, column] += newFactor + leadLeft[row][column] + leadRight[row][column]
    return arrayNew

def GreensFunc (omega, totalSystem):
    array = np.identity(len(totalSystem))*omega - totalSystem
    return np.linalg.inv(array)

def semiInfiniteLeads (omega, Hamiltonian, hoppingInter, hoppingBath, size, eta=0.000000001):
    sigmaL = np.zeros((size, size), dtype=complex)
    sigmaR = np.zeros((size, size), dtype=complex)
    
    #sigmaL(1,1) = t_c^2 *(w +eta*1i)/(2*v^2)*(1-sqrt(1-(4*v^2)/(w+eta *1i)^2));
    #sigmaR(N,N) =t_c^2 *(w +eta*1i)/(2*v^2)*(1-sqrt(1-(4*v^2)/(w+eta *1i)^2));
    
    #G = (w +eta*1i)/(2*v^2)*(1-sqrt(1-(4*v^2)/(w+eta *1i)^2))
    G = (omega+1j*eta)/(2*hoppingBath**2)*(1-np.sqrt(1-(4*hoppingBath**2)/(omega+1j*eta)**2))
    
    sigmaL[0,0] = hoppingInter**2 * G
    sigmaR[size-1,size-1] =hoppingInter**2 * G
    
    sigmaLdagger = np.conj(sigmaL.T)
    sigmaRdagger = np.conj(sigmaR.T)
    
    gammaR = 1j*(sigmaR - sigmaRdagger);
    gammaL = 1j* (sigmaL - sigmaLdagger);
    
    Diag = (omega+1j*eta)*np.identity(size, dtype=complex)
    GreensInv = Diag - (Hamiltonian + sigmaL + sigmaR)
    Greens = np.linalg.inv(GreensInv);
    
    return Greens, gammaL, gammaR

def Transmission (GreensFunc, couplingLeft, couplingRight, value):
    # calculate the matrix product
    transport = np.dot(np.dot(GreensFunc, (couplingLeft - couplingRight)), np.conj(GreensFunc.T))
    
    # calculate the current
    Transmission = np.imag(transport[value, value+1])
    
    return Transmission

#%% Variables
# variables for the sample
lengthSample = 96
eigenenergy = 0.1
hopping = 1

# variables for the leads
lengthTotal = 256
maxVal = 1
decay = 0.2 #0.3
offset = 32
hoppingLead = hopping

# variables for the final calculation
step = 0.005
omegas = np.arange(-2, 2+step, step)

#%% preperatory computation

#%%% prepare the Hamiltonians and coupling strengths of the leads
# compute the fermi functions
sites = range(lengthTotal)
fermifuncLeft = fermifuncLeft(sites, maxVal, decay, offset)
fermifuncRight = fermifuncRight(sites, maxVal, decay, offset, lengthTotal)

Plot = False
if Plot == True:
    FigFermi, AxFermi = plt.subplots(1, 1)
    AxFermi.plot(sites, fermifuncLeft)
    AxFermi.plot(sites, fermifuncRight)
    AxFermi.plot(sites, np.flip(fermifuncRight))

# generate the Hamiltonians of the leads
length = int((lengthTotal-lengthSample)/2)
leadLeft = makeHamiltonian(1j*fermifuncLeft, hoppingLead, size=length, alignment = 'top left')
leadRight = makeHamiltonian(1j*fermifuncRight, hoppingLead, size=length, alignment = 'bottom right')

# compute the coupling strengths
couplingLeft = 1j*(leadLeft-leadLeft.conj().T)
couplingRight = 1j*(leadRight-leadRight.conj().T)

averageTimes = 1#20
avTransmissionEM = list()
avTransmissionEM.append(list())
avTransmissionSI = list()
avTransmissionSI.append(list())
for i in range(averageTimes):
    #%%% compute the Hamiltonian of the Sample
    sample = makeHamiltonian(randomNum(eigenenergy, lengthSample), hopping)

    # combine the Hamiltonians of the sample and the leads
    hoppingInter = hopping
    totalSystem = combineH(sample, leadLeft, leadRight, hoppingInter)

    #%% compute the Transmission
    TransmissionsEM = list()
    TransmissionsSI = list()
    for omega in omegas:
        value = int((lengthTotal-lengthSample)/2)
        
        GreensFunctionEM = GreensFunc(omega, totalSystem)
        T_EM = Transmission(GreensFunctionEM, couplingLeft, couplingRight, value)
        
        GreensFunctionSI, couplingLeftSI, couplingRightSI = semiInfiniteLeads(omega, sample, hoppingInter, hoppingLead, lengthSample)
        T_SI = Transmission(GreensFunctionSI, couplingLeftSI, couplingRightSI, 0)
        
        print(i+1, 'done', omega, T_EM, T_SI)
        TransmissionsEM.append(T_EM)
        TransmissionsSI.append(T_SI)
    
    avTransmissionEM.append(TransmissionsEM)
    avTransmissionSI.append(TransmissionsSI)

for i in range(max(len(avTransmissionEM[1]), len(avTransmissionSI[1]))):
    valuesEM = [avTransmissionEM[j][i] for j in np.arange(1, len(avTransmissionEM))]
    averageEM = sum(valuesEM)/len(valuesEM)
    avTransmissionEM[0].append(averageEM)
    
    valuesSI = [avTransmissionSI[j][i] for j in np.arange(1, len(avTransmissionSI))]
    averageSI = sum(valuesSI)/len(valuesSI)
    avTransmissionSI[0].append(averageSI)

#%% plot the final Graph
def size(factor = 1, offset = 0):
    standardSize = 14.53
    size = factor*standardSize + offset
    return size/2.54 # cm -> in

def frac(standardSize, factor = 1, offset = 0):
    size = factor*standardSize + offset/2.54
    return size/standardSize

#%%% prepare the plot
height = size(factor=1/2) #cm
width = size(factor=1) #cm

plt.rcParams["figure.figsize"] = (height,width)
plt.rcParams["font.family"] = "Arial"
plt.rcParams.update({'font.size': 20})

Fig, Ax = plt.subplots(1, 1)
Fig.set_figwidth(width)
Fig.set_figheight(height)

plt.rcParams['text.usetex'] = False

layout = False
#plt.tight_layout()
if layout == True:
    leftFFT = frac(width, factor=0, offset=1.1)
    bottomFFT = frac(height, factor=0, offset=1.25)
    rightFFT = frac(width, offset=-0.1)
    topFFT = frac(height, offset=-0.1)
    plt.subplots_adjust(left = leftFFT, 
                        bottom = bottomFFT, 
                        right = rightFFT, 
                        top = topFFT, 
                        wspace=None, hspace=None)
    plt.minorticks_on()
    
    widthAx = (rightFFT - leftFFT)*width
    heightAx = (topFFT - bottomFFT)*height
    
#%%%% plot the data
Ax.plot(omegas, avTransmissionEM[0])
Ax.plot(omegas, avTransmissionSI[0])
#Ax.set_yscale("log")

#%%%% change Axes
Title = 'Eigenenergy: w='+str(eigenenergy)+'; hopping: t='+str(hopping)
Ax.set_title(Title)
Ax.set_xlabel('omega')
Ax.set_ylabel('Transmission')

#%%%% show the plot
plt.show()
