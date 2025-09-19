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

def Hamiltonian (eigenenergies, hopping, size=None, alignment = 'top left'):
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

def combineH (sample, leadLeft, leadRight, hopping):
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
                newFactor = hopping
            elif row == values[0] and column == values[0]-1:
                newFactor = hopping
            if row == values[1]-1 and column == values[1]:
                newFactor = hopping
            if row == values[1] and column == values[1]-1:
                newFactor = hopping
            
            arrayNew[row, column] += newFactor + leadLeft[row][column] + leadRight[row][column]
    return arrayNew

def Transmission (omega, lengthTotal, lengthSample, GreensFunc, couplingLeft, couplingRight):
    # calculate the matrix product
    transport = np.dot(np.dot(GreensFunc, (couplingLeft - couplingRight)), np.conj(GreensFunc.T))
    
    
    # calculate the current
    value = int((lengthTotal-lengthSample)/2)
    Transmission = np.imag(transport[value, value+1])
    
    return Transmission

def GreensFunc (omega, totalSystem):
    array = np.identity(len(totalSystem))*omega - totalSystem
    return np.linalg.inv(array)

#%% Variables
# variables for the sample
lengthSample = 96
eigenenergy = 0
hopping = 1

# variables for the leads
lengthTotal = 256
maxVal = 1
decay = 0.3 #0.3
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
leadLeft = Hamiltonian(1j*fermifuncLeft, hoppingLead, size=length, alignment = 'top left')
leadRight = Hamiltonian(1j*fermifuncRight, hoppingLead, size=length, alignment = 'bottom right')

# compute the coupling strengths
couplingLeft = 1j*(leadLeft-leadLeft.conj().T)
couplingRight = 1j*(leadRight-leadRight.conj().T)

averageTimes = 1#20
avTransmission = list()
avTransmission.append(list())
for i in range(averageTimes):
    #%%% compute the Hamiltonian of the Sample
    sample = Hamiltonian(randomNum(eigenenergy, lengthSample), hopping)

    # combine the Hamiltonians of the sample and the leads
    hoppingInter = hopping
    totalSystem = combineH(sample, leadLeft, leadRight, hoppingInter)

    #%% compute the Transmission
    Transmissions = list()
    for omega in omegas:
        GreensFunction = GreensFunc(omega, totalSystem)
        T = Transmission(omega, lengthTotal, lengthSample, GreensFunction, couplingLeft, couplingRight) 
        print(i+1, 'done', omega, T)
        Transmissions.append(T)
    
    avTransmission.append(Transmissions)

for i in range(len(avTransmission[1])):
    values = [avTransmission[j][i] for j in np.arange(1, len(avTransmission))]
    average = sum(values)/len(values)
    avTransmission[0].append(average)

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
Ax.plot(omegas, Transmissions)
#Ax.set_yscale("log")

#%%%% change Axes
Title = 'Eigenenergy: w='+str(eigenenergy)+'; hopping: t='+str(hopping)
Ax.set_title(Title)
Ax.set_xlabel('omega')
Ax.set_ylabel('Transmission')

#%%%% show the plot
plt.show()
