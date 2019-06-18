#!/usr/bin/python
import os
import sys
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
plt.ion()
import pickle
sys.path.append('/home/ntv/integrate/analysis/')
import pySlice
from tqdm import tqdm
from skimage import measure
from scipy.ndimage import convolve, rotate
from timeit import default_timer as timer
import sys
import pandas as pd
from tqdm import tqdm
import tensorflow as tf
from keras.backend.tensorflow_backend import set_session
config = tf.ConfigProto(device_count = {'GPU':0})
#config.gpu_options.per_process_gpu_memory_fraction = 0.15

set_session(tf.Session(config=config))
from keras.models import Model, load_model, save_model
import mltools
from mantid.geometry import SpaceGroupFactory, PointGroupFactory
import itertools

popList = []
for i in range(len(sys.path))[::-1]:
    if 'antid' in sys.path[i]:
        sys.path.pop(i)
import socket
if 'sns' in socket.gethostname():
    sys.path.append('/SNS/users/ntv/mantid/mantid/release/bin')
else: 
    sys.path.append('/home/ntv/workspace/mantid/release/bin/')

from mantid.simpleapi import *
import ICCFitTools as ICCFT


np.random.seed = 42
nX = 32
nY = 32
nZ = 32
nChannels = 1

#Do some initial stuff for tensorflow
#model_file = 'model_keras.h5' #First pass
#model_file = '/home/ntv/ml_peak_integration/models/model_withQMask_fromdfpeaks_relu_halfrot_strongonly_0p5dropout.h5'
#trainedOnHKL = False; thresh=0.15
#model_file = '/home/ntv/ml_peak_integration/models/model_withQMask_fromdfpeaks_relu_halfrot_allruns_limitNoise_hklFull.h5'
#trainedOnHKL = True; thresh=0.15
#model_file = '/home/ntv/ml_peak_integration/models/beta_lac_secondxtal_0p15peakThreshold_2.h5' #good one
model_file = '/home/ntv/ml_peak_integration/models/folded_second_betalac_xtal.h5' #combined
trainedOnHKL = False; thresh=0.15
model = load_model(model_file, custom_objects={'bce_dice_loss':mltools.bce_dice_loss, 'dice_coeff':mltools.dice_coeff, 
                                               'dice_loss':mltools.dice_loss, 'mean_iou':mltools.mean_iou})

#Load our mantid data

#beta lac july 2018 second xtal
peaksFile = '/SNS/users/ntv/dropbox/beta_lac_july2018_secondxtal_mbvg_2/beta_lac_secondcrystal_combined_pf.integrate'
UBFile = '/home/ntv/mandi_preprocessing/beta_lac_july2018/beta_lac_july2018_secondxtal.mat'
DetCalFile = '/home/ntv/mandi_preprocessing/MANDI_June2018.DetCal'
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
nxsTemplate = '/SNS/MANDI/IPTS-8776/nexus/MANDI_%i.nxs.h5'
dQPixel=0.0035
q_frame = 'lab'
pplmin_frac=0.9; pplmax_frac=1.1; mindtBinWidth=15; maxdtBinWidth=50
moderatorFile = '/home/ntv/integrate/bl11_moderatorCoefficients_2018.dat'

# Some parameters
importPeaks = True
print('Loading peaks_ws')
for ws in mtd.getObjectNames():
    if mtd[ws].getComment() == '%s'%peaksFile:
        print '    using already loaded peaks file'
        importPeaks = False
        peaks_ws = mtd[ws]
if importPeaks:
    peaks_ws = LoadIsawPeaks(Filename = peaksFile)
    peaks_ws.setComment(peaksFile)

LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()
dQ = np.abs(ICCFT.getDQFracHKL(UBMatrix, frac=0.5))

df = pd.DataFrame(peaks_ws.toDict())
df = df[df['DSpacing']>1.6]

def getReflectionFamily(row):
    sg = SpaceGroupFactory.createSpaceGroup('P 32 2 1')
    pg = sg.getPointGroup()
    h = float(row['h'])
    k = float(row['k'])
    l = float(row['l'])
    hklPlus = np.array([h,k,l])
    hklMinus = np.array([-h, -k, -l])
    refFamPlus = pg.getReflectionFamily(hklPlus)
    refFamMinus = pg.getReflectionFamily(hklMinus)
    checkIDX = (~np.equal(refFamPlus, refFamMinus)).argmax()
    if refFamPlus[checkIDX] >= refFamMinus[checkIDX]:
        return tuple(np.array(refFamPlus,dtype=int))
    else: # refFamPlus[checkIDX] < refFamMinus[checkIDX]:
        return tuple(np.array(refFamMinus,dtype=int))

df['reflectionFamily'] = df.apply(getReflectionFamily, axis=1)
runNumbersAll = range(9113,9117+1)

runNumberHashes = {}
for i in range(1,len(runNumbersAll)+1):
    combos = itertools.combinations(runNumbersAll, i)
    for combo in combos:
        strKey = str(np.array(combo))
        runNumberHashes[strKey] = hash(strKey)

dfRefFam = pd.DataFrame(df['reflectionFamily'].value_counts())
dfRefFam['meanI'] = df.groupby('reflectionFamily')['Intens'].mean()
dfRefFam['meanSig'] = df.groupby('reflectionFamily')['SigInt'].mean()
dfRefFam['RunNumbers'] = df.groupby('reflectionFamily')['RunNumber'].apply(np.unique)
dfRefFam['RunNumbersHashed'] = dfRefFam['RunNumbers'].apply(lambda x: hash(str(x)))
dfRefFam['NumRuns'] = dfRefFam['RunNumbers'].apply(len)
dfRefFam['PeakNumbers'] = df.groupby('reflectionFamily')['PeakNumber'].apply(np.unique)
dfRefFam['DSpacing'] = df.groupby('reflectionFamily')['DSpacing'].apply(np.unique)
df['hkl'] = df[['h', 'k', 'l']].apply(tuple,axis=1)
dfRefFam['hkl'] = df.groupby('reflectionFamily')['hkl'].apply(list)

dfRefFam = dfRefFam.rename(index=str, columns={"reflectionFamily":"counts"})
dfRefFam['IntensML'] = np.zeros(len(dfRefFam))
dfRefFam['SigIntML'] = np.ones(len(dfRefFam),dtype=float)
dfRefFam['meanBG'] = np.zeros(len(dfRefFam))
dfRefFam['numVoxelsInPeak'] = np.zeros(len(dfRefFam))


if len(sys.argv) == 1:
    runNumbers = df['RunNumber'].unique()
else:
    runNumbers = map(int,sys.argv[1:])

currentHash = hash(str(np.array(runNumbers)))

#Need to load all of the MD workspaces into memory
for rn in runNumbers:
    print 'Loading MD for {}'.format(rn)
    LoadMD(Filename='/home/ntv/Desktop/data_folding/mde/MANDI_mde_sample_{}.nxs'.format(rn), OutputWorkspace='md_{}'.format(rn))
    LoadIsawUB(InputWorkspace='md_{}'.format(rn), Filename=UBFile)
    
#Solid angle and spectrum workspaces
flux_cdf = Load(Filename='/home/ntv/Desktop/data_folding/spectra.nxs', OutputWorkspace='flux_cdf')
sa = Load(Filename='/home/ntv/Desktop/data_folding/MANDI_solidAngle_nomask.nxs', OutputWorkspace='sa')    

print('Integrating run numbers:', runNumbers)
qMask = pickle.load(open('/data/ml_peak_sets/beta_lac_secondcrystal_0p4qMask_0p15peakThreshold_folded/qMask.pkl', 'rb'))

cX, cY, cZ = np.array(qMask.shape)//2
dX, dY, dZ = nX//2, nY//2, nZ//2
qMaskSimulated = qMask[cX-dX:cX+dX, cY-dY:cY+dY, cZ-dZ:cZ+dZ]

numBad = 0
numFit = 0
paramList = []
integrateIDX = dfRefFam['RunNumbersHashed'] == currentHash
dfRefFam = dfRefFam[integrateIDX]
hklsToIntegrate = np.array(dfRefFam.index)
tmp = range(len(hklsToIntegrate))
np.random.shuffle(tmp)
#hklsToIntegrate = hklsToIntegrate[tmp]
for hklToMerge in tqdm(hklsToIntegrate):
    try:
        #---Determine which peaks to try and then run them.
        dfRefKey = hklToMerge
        hklToMerge = tuple(np.array(hklToMerge.strip('()').split(','), dtype=float)) #Convert back to tuple
        peaksToSum = df[(df['reflectionFamily'] == hklToMerge)].index

        #Clear the accumulation workspaces for the peaks
        if mtd.doesExist('dataMD'):
            DeleteWorkspace('dataMD')
        if mtd.doesExist('normMD'):
            DeleteWorkspace('normMD')
        if mtd.doesExist('result'):
            DeleteWorkspace('result')    
        #Add all of our peaks    
        #d = np.zeros([39,59,37])
        #n = np.zeros([39,59,37])
        #errorsq  = np.zeros([39,59,37])

        for peakInSeries, peakNumber in enumerate(peaksToSum):
            with open('/home/ntv/Desktop/data_folding/progress_reports/%s.txt'%runNumbers, 'w') as f:
                f.write(str(peakNumber))
            peak = peaks_ws.getPeak(peakNumber);
            print('==================================',peakNumber)
            if mtd.doesExist('result'):
                DeleteWorkspace('result')   
            qx, qy, qz = peak.getQSampleFrame()
            runNumber = peak.getRunNumber()
            MDdata = mtd['md_{}'.format(runNumber)]
            MDNorm(InputWorkspace=MDdata,
                   SolidAngleWorkspace='sa',
                   RLU=False,
                   FluxWorkspace='flux_cdf',
                   QDimension0='1,0,0',
                   QDimension1='0,1,0',
                   QDimension2='0,0,1',
                   Dimension0Name='QDimension0',
                   Dimension0Binning='{0},{1},{2}'.format(qx-dQ[0][0], dQPixel, qx+dQ[0][1]),
                   Dimension1Name='QDimension1',
                   Dimension1Binning='{0},{1},{2}'.format(qy-dQ[1][0], dQPixel, qy+dQ[1][1]),
                   Dimension2Name='QDimension2',
                   Dimension2Binning='{0},{1},{2}'.format(qz-dQ[2][1], dQPixel, qz+dQ[2][1]),
                   TemporaryDataWorkspace=None, 
                   TemporaryNormalizationWorkspace=None, 
                   OutputWorkspace='result',
                   OutputDataWorkspace='dataMD',
                   OutputNormalizationWorkspace='normMD')   
            if peakInSeries == 0:
                data = mtd['dataMD'].getSignalArray().copy()
                errorsq = mtd['dataMD'].getErrorSquaredArray().copy()
                norm = mtd['normMD'].getSignalArray().copy()
            else:
                data += mtd['dataMD'].getSignalArray()
                errorsq += mtd['dataMD'].getErrorSquaredArray()
                norm += mtd['normMD'].getSignalArray()
            print('Finished MDNorm. Peak number %i'%peakNumber)

        normData = 1.*data/norm
        edgeConvBox = np.ones([3,3,3])
        newMask = convolve(~np.isfinite(norm), edgeConvBox).copy()
        normData[newMask] = 0.

        mtd['dataMD'].setSignalArray(data)
        mtd['dataMD'].setErrorSquaredArray(errorsq)
        mtd['result'].setSignalArray(normData)
        normErrorSq = 1.*errorsq/norm/norm
        normErrorSq[newMask] = 0.
        mtd['result'].setErrorSquaredArray(normErrorSq)

        box = mtd['dataMD']
        n_events_cropped, n_errorsq_cropped, image = mltools.getImageFromBox(box, UBMatrix, peak, rebinToHKL=trainedOnHKL, qMaskSimulated=qMaskSimulated, returnErrorSq=True)
        peakMask, testim, blobs = mltools.getPeakMask(image, model,thresh=thresh)
        box = mtd['result']

        n_events_cropped, n_errorsq_cropped, image2 = mltools.getImageFromBox(box, UBMatrix, peak, rebinToHKL=trainedOnHKL, qMaskSimulated=qMaskSimulated, returnErrorSq=True)
        peakIDX = np.logical_and(peakMask, np.isfinite(n_events_cropped))
        n_events_cropped = np.ma.masked_invalid(n_events_cropped)

        #Integration
        intensScalFactor = 1000.0
        nPeak = peakIDX.sum()
        bgIDX = np.logical_and(~peakIDX, np.isfinite(n_events_cropped))

        bgConvBox = np.ones([7,7,7])
        bgConvBox2 = np.ones([5,5,5])
        bgIDX = convolve(peakIDX, bgConvBox)
        notBGIDX = convolve(peakIDX, bgConvBox2)
        bgIDX = reduce(np.logical_and, [bgIDX, ~notBGIDX, ~peakIDX, np.isfinite(n_events_cropped), qMaskSimulated])

        nBG = bgIDX.sum()
        countsInPeak = np.sum(n_events_cropped[peakIDX])
        countsInBG = np.sum(n_events_cropped[bgIDX])
        intensity = countsInPeak - 1.0*countsInBG/nBG*nPeak
        sigSqr = np.sum(n_errorsq_cropped[peakIDX]) + np.mean(n_errorsq_cropped[bgIDX])*nPeak
        sigma = np.sqrt(sigSqr)
        intensity = intensity*intensScalFactor
        sigma = sigma*intensScalFactor
       
        #Record the results
        peakDict = {}
        peakDict['PeakNumber'] = peaksToSum
        peakDict['Intens'] = [peaks_ws.getPeak(i).getIntensity() for i in peaksToSum]
        peakDict['SigInt'] = [peaks_ws.getPeak(i).getSigmaIntensity() for i in peaksToSum]
        peakDict['theta'] = [peaks_ws.getPeak(i).getScattering()*0.5 for i in peaksToSum]
        peakDict['numVoxelsSimulated'] = peakIDX.sum()
        paramList.append(peakDict)
        
        dfRefFam.at[dfRefKey,'IntensML'] = intensity
        dfRefFam.at[dfRefKey,'SigIntML'] = sigma
        dfRefFam.at[dfRefKey,'meanBG'] = np.mean(intensScalFactor*n_events_cropped[bgIDX])
        dfRefFam.at[dfRefKey,'numVoxelsInPeak'] = np.sum(peakMask)
        try: 
            print('old: %4.2f +- %4.2f (%4.2f)'%(np.mean(peakDict['Intens']),np.mean(peakDict['SigInt']), np.mean(peakDict['Intens'])/np.mean(peakDict['SigInt'])))
        except:
            print('old: %4.2f +- %4.2f (%4.2f)'%(np.mean(peakDict['Intens']),np.mean(peakDict['SigInt']), np.inf))
        try: 
            print('new (%i voxels): %4.2f +- %4.2f (%4.2f)'%(nPeak, intensScalFactor*intensity, intensScalFactor*sigma, intensity/sigma))
        except:
            print('new (%i voxels): %4.2f +- %4.2f (%4.2f)'%(nPeak, intensScalFactor*intensity, intensScalFactor*sigma, np.inf))
        print('****Finished fitting hkl set', hklToMerge)
        
        numFit += 1
        if numFit%50 == 0:
            print('Saving!')
            paramFileName = '/home/ntv/Desktop/data_folding/ml_results/beta_lac_secondcrystal_combined_folded_paramList_%s.pkl'%runNumbers
            paramFileName = paramFileName.replace(' ', '_').replace('[','').replace(']','')
            dfFileName = '/home/ntv/Desktop/data_folding/ml_results/beta_lac_secondcrystal_combined_folded_%s.pkl'%runNumbers
            dfFileName = dfFileName.replace(' ', '_').replace('[','').replace(']','')
            pickle.dump(paramList, open(paramFileName, 'wb'))
            dfRefFam.to_pickle(dfFileName)            
    except KeyboardInterrupt:
        0/0
    except:
        print('Error fitting hkl set', hklToMerge)
        numBad += 1        
print('Saving!')
paramFileName = '/home/ntv/Desktop/data_folding/ml_results/beta_lac_secondcrystal_combined_folded_paramList_%s.pkl'%runNumbers
paramFileName = paramFileName.replace(' ', '_').replace('[','').replace(']','')
dfFileName = '/home/ntv/Desktop/data_folding/ml_results/beta_lac_secondcrystal_combined_folded_%s.pkl'%runNumbers
dfFileName = dfFileName.replace(' ', '_').replace('[','').replace(']','')
pickle.dump(paramList, open(paramFileName, 'wb'))
dfRefFam.to_pickle(dfFileName)           

    




