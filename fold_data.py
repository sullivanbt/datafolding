import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from mantid.simpleapi import *
from mantid.geometry import SpaceGroupFactory, PointGroupFactory
import ICCFitTools as ICCFT


if False:
    #LoadNexus(OutputWorkspace='flux', Filename='/home/ntv/Desktop/data_folding/spectra_wavelength.nxs')
    flux_cdf = Load(Filename='/home/ntv/Desktop/data_folding/spectra.nxs', OutputWorkspace='flux_cdf') # 2019-02-20T20:25:48.718948000
    sa = Load(Filename='/home/ntv/Desktop/data_folding/MANDI_solidAngle.nxs', OutputWorkspace='sa') # 2019-02-20T20:26:07.526985000
    runNumbers = [9113, 9114, 9115, 9116, 9117]
    nxsFormat = '/data/betalac_secondxtal/MANDI_%i.nxs.h5'

    #Set up the space group
    sg = SpaceGroupFactory.createSpaceGroup('P 32 2 1')
    pg = sg.getPointGroup()
    peaks_ws = LoadIsawPeaks('/SNS/users/ntv/dropbox/beta_lac_july2018_secondxtal_mbvg_2/beta_lac_secondcrystal_combined_pf.integrate')
    LoadIsawUB(InputWorkspace=peaks_ws,Filename='/home/ntv/mandi_preprocessing/beta_lac_july2018/beta_lac_july2018_secondxtal.mat')

    #Set up the dataframes for peaks so we can merge equivalent peaks
    df = pd.DataFrame(peaks_ws.toDict())
    def getReflectionFamily(row):
        sg = SpaceGroupFactory.createSpaceGroup('P 32 2 1')
        pg = sg.getPointGroup()
        return tuple(pg.getReflectionFamily(np.array([row['h'], row['k'], row['l']],dtype=int)))
    df['reflectionFamily'] = df.apply(getReflectionFamily, axis=1)

    def getHKL(row):
        return tuple(np.array([row['h'], row['k'], row['l']))

    #A dataframe for the reflection families so we know how bright peaks are on average
    dfRefFam = pd.DataFrame(df['reflectionFamily'].value_counts())
    dfRefFam['meanI'] = df.groupby('reflectionFamily')['Intens'].mean()
    dfRefFam['meanSig'] = df.groupby('reflectionFamily')['SigInt'].mean()
    dfRefFam['hkl'] = df.groupby('reflectionFamily')
    dfRefFam = dfRefFam.rename(index=str, columns={"reflectionFamily":"counts"})

    for rn in runNumbers:
        print 'Loading MD for {}'.format(rn)
        LoadMD(Filename='/home/ntv/Desktop/data_folding/mde/MANDI_mde_sample_{}.nxs'.format(rn), OutputWorkspace='md_{}'.format(rn))
        LoadIsawUB(InputWorkspace='md_{}'.format(rn),Filename='/home/ntv/mandi_preprocessing/beta_lac_july2018/beta_lac_july2018_secondxtal.mat')

def cleanMDBoxes():
    for wsName in mtd.getObjectNames():
        if 'MDbox_' in wsName:
            mtd.remove(wsName)



q_frame='sample'
cleanMDBoxes()

#hklToMerge = (22,-4,-1)
#hklToMerge = (2,1,26) #Great for demonstrations!
hklToMerge = (1,1,29)
peaksToSum = df[df['reflectionFamily'] == hklToMerge].index
meanI = df.loc[peaksToSum]['Intens'].mean()
meanSigI = df.loc[peaksToSum]['SigInt'].mean()



#Clear the output workspaces for the peaks
if mtd.doesExist('dataMD'):
    DeleteWorkspace('dataMD')
if mtd.doesExist('normMD'):
    DeleteWorkspace('normMD')
for peakNumber in peaksToSum:
    peak = peaks_ws.getPeak(peakNumber)
    qx, qy, qz = peak.getQSampleFrame()
    dqPixel = 0.003
    nVoxels = 30
    dQ = dqPixel*nVoxels/2
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
           Dimension0Binning='{0},{1},{2}'.format(qx-dQ, dqPixel, qx+dQ),
           Dimension1Name='QDimension1',
           Dimension1Binning='{0},{1},{2}'.format(qy-dQ, dqPixel, qy+dQ),
           Dimension2Name='QDimension2',
           Dimension2Binning='{0},{1},{2}'.format(qz-dQ, dqPixel, qz+dQ),
           TemporaryDataWorkspace='dataMD' if mtd.doesExist('dataMD') else None,
           TemporaryNormalizationWordspace='normMD' if mtd.doesExist('normMD') else None,
           OutputWorkspace='result',
           OutputDataWorkspace='dataMD',
           OutputNormalizationWorkspace='normMD')           
    Box = BinMD(InputWorkspace=MDdata,
                AlignedDim0='Q_%s_x,' % q_frame +
                str(qx-dQ)+','+str(qx+dQ)+','+str(nVoxels),
                AlignedDim1='Q_%s_y,' % q_frame +
                str(qy-dQ)+','+str(qy+dQ)+','+str(nVoxels),
                AlignedDim2='Q_%s_z,' % q_frame +
                str(qz-dQ)+','+str(qz+dQ)+','+str(nVoxels),
                OutputWorkspace='MDbox_{}'.format(peakNumber))

rA = 0.023
rB = 0.025
rC = 0.027
nScaled = mtd['result'].getSignalArray()
errorSq = mtd['result'].getErrorSquaredArray()

QX, QY, QZ = ICCFT.getQXQYQZ(Box)
rSqr = (QX-qx)*(QX-qx) + (QY-qy)*(QY-qy) + (QZ-qz)*(QZ-qz)
peakIDX = np.logical_and(rSqr <= rA*rA, ~np.isnan(nScaled))
bgIDX = np.logical_and(rSqr < rC*rC, rSqr > rB*rB)
bgIDX = np.logical_and(bgIDX, ~np.isnan(nScaled))
nPeak = peakIDX.sum()
nBG = bgIDX.sum()

countsInPeak = np.sum(nScaled[peakIDX])
countsInBG = np.sum(nScaled[bgIDX])
intensity = countsInPeak - 1.0*countsInBG/nBG*nPeak
sigSqr = np.sum(errorSq[peakIDX]) + np.sum(errorSq[bgIDX])
sigma = np.sqrt(sigSqr)
print('Original {:d} peaks with meanI of {:4.2f} +- {:4.2f} ({:4.2f})'.format(len(peaksToSum), meanI, meanSigI, meanI/meanSigI))
print('Intensity: {:4.2f} +- {:4.2f} ({:4.2f})'.format(1.e5*intensity, 1.e5*sigma, intensity/sigma))

plt.figure(1)
plt.clf()
numRows = (len(peaksToSum)+1) // 3
if (len(peaksToSum)+1) % 3 > 0:
    numRows += 1
for i, peakNumber in enumerate(peaksToSum[:11]):
    peak = peaks_ws.getPeak(peakNumber)
    plt.subplot(numRows,3,i+1)
    plt.imshow(mtd['MDbox_{}'.format(peakNumber)].getNumEventsArray()[:,:,15])
    plt.axis('off')
    plt.title('Run {:d} ({:d},{:d},{:d})'.format(int(peak.getRunNumber()), int(peak.getH()), int(peak.getK()), int(peak.getL())))
    
plt.subplot(numRows,3,numRows*3)
plt.imshow(mtd['result'].getSignalArray()[:,:,15])
plt.axis('off')
plt.title('Normalized Sum')
plt.suptitle('Folded peaks for reflection family ({:d},{:d},{:d})'.format(int(hklToMerge[0]), int(hklToMerge[1]), int(hklToMerge[2])))



