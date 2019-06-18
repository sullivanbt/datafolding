from tqdm import tqdm
from mantid.simpleapi import *
import numpy as np
import matplotlib.pyplot as plt

#Load
LoadNexus(OutputWorkspace='flux', Filename='/home/ntv/Desktop/data_folding/spectra_wavelength.nxs')
#Load(Filename='/home/ntv/Desktop/data_folding/MANDI_solidAngle.nxs', OutputWorkspace='sa') # 2019-02-20T20:26:07.526985000
LoadIsawPeaks(Filename='/home/ntv/Desktop/beta_lac_july2018_firstxtal/peaks_profileFitted_ws_9078_mandi_autoreduced_jan19.integrate',
              OutputWorkspace='peaks_ws')

flux = mtd['flux']
peaks_ws = mtd['peaks_ws']
numDetectors = flux.getNumberHistograms()
numPeaks = peaks_ws.getNumberPeaks()
spec_wl = flux.readX(0)


#Map bank names to spec numbers
bankDict = {}
inst = flux.getInstrument()
for i in range(numDetectors):
    trial_detid = flux.getSpectrum(i).getDetectorIDs()[500]
    det = inst.getDetector(trial_detid)
    detName = det.getName().split('(')[0]
    bankDict[detName] = i

peaks_ws_clone = CloneWorkspace(InputWorkspace='peaks_ws', OutputWorkspace='peaks_ws_clone')
normFactorList = []
for peakNumber in tqdm(range(numPeaks)):
    try:
        peak = peaks_ws_clone.getPeak(peakNumber)
        peakDetID = peak.getDetectorID()
        peak_wl = peak.getWavelength()
        peak_intens = peak.getIntensity()
        peak_sigint = peak.getSigmaIntensity()

        specNumber = bankDict[peak.getBankName()]
        spec_intens = flux.readY(specNumber)
        spec_wl_idx = np.argmin(np.abs(spec_wl-peak_wl))
        normFactor = spec_intens[spec_wl_idx] / 194793.
        
        if normFactor > 0:
            peak.setIntensity(peak_intens / normFactor)
            peak.setSigmaIntensity(peak_sigint / normFactor)
            normFactorList.append(normFactor)
    except KeyboardInterrupt:
        0/0
    except:
        print 'Error with peak {}'.format(peakNumber)
        pass

SaveLauenorm(InputWorkspace=peaks_ws_clone, Filename='/home/ntv/Desktop/laue_lorentz/using_spectrum/laueNorm', SortFilesBy='RunNumber', MinIsigI=1.0, MinWavelength=2.0, MaxWavelength=4.0, ScalePeaks=3.0, WidthBorder=3)   





