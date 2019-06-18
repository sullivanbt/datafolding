from mantid.simpleapi import *
from tqdm import tqdm
event_ws = Load(Filename='/SNS/MANDI/IPTS-20460/nexus/MANDI_9300.nxs.h5')

# Create an xml that maps spectrum number to banks

# 1) map detector ids to banks
bankDict = {}
specDict = {}
for i in range(100):
    inst = event_ws.getInstrument()
    bank = inst.getComponentByName("bank%i"%i)
    if bank is not None:
        lowID = bank.minDetectorID()
        highID = bank.maxDetectorID()
        bankDict['bank%i'%i] = (lowID, highID)
        specDict['bank%i'%i] = [1e10, -1]
# Now spectrum ids to detector ids, which are mapped to banks

numSpectra = event_ws.getNumberHistograms()
for i in tqdm(range(numSpectra)):
    spec = event_ws.getSpectrum(i)
    detID = spec.getDetectorIDs()[0]
    for key in bankDict.keys():
        lowID = bankDict[key][0]
        highID = bankDict[key][1]
        if lowID <= detID and detID <= highID:
            if i < specDict[key][0]:
                specDict[key][0] = i
            elif i > specDict[key][1]:
                specDict[key][1] = i
                
# Write the file
with open('/home/ntv/Desktop/data_folding/mandi_grouping.xml', 'w') as f:
    f.write("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n")
    f.write("<detector-grouping>\n")
    for i in range(100):
        if 'bank%i'%i in specDict.keys():
            lowID = specDict['bank%i'%i][0]
            highID = specDict['bank%i'%i][1]
            f.write("<group name=\"%s\"> <ids val=\"%i-%i\"/> </group>\n"%("bank%i"%i, lowID, highID))
    f.write("</detector-grouping>")

""" # Creates an xml that maps detector ids to banks
with open('/home/ntv/Desktop/data_folding/mandi_grouping.xml', 'w') as f:
    f.write("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n")
    f.write("<detector-grouping>\n")
    for i in range(100):
        inst = event_ws.getInstrument()
        bank = inst.getComponentByName("bank%i"%i)
        if bank is not None:
            lowID = bank.minDetectorID()
            highID = bank.maxDetectorID()
            f.write("<group name=\"%s\"> <ids val=\"%i-%i\"/> </group>\n"%("bank%i"%i, lowID, highID))
    f.write("</detector-grouping>")        
"""    
