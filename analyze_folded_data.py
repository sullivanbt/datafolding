import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import pickle
import glob
import itertools
from tqdm import tqdm

#==========================================
# Load the data and create a dataframe
#==========================================
data_dir = '/home/ntv/Desktop/ml_results/'
run_nums = [9113,9114,9115,9116,9117]
#run_nums = range(9382,9390+1)


runNumbersList = []
for i in range(1,len(run_nums)+1):
    for v in list(itertools.combinations(run_nums, i)):
        runNumbersList.append(np.array(v))


numBad = 0
totalMissing = 0
dfTmpList = []
badRuns = []
for runNumbers in runNumbersList:
    try:
        paramFileName = '/SNS/users/ntv/Desktop/data_folding/ml_results/beta_lac_secondcrystal_combined_folded_paramList_%s_v2.pkl'%runNumbers
        #paramFileName = '/SNS/users/ntv/Desktop/data_folding/ml_results/nak2019_combined_folded_paramList_%s_v2.pkl'%runNumbers
        paramFileName = paramFileName.replace(' ', ',_').replace('[','').replace(']','')
        dfFileName = '/SNS/users/ntv/Desktop/data_folding/ml_results/beta_lac_secondcrystal_combined_folded_%s_v2.pkl'%runNumbers
        #dfFileName = '/SNS/users/ntv/Desktop/data_folding/ml_results/nak2019_combined_folded_%s_v2.pkl'%runNumbers
        dfFileName = dfFileName.replace(' ', ',_').replace('[','').replace(']','')
        paramList = pickle.load(open(paramFileName, 'rb'))
        dfTmp = pd.read_pickle(dfFileName)   
        print('%s: len(df): %i, len(paramList): %i'%(str(runNumbers), len(dfTmp), len(paramList)))
        if len(dfTmp)!=len(paramList):
            badRuns.append(runNumbers)
            totalMissing += (len(dfTmp)-len(paramList))
        #for key in dfTmp2.keys():
        #    dfTmp[key] = dfTmp2[key]        
        dfTmpList.append(dfTmp)
    except:
        numBad += 1
        print 'Error with %s'%str(runNumbers)
        print '    %s'%paramFileName
        print '    %s'%dfFileName
df = pd.concat(dfTmpList)

#==========================================
# write an output file
#==========================================
numRows = df.shape[0]

"""
nanIDX = pd.isnull(df['IntensML'])
df.loc[nanIDX,'IntensML']=0.
df.loc[nanIDX,'SigIntML']=1.
"""

d = {}
d['h'] = []
d['k'] = []
d['l'] = []
d['Intens'] = []
d['SigInt'] = []


gIDX = ~(pd.isnull(df['IntensML'])) &(df['IntensML']/df['SigIntML']>-3.)
print('gIDX.sum(): ', gIDX.sum())
for i in tqdm(range(numRows)):
    if gIDX[i]:
        ind = df.index[i]
        ind = tuple(np.array(ind.strip('()').split(','), dtype=float))
        d['h'].append(ind[0])
        d['k'].append(ind[1])
        d['l'].append(ind[2])
        
        d['Intens'].append(df.iloc[i]['IntensML'])
        d['SigInt'].append(df.iloc[i]['SigIntML'])

pickle.dump(d, open('/home/ntv/Desktop/data_folding/integrated_intensities.pkl','wb')) 
    





