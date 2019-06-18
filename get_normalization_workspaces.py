import numpy as np
# Setting up the workspaces containing information about the flux and the solid angle (from a vanadium run)
rawVan=Load(Filename='/SNS/MANDI/IPTS-8776/nexus/MANDI_8905.nxs.h5')
LoadIsawDetCal(Inputworkspace=rawVan, Filename='/home/ntv/mandi_preprocessing/MANDI_June2018.DetCal')

lam_min = 1.9
lam_max = 4.05
p_min = 2*np.pi/lam_max
p_max = 2*np.pi/lam_min

rawVan=ConvertUnits(InputWorkspace=rawVan,Target='Momentum')
MaskBTP(Workspace=rawVan,Bank="13")
#MaskBTP(Workspace=rawVan,Pixel="0-9,246-255")
#MaskBTP(Workspace=rawVan,Tube="0-9,246-255")
rawVan=CropWorkspace(InputWorkspace=rawVan,XMin=p_min,XMax=p_max)

#Solid angle
sa=Rebin(InputWorkspace=rawVan,Params='{},{},{}'.format(p_min, (p_max-p_min), p_max),PreserveEvents='0')
SaveNexus(InputWorkspace=sa, Filename="/home/ntv/Desktop/data_folding/MANDI_solidAngle_nomask.nxs")


#flux
rawVan = Rebin(InputWorkspace=rawVan, Params='{},{},{}'.format(p_min, (p_max-p_min)/100., p_max))
flux=GroupDetectors(InputWorkspace=rawVan,MapFile='/home/ntv/Desktop/data_folding/mandi_grouping.xml')
#DeleteWorkspace(rawVan)
flux=CompressEvents(InputWorkspace=flux,Tolerance=1e-4)
flux=Rebin(InputWorkspace=flux,Params='{},{},{}'.format(p_min, (p_max-p_min), p_max))
for i in range(flux.getNumberHistograms()):
  try:
      el=flux.getSpectrum(i)
      el.divide(flux.readY(i)[0],0)
  except:
      pass
flux=Rebin(InputWorkspace=flux,Params='{},{},{}'.format(p_min, (p_max-p_min), p_max))
flux_cdf=IntegrateFlux(flux)
SaveNexus(InputWorkspace=flux_cdf, Filename="/home/ntv/Desktop/data_folding/spectra.nxs")

