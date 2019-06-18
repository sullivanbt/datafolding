import iotbx.mtz
import cctbx.sgtbx
import numpy as np
import cctbx.array_family as af
import pickle

#cob = cctbx.sgtbx.change_of_basis_op('x,y,-z')
#sg = cctbx.sgtbx.space_group("P 32 2 1")
#sg = cctbx.sgtbx.space_group('P 2ac 2ab')
#oldops = [sg.all_ops()[i].as_xyz() for i in range(sg.n_smx())]
#sg = sg.change_basis(cob)
#newops = [sg.all_ops()[i].as_xyz() for i in range(sg.n_smx())]
#print oldops
#print newops

"""
hIn = []; kIn = []; lIn = []
Iin = []; SigIn = [];
laueNormPath = '/home/ntv/Desktop/second_xtal_normalization/'
numFiles = 5
for fileNumber in range(numFiles):
    with open(laueNormPath+'laueNorm%03i'%(fileNumber+1), 'r') as f:
        lines = f.readlines()
    for line in lines:
        h = int(line[:5])
        k = int(line[5:10])
        l = int(line[10:15])
        lam = float(line[15:25])
        theta = float(line[25:35])
        intens = int(line[35:45])
        sigInt = int(line[45:55])
        if sigInt > 0:
            if intens/sigInt > -3.0:
                hIn.append(h)
                kIn.append(k);
                lIn.append(l);
                Iin.append(intens);
                SigIn.append(sigInt)
"""
d = pickle.load(open('/home/ntv/Desktop/data_folding/integrated_intensities.pkl','rb'))
hIn = d['h']
kIn = d['k']
lIn = d['l']
Iin = d['Intens']
SigIn = d['SigInt']

#define the mtz object and set the space group properties
"""
sg = cctbx.sgtbx.space_group("P 1")
m = iotbx.mtz.object()
m.add_crystal("BS XTAL", "TEST", (73.4, 73.4, 99.39, 90, 90, 120))
m.set_title("TEST MTZ")
m.set_space_group_name('P 1')
m.set_point_group_name('1')
m.set_lattice_centring_type('P')
m.set_space_group_number(1)
m.set_space_group(sg)
"""



#beta lac second xtal
sg = cctbx.sgtbx.space_group("P 32 2\"")
m = iotbx.mtz.object()
m.add_crystal("BS XTAL", "TEST", (73.4, 73.4, 99.39, 90, 90, 120))
m.set_title("TEST MTZ")
m.set_space_group_name('P 32 2 1')
m.set_point_group_name('321')
m.set_lattice_centring_type('P')
m.set_space_group_number(154)
m.set_space_group(sg)
"""
#nak 2019
sg = cctbx.sgtbx.space_group("I 4")
m = iotbx.mtz.object()
m.add_crystal("BS XTAL", "TEST", (69.1617, 69.1617, 93.19, 90, 90, 90))
#m.add_crystal("BS XTAL", "TEST", (68.47, 69.33, 93.20, 90, 90, 90))
m.set_title("TEST MTZ")
m.set_space_group_name('I 4')
m.set_point_group_name('4')
m.set_lattice_centring_type('I')
m.set_space_group_number(79)
m.set_space_group(sg)
"""


m.set_n_reflections(len(hIn))

cryst = iotbx.mtz.crystal(m, 0)
cryst.add_dataset('ds1', 0.0)
ds = iotbx.mtz.dataset(cryst, 0)
hcol = ds.add_column("H", "H")
kcol = ds.add_column("K", "H")
lcol = ds.add_column("L", "H")
Icol = ds.add_column("I", "J")
Sigcol = ds.add_column("SIGI", "Q")


h = af.flex.float(hIn)
k = af.flex.float(kIn)
l = af.flex.float(lIn)
I = af.flex.float(Iin)
Sig = af.flex.float(SigIn)
ff = af.flex.bool([True]*len(hIn))
hcol.set_values(h, ff)
kcol.set_values(k, ff)
lcol.set_values(l, ff)
Icol.set_values(I, ff)
Sigcol.set_values(Sig, ff)
m.write("/home/ntv/Desktop/writemtz/out.mtz")

