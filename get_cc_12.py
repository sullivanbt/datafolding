import numpy as np
import matplotlib.pyplot as plt
from iotbx.file_reader import any_file
mtzfile1 = '/home/ntv/Desktop/writemtz/out.mtz'
mtzfile2 = '/home/ntv/Desktop/second_xtal_normalization/beta_lac_merged.mtz'
m1 = any_file(mtzfile1).file_content.as_miller_arrays()[0]
m2 = any_file(mtzfile2).file_content.as_miller_arrays()[0]
m1.sort(by_value='packed_indices')
m2.sort(by_value='packed_indices')
indices_1 = list(m1.indices())
indices_2 = list(m2.indices())

indices_total = indices_1+indices_2
print(len(indices_1), len(indices_2), len(np.unique(indices_total)))


