# G12: [1,2]; G15: [8,9,10, 17]
import numpy as np

desi_fields     = np.array(['R1', 'R2', 'R8', 'R9', 'R10', 'R17'])
desi_allfields  = np.array(['R{}'.format(x) for x in range(20)])