import numpy as np
from wdata.io import WData, Var

fwtxt="test.wtxt"
data = WData.load(fwtxt)

# read from wdat
density_a = np.array(data.density_a)
current_a = np.array(data.current_a)
delta = np.array(data.delta)

# save
np.save(data.prefix+"_rho_a", density_a)
np.save(data.prefix+"_j_a", current_a)
np.save(data.prefix+"_psi", delta)

# add to wtxt new entries
f = open(data.prefix+".wtxt", "a")
f.write("var               rho_a                    real                    none                    npy\n")
f.write("var                 j_a                  vector                    none                    npy\n")
f.write("var                 psi                 complex                    none                    npy\n")
f.close()
