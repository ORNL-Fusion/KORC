import h5py
import numpy as np
import pandas as pd

mu0=4*np.pi*10**(-7)
c = 2.99792458E8 #% Speed of light (m/s)
qe = 1.60217662E-19 #% Electron charge (C)
me = 9.10938356E-31 #% Electron mass (kg)
ep0 = 8.854E-12 #% Electric permitivity [C**2/(N*m**2)]
re = qe**2/(4*np.pi*ep0*me*c**2) #% classical electron radius (m)

csv_filepath="start_pts (1).csv"
template_h5_path="restart_file_tmp.h5"
output_h5_path="restart_file.h5"


# 1. Load particle positions from CSV
print(f"Reading particle locations from {csv_filepath}...")
df = pd.read_csv(csv_filepath)

R = df["R"].values
Z = df["Z"].values
PHI = df["phi"].values
num_particles = len(R)
print(f"Loaded {num_particles} particles from CSV.")

vmag=c*0.6
g=1/np.sqrt(1-(vmag/c)**2)
eta=0. #in degrees

# Assuming nmpi_prev (number of MPI processes the restart file was saved with)
nmpi_prev = 1  # Match whatever nmpi was used to generate the restart file
ppp = num_particles // nmpi_prev

if num_particles % nmpi_prev != 0:
    raise ValueError("Total particles must be divisible by nmpi_prev!")

# Build X_data and V_data rank-by-rank with component-major ordering:
# [R_rank0, PHI_rank0, Z_rank0, R_rank1, PHI_rank1, Z_rank1, ...]
X_list = []
V_list = []

for ii in range(nmpi_prev):
    start_idx = ii * ppp
    end_idx = (ii + 1) * ppp
    
    # Coordinates for rank ii
    X_list.extend(R[start_idx:end_idx])
    X_list.extend(PHI[start_idx:end_idx])
    X_list.extend(Z[start_idx:end_idx])
    
    # Velocities/Momenta placeholders for rank ii (zeros or actual values)
    V_list.extend(np.zeros(ppp))  # V1 component
    V_list.extend(np.zeros(ppp))  # V2 component
    V_list.extend(np.zeros(ppp))  # V3 component

X_data = np.array(X_list, dtype=np.float64)
V_data = np.array(V_list, dtype=np.float64)

# 2. Open template H5 to copy global attributes/metadata, or create a new one
with h5py.File(template_h5_path, "r") as f_templ:
  with h5py.File(output_h5_path, "w") as f_out:
    # Copy global simulation datasets
    for key in f_templ.keys():
      if key != "spp_1" and key != "num_mpi" :
        f_templ.copy(key, f_out)

    f_out.create_dataset(
        "num_mpi", data = nmpi_prev, dtype=np.int32
    )

    # Create spp_1 group
    spp_out = f_out.create_group("spp_1")

    # Write X and V datasets in (3, N) shape
    spp_out.create_dataset("X", data=X_data, dtype=np.float64)
    spp_out.create_dataset("V", data=V_data, dtype=np.float64)
    
    # Create default flags and weights arrays of length num_particles
    spp_out.create_dataset(
        "flagCon", data=np.ones(num_particles, dtype=np.int32)
    )
    spp_out.create_dataset(
        "flagCol", data=np.ones(num_particles, dtype=np.int32)
    )
    spp_out.create_dataset(
        "flagRE", data=np.ones(num_particles, dtype=np.int32)
    )
    spp_out.create_dataset(
        "g", data=np.ones(num_particles, dtype=np.float64)*g
    )
    spp_out.create_dataset(
        "eta", data=np.ones(num_particles, dtype=np.float64)*eta
    )

print(
    f"Successfully generated KORC restart file at: {output_h5_path} with"
    f" {num_particles} particles."
    )