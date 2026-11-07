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

# Depending on how KORC is configured (cylindrical vs Cartesian restart),
# X can pack [R, Z, PHI] or cartesian [X, Y, Z]. Assuming cylindrical/flux coordinates
# or direct mapping into the (3, N) layout expected by X:
X_data = np.vstack((R, Z, PHI))  # Shape: (3, num_particles)


vmag=c*0.6
g=1/np.sqrt(1-(vmag/c)**2)
eta=0. #in degrees

# Placeholder velocities/momenta array of shape (3, num_particles)
V_data = np.zeros((3, num_particles), dtype=np.float64)

# 2. Open template H5 to copy global attributes/metadata, or create a new one
with h5py.File(template_h5_path, "r") as f_templ:
  with h5py.File(output_h5_path, "w") as f_out:
    # Copy global simulation datasets
    for key in f_templ.keys():
      if key != "spp_1":
        f_templ.copy(key, f_out)

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