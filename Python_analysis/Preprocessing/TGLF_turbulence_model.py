import numpy as np

# Parameters which were fitted to the (delta_Br/B)**2 from TGLF and
# normalized to yield a mid-plane averaged fluctuation amplitude of 6.5e-8 (in agreement with the RIP diagnostic measurements)
A    = 1.620386820342208e-06
mu   = 0.32309683251166604
sigma= 0.017691893772250868

# The following function is the minimal eample to describe (the delta_Br/B)**2 fluctuations 
# as a function of minor radius [m] and poloidal angle [rad]
def dBr_norm_squared(r, theta, A, mu, sigma):
    g_r = A * np.exp(-0.5 * ((r - mu)/sigma)**2)
    ballooning = 0.25 * (1.0 + np.cos(theta))**2
    return g_r * ballooning


######################################################################################

# All the following is to visualize this analytic description of δBr

import matplotlib.pyplot as plt

# First, the normalized toroidal flux coordinate ρ_t (for reference only here, not used)
rhot = np.arange(0.1, 0.95, 0.05)

# Corresponding minor radii r(ρ_t) in meters from the equilibrium.
radii = np.array([
    0.06591262, 0.09884386, 0.13175743, 0.16459691, 0.19729905, 0.22977723,
    0.26193569, 0.29365476, 0.32480509, 0.3552533 , 0.38487346, 0.41351341,
    0.44098735, 0.46713003, 0.49172032, 0.51454623, 0.5349404
])

# From the equilibrium reconstruction, the elongation κ(r). We use this to
# construct R–Z contours that reproduce the actual plasma shape better
kappa_r = np.array([
    1.22291, 1.22609, 1.2297 , 1.23414, 1.23984, 1.24711, 1.25581, 1.26617,
    1.27842, 1.29304, 1.31025, 1.33048, 1.35444, 1.38288, 1.417  , 1.45834,
    1.51087
])


# Now, ronstruct a (r, θ) grid that will be mapped to (R, Z)

# Radial extent of the visualization in minor radius [m], 
# taken from the experimental equilibrium reconstruction
r_min = 0
r_max = 0.55

# Grid resolution:
nr_fine  = 200
ntheta   = 200

# Radial grid (1D) in minor radius r.
r_fine   = np.linspace(r_min, r_max, nr_fine)

# Poloidal angle grid (1D), measured from the outboard midplane (θ = 0).
theta    = np.linspace(0.0, 2.0*np.pi, ntheta)

# Interpolate κ(r) onto the fine radial grid r_fine
kappa_fine = np.interp(r_fine, radii, kappa_r)

# Create 2D arrays:
#   R  : minor radius r on the 2D grid, shape (nr_fine, ntheta)
#   TH : poloidal angle θ on the 2D grid, same shape
R, TH = np.meshgrid(r_fine, theta, indexing="ij")

# Build κ(r, θ) by simply calculating κ(r) along the θ direction:
K = np.tile(kappa_fine[:, None], (1, ntheta))           # κ(r, θ)

# Calculate the analytic δBr^2(r, θ) on the (r, θ) grid
deltaBr_2D = dBr_norm_squared(R, TH, A, mu, sigma)

# Map the (r, θ) grid into cylindrical coordinates (R_cart, Z_cart)
R0 = 1.7   # major radius [m]
R_cart = R0 + R * np.cos(TH)
Z_cart =      K * R * np.sin(TH)

# Plot δBr(r, θ) in the poloidal (R, Z) plane
fig1, ax1 = plt.subplots(figsize=(5, 5))

mesh = ax1.pcolormesh(
    R_cart,
    Z_cart,
    deltaBr_2D,
    shading='auto'
)

ax1.set_aspect('equal') 
ax1.set_xlabel('R (m)')
ax1.set_ylabel('Z (m)')
ax1.set_title(r'$(\delta B_r/B)^2$ in poloidal plane')

# Colorbar to show the magnitude of δBr on the plot.
cbar = fig1.colorbar(mesh, ax=ax1)
cbar.set_label(r'$(\delta B_r/B)^2$ (arb. units)')

plt.tight_layout()
plt.show()