"""

jason@marathonfusion.com

This code is a demo solver for steady state neutron transport equations 
that solve for the neutron distribution function as a function of x

Neutron transport equations for three energy bins in a slab and a single 198-Hg layer with a fixed down scattering rate ξ

definitions of stuff:

G: number of energy bins for neutron distribution function

Σ_abs: neutron absorption cross section of energy group g. Units: cm^{-1}.  counts every reaction that removes a neutron from the chain ( n,gamma) , (n,alpha), etc.).

Σ_tot: total macroscopic cross section of group g

Σ_t: n,t macroscopic cross section

ξ: (TOY MODEL ONLY) chosen fraction of 
non-absorptive collisions that scatter down 
to the next lower energy group. In future we
will omit ξ entirely and calculate / read in 
the full Σ_s (g'-->g) scattering matrix from the 
cross section libraries

Σ_s: scattering matrix (size GxG). each element is the macroscopic cross-section for a neutron that was born in energy group src to re-appear in energy group dest after an elastic or inelastic collision. Units: cm^{-1}.convention used in code: row = destination, column = source

σ_scat: per-group scattering loss used inside the loop, σ_scat = Σ_t - Σ_abs

Σ_2n: macro (n,2n) or fission multiplication cross section of group g. Units: cm^{-1}. zero in this Hg demo, but we can easily include. If non-zero, each collision produces ν-1 extra neutrons.

ν: integer neutron multiplication value for (n,νn) reactions. e.g. for (n,2n) ν = 2

A: the streaming/removal coefficient matrix in the spatial ODE system d phi/d x + A phi = 0. A = diag(Σ_t) - Σ_s.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

AVOGADRO = 6.02214076e23          # mol‑1
BARN     = 1e-24                  # cm²

# 1. energy bins: choose the number of bins G and make it logarithimic
G = 6 # Number of bins
E_edges = np.logspace(np.log10(14.0), np.log10(1e-8), G+1)  # 14 MeV to 0.01 eV, ends up at 0.01 eV, 0.33eV, 11eV, 370eV, 1.25keV, 0.42MeV, 14MeV for 6 bins 
                                                            # 116K, 3800K, 4.2e6K, etc.      
E_mid   = 0.5*(E_edges[:-1] + E_edges[1:])                  # bin centers

L = 35.0
dx = 0.05 # length, grid spacing
Nx = int(L/dx) + 1 # number of grid points
x = np.linspace(0, L, Nx) #array of x values to use

# 2. macroscopic cross sections; for now these are TOY cross sections for each energy bin but in future we will actually add database cross sections

#MULTIPLIER LAYER##

rho_Hg198 = 13.6                     # g cm⁻³
Mr_Hg198   = 198
n_Hg198   = rho_Hg198 * AVOGADRO / Mr_Hg198

Hg198_Σ_tot = np.array([6, 10.0, 30.0, 20.0, 13.0, 13.0]) * n_Hg198 * BARN
Hg198_Σ_abs = np.array([2.1, 0.15, 0.3, 0.5, 0.3, 1.0]) * n_Hg198 * BARN
Hg198_Σ_2n = np.array([2.0, 0.0, 0.0, 0.0, 0.0, 0.0]) * n_Hg198 * BARN
Hg198_Σ = np.vstack((Hg198_Σ_tot, Hg198_Σ_abs, Hg198_Σ_2n))

##MODERATOR LAYER##

rho_C12 = 2.20
Mr_C12   = 12
n_C12   = rho_C12 * AVOGADRO / Mr_C12 

C12_Σ_tot = np.array([2.0, 4.0, 5.0, 5.0, 5.0, 5.0]) * n_C12 * BARN
C12_Σ_abs = np.array([0.1, 0.0, 0.0, 0.0, 0.0, 0.0]) * n_C12 * BARN
C12_Σ = np.vstack((C12_Σ_tot, C12_Σ_abs))

##BREEDER LAYER##

n_Li6 = 4.0e22

Li6_Σ_tot = np.array([3.0, 3.0, 4.0, 20.0, 100.0, 500.0]) * n_Li6 * BARN #first group is highest energy
Li6_Σ_abs = np.array([0.4, 1.0, 3.0, 20.0, 100.0, 500.0]) * n_Li6 * BARN
Li6_Σ_t = np.array([0.2, 1.0, 3.0, 20.0, 100.0, 500.0]) * n_Li6 * BARN
Li6_Σ = np.vstack((Li6_Σ_tot, Li6_Σ_abs, np.zeros((1, G)),Li6_Σ_t)) #row of zeroes to represent no n,2n reactions

ξ = 0.2          # fraction of non-absorbed collisions that scatter DOWN (nothing can scatter up)
ν = 2.0 # (neutron multiplication factor for (n,2n))

# build the scattering matrix (dest row, src col)
def buildscatterm(cross_sections):
    Σ_s = np.zeros((G, G)) # scattering matrix
    for src in range(G):
        Σ_scat = cross_sections[0][src] - cross_sections[1][src] # all scattering out of the source group into another group

        if src < G-1:
            down = ξ * Σ_scat # src -> lower group# TOY MODEL, this is our simple model for the scattering. to replace with actually correct libraries
            stay = Σ_scat - down # remainder stays
            Σ_s[src+1, src] += down
        else:
            stay = Σ_scat # last group: all stay
        Σ_s[src, src] += stay
    return Σ_s

def buildn2nm(cross_sections):
    Σ_2n = np.zeros((G, G))
    if np.shape(cross_sections)[0] > 1: #if n,2n reactions present
        Σ_2nxs = cross_sections[2]
        for src in range(G):
            if src == 0: #only highest energy bin can do n,2n in Hg198
                # down = ν * Σ_2nxs[src]
                # Σ_2n[src+1, src] = down #n,2n reaction spawns 2 neutrons at the 2nd highest bin
                Σ_2n[src, src] -= Σ_2nxs[src]
                Σ_2n[src+1, src] += ν*Σ_2nxs[src]
    return Σ_2n

# 3. streaming/removal matrx
    # A_Hg198 = np.diag(Hg198_Σ[0]) - buildscatterm(Hg198_Σ) - np.diag(ν*Hg198_Σ[2])#−buildscatterm gives +source on RHS, (n,2n) term subtracts (ν)Σ_m on the diagonal
A_Hg198 = np.diag(Hg198_Σ[0]) - buildscatterm(Hg198_Σ) - buildn2nm(Hg198_Σ)
A_C12 = np.diag(C12_Σ[0]) - buildscatterm(C12_Σ)
A_Li6 = np.diag(Li6_Σ[0]) - buildscatterm(Li6_Σ)

# sanity check: column balance (absorbed + all scattered out = total) this has to be enforced...
#assert np.allclose(Σ_abs + Σ_s.sum(axis=0), Σ_tot)
def transport_eqs(x, boundaries=np.array([10.0, 20.0, 5.0])):

    # calculating flux at next timestep with analytic solution
    prop = lambda φ, Δx, A: expm(-A*Δx) @ φ # this is our model for calculating flux at the next spatial grid point φ(x+Δx) = e^(−AΔx) φ(x), only possible because A is constant w.r.t. x!

    #  4. slab & mesh parameters.
    x1 = boundaries[0]
    x2 = boundaries[1]
    x3 = boundaries[2]

    def region(x):
        if x <= x1:
            return 'multiplier'
        elif x1 < x and x <= x1 + x2:
            return 'moderator'
        elif x1 + x2 < x:
            return 'breeder'

    # intializing the flux
    φ = np.zeros((Nx, G))
    φ[0, 0] = 1.0 # 14 MeV source only at the outer boundary!

    # find the flux in each region.
    for i in range(1, len(x)):
        if region(x[i]) == 'multiplier':
            φ[i] = prop(φ[i-1], dx, A_Hg198)
        elif region(x[i]) == 'moderator':
            φ[i] = prop(φ[i-1], dx, A_C12)
        else:
            φ[i] = prop(φ[i-1], dx, A_Li6)

    sum_n_out = 0
    for i in range(G):
        for j in range(Nx):
            if region(x[j]) == 'breeder':
                sum_n_out += φ[j][i]*Li6_Σ[3][i]*dx

    TBR = sum_n_out / sum(φ[0])

    return φ, TBR

def find_x_max(x, φ):
    φ_total = np.sum(φ, axis=1)
    φ_max = np.amax(φ_total, axis=0)
    for i in range(np.shape(φ)[0]):
        if φ_total[i] == φ_max:
            x_max = x[i]

    return x_max, φ_max

def plot_fluxes_and_energy(x, φ):
    # 5. plot stuff
    # plot the total flux

    plt.figure(figsize=(6,4))
    for g in range(G):
        plt.semilogy(x, φ[:,g], label=f"group {g+1}")
    plt.semilogy(x, np.sum(φ, axis=1), label='total', color='black')
    plt.xlabel("x [cm]"); plt.ylabel("$\\varphi_g$ [n cm⁻² s⁻¹]")
    plt.ylim(bottom = 1e-3, top = 1e1)
    plt.title("group fluxes in a combined slab")
    plt.grid(True, which='both')
    plt.legend()

    # neutron energy
    plt.figure(figsize=(6,4))
    plt.plot(x, (φ*E_mid).sum(axis=1))
    plt.xlabel("x [cm]"); plt.ylabel("energy density [MeV cm⁻² s⁻¹]")
    plt.title("neutron energy vs depth")
    plt.grid(True)
    plt.tight_layout()
    #plt.show()

def main():
    breedwidths = np.array([2.0])
    TBRarray = np.zeros((len(breedwidths), Nx))
    for j in range(len(breedwidths)):
        print(f"Scan reached breeder width of {breedwidths[j]}cm")
        multwidth = 0
        for i in range(Nx):
            if multwidth <= L - breedwidths[j]:
                boundaries = np.array([multwidth, L-breedwidths[j]-multwidth, breedwidths[j]])
                φ, TBR = transport_eqs(x, boundaries)
                if np.isclose(multwidth, 10.0, atol=0.001):
                    φ_at_10 = φ
                # x_max, φ_max = find_x_max(x, φ)
                # print(f"Highest relative flux reached: {φ_max} at x = {x_max}cm")
                # print(f"TBR ≈ {np.around(TBR, 2)}")
                TBRarray[j][i] = TBR
                multwidth += dx #evaluation steps spaced 0.5cm apart
                if x[i].is_integer():
                    print(f"Position {x[i]}cm")
        plot_fluxes_and_energy(x, φ_at_10)
        plt.figure(figsize=(6,4))
        plt.plot(x, TBRarray[j], label=f'{breedwidths[j]}')
    plt.xlabel("Mult/mod boundary position (cm)")
    plt.ylabel("TBR")
    plt.legend(title="Breeder layer thickness (cm)")
    plt.title(f"TBR vs. blanket layer thickness")
    #plt.xlim(right=L-breedwidth)
    plt.show()

if __name__ == '__main__':
    main()