import numpy as np
import matplotlib.pyplot as plt
from simple_pickle_loader import SimplePickleNuclearDatabase as nucleardb

def get_cross_section_at_energy(reaction, target_energy: float) -> float:
    """
    get cross section at specific energy using linear interpolation, could improve in future w/ better interpolation suited for logarithmic spacing etc
    """
    if not getattr(reaction, 'has_tabulated_data', False):
        raise RuntimeWarning("No tabulated data")
        return 0.0
    tabulated_xs = getattr(reaction, 'tabulated_xs', None)
    if not tabulated_xs:
        raise RuntimeWarning("No tabulated cross sections")
        return 0.0
        
    for temp, (energies, xs_values) in tabulated_xs.items():
        if not energies or not xs_values:
            continue
        min_energy, max_energy = min(energies), max(energies)
        if target_energy < min_energy or target_energy > max_energy:
            return 0.0
        for i in range(len(energies) - 1):
            if energies[i] <= target_energy <= energies[i + 1]:
                e1, e2 = energies[i], energies[i + 1]
                xs1, xs2 = xs_values[i], xs_values[i + 1]
                if e2 == e1:
                    return xs1
                xs_interpolated = xs1 + (xs2 - xs1) * (target_energy - e1) / (e2 - e1)
                return max(0.0, xs_interpolated)
        if target_energy == energies[0]:
            return xs_values[0]
        elif target_energy == energies[-1]:
            return xs_values[-1]
        break
    return 0.0

db=nucleardb(r"C:\Users\lowin\Documents\Marathon\nuclear-parser-temp\nuclear-parser-temp\My crap")

isotope = db.find_nuclide("C", 12)
allreactions = db.get_reactions_from(isotope)
#print(allreactions)
reactions = [reaction for reaction in allreactions if hasattr(reaction, 'reaction_type') and reaction.reaction_type.value=='n,non']
print(f"No. of {reactions[0].reaction_type.value} reactions: " + str(len(reactions)))

energy_min = 1e-1  # in eV
energy_max = 2e7   # eV (20 MeV)
steps = 300
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), steps)

xs = np.zeros(steps)
for reaction in reactions:
    for i in range(steps):
        xs[i] = get_cross_section_at_energy(reaction, energies[i])

plt.plot(energies, xs)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Neutron Energy (eV)', fontsize=14)
plt.ylabel('Cross Section (barns)', fontsize=14)
plt.title(f'{reactions[0].reaction_type.value} Reaction Cross Sections for {reactions[0].reactant}', fontsize=16)
plt.grid(True, which='both', ls='--', alpha=0.5)
plt.tight_layout()
plt.show()