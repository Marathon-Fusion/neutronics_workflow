import os
from pprint import pprint
import shutil
import subprocess
import urllib.request

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.patches import Rectangle

import openmc.data

# url = 'https://anl.box.com/shared/static/kxm7s57z3xgfbeq29h54n7q6js8rd11c.ace'
# filename, headers = urllib.request.urlretrieve(url, 'c12.ace')
isotope = openmc.data.IncidentNeutron.from_ace(r'HG198.ace')
nel = isotope[2] #MT2 = elastic neutron collision
neutron = nel.products[0]
distribution = neutron.distribution[0]
#print([angle for angle in distribution.angle])

# energies = c12.energy['300K']
# c12nelxs = c12nel.xs['300K'](energies)
# plt.loglog(energies, c12nelxs)
# plt.xlabel('Energy (eV)')
# plt.ylabel('Cross section (b)')
# plt.savefig("/opt/marathon/workspace/Matplotlib_output_dump/latest_data_reader_plot.png")

# # print the reactions
# pprint(list(c12.reactions.values())[:10])

# (n,2n)
n2n = isotope[16]
print('Threshold = {} eV'.format(n2n.xs['294K'].x[0]))
xs = n2n.xs['294K']
plt.plot(xs.x, xs.y)
plt.xlabel('Energy (eV)')
plt.ylabel('Cross section (b)')
plt.xlim((xs.x[0], xs.x[-1]))
plt.savefig("/opt/marathon/workspace/Matplotlib_output_dump/latest_data_reader_plot_xs.png")


#energy distribution out
neutron = n2n.products[0]
dist = neutron.distribution[0]
# for e_in, e_out_dist in zip(dist.energy, dist.energy_out):
#     plt.semilogy(e_out_dist.x, e_out_dist.p, label='E={:.2f} MeV'.format(e_in/1e6))
distno = 0
for energy_out_dist in dist.energy.energy_out:
    #plt.semilogy(energy, energy_out, label='E={:.2f} MeV')#.format(e_in/1e6))
    distno += 1
    print(distno)
    samples = energy_out_dist.sample(n_samples=100)
    plt.hist(samples, bins=100, density=True, alpha=0.7, label=f'dist {distno}')
plt.xlabel('Value')
plt.ylabel('Density')
plt.title(' Distribution (Sampled)')
plt.legend()
plt.savefig("/opt/marathon/workspace/Matplotlib_output_dump/latest_data_reader_plot_output_energies.png")
# plt.ylim(top=1e-6)
# plt.legend()
# plt.xlabel('Outgoing energy (eV)')
# plt.ylabel('Probability/eV')
# plt.savefig("/opt/marathon/workspace/Matplotlib_output_dump/latest_data_reader_plot_output_energies.png")