import pickle
import numpy as np
import os
import nuclear_chain_analyzer as NCA # type: ignore
import copy

datareadpath = r'C:\Users\lowin\Documents\Marathon\nuclear-parser-temp\nuclear-parser-temp\My crap\jendl_nuclear_data_appended.pkl'
datawritepath = r'C:\Users\lowin\Documents\Marathon\nuclear-parser-temp\nuclear-parser-temp\My crap\jendl_nuclear_data_appended.pkl'
txtdatawritepath = r'C:\Users\lowin\Documents\Marathon\nuclear-parser-temp\nuclear-parser-temp\My crap\jendl_nuclear_data_appended.txt'

##list of relevant nuclides
Li6 = NCA.Nuclide.from_identity(z=3, a=6, m=0)
T = NCA.Nuclide.from_identity(z=1, a=3, m=0)
C12 = NCA.Nuclide.from_identity(z=6, a=12, m=0)
Hg198 = NCA.Nuclide.from_identity(z=80, a=198, m=0)

# Open the file in binary mode and load the data
with open(datareadpath, 'rb') as file:
    data = pickle.load(file)

##CHANGE THIS EACH TIME##
xsdatapath = r'C:\Users\lowin\Documents\Marathon\nuclear-parser-temp\nuclear-parser-temp\My crap\li6nt.txt'
reactant = Li6
product = T
reactiontype = NCA.ReactionType.NEUTRON_TRITON
reactionvalue = 'n,t' #not used in reaction data, but important for checking for duplicate reactions before data write

first_column = []
second_column = []

with open(xsdatapath, 'r') as file:
    for line in file:
        if line.strip():  # skip empty lines
            first_value = line.split()[0]
            second_value = line.split()[1]
            first_column.append(float(first_value)*1e6)  
            second_column.append(float(second_value))

energies = first_column
xs = second_column

moddata = copy.deepcopy(data)

existsalready = False
for reaction in moddata['reactions']:
    if reaction.reaction_type.value==reactionvalue and reaction.reactant==reactant:
        existsalready = True
        print(reaction.reactant, reaction.reaction_type)
        raise RuntimeError(f"Desired reactant {reactant} already has reaction of desired type: {reactiontype}")
        #same reactant, same reaction type means no data is changed and error raised - trying to avoid writing over existing data/creating duplicates

moddata['reactions'].append(NCA.Reaction(reactant=reactant, 
                                         product=product, 
                                         reaction_type = reactiontype,
                                         tabulated_xs = {'0K': [energies, xs]}, #temp of 0K is not accurate, just mimics existing data
                                         has_tabulated_data = True))

with open(datawritepath, "wb") as f:
    pickle.dump(moddata, f)
with open(txtdatawritepath, "w") as f: #outputs human readable text for sanity checking
    f.write(str(moddata))