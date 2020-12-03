###################      MANIPULATING DATA     ###################
# the excited states energies are in relation to the ground state in
# its own geometry, for the plot is necessary to find the lowest gound
# state energy, set it to zero and calculate all energies in relation
# to the lowest ground state

def energies(args, energy):
    # sort excitation energies, since orca gives it mixed
    for i in energy.index:
        energy.loc[i, :] = energy.loc[i, :].sort_values().tolist()

    minimal = energy[0].min()

    # dislocates all ground states energies
    energy[0] = energy[0].subtract(minimal)

    test = energy.add(energy[0], axis=0)
    test[0] = energy[0]

    return test.sort_index()