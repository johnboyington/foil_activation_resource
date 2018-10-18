import numpy as np
from multigroup_utilities import energy_groups
from reaction import build_foil_library


def generate_response_matrix(foil_names, structure='wims69', source='trigaC', filename='response_matrix'):
    """This function generates a response matrix given a set of foil ids."""

    # grab energy groups
    eb = energy_groups(structure)[::-1]

    # initialize response matrix
    response_matrix = np.empty((len(foil_names), len(eb[1:])))

    # build foil library
    foil_lib = build_foil_library('trigaC')

    # calculate individual foil response functions
    for i, foil_name in enumerate(foil_names):
        response_matrix[i], binned_flux = foil_lib[foil_name].discretize(eb)

    # store data
    np.save(filename + '.npy', response_matrix)
    np.save(source + '.npy', binned_flux)

    return response_matrix


if __name__ == '__main__':
    foil_names = ['In', 'InCd', 'Au', 'AuCd', 'Mo', 'MoCd', 'Al']
    generate_response_matrix(foil_names)
