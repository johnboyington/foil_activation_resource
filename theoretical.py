'''
An object used to simulate a theoretical experiment.
'''
import subprocess
from wand import Wand
import pickle
from reaction import build_foil_library
from tex_experiment import experiment_template
from amalgamated_plot import amalgamate
import sys
import numpy as np


class Theoretical(object):
    def __init__(self, experimentname, irradiations, detector='ksu', source='trigaC'):
        self.experimentname = experimentname
        self.dataname = 'theoretical_data/theoretical_' + experimentname + '.txt'
        self.detector = detector
        self.source = source
        self.foil_library = build_foil_library(self.source)
        try:
            with open(self.dataname, 'rb') as F:
                self.data = pickle.load(F)
        except:
            self.data = {}

        self.foils_irradiated = []
        for data in irradiations:
            wand = Wand(*data, detector, source, self.foil_library, experimentname)
            wand.irradiate(True)
            self.data.update(wand.package_data())
            self.foils_irradiated.append(data[1])

        rep_power = irradiations[0][-1]  # grabs the first power level used
        amalgamate(self.experimentname, self.foils_irradiated, self.foil_library, self.source, rep_power)
        experiment_temp = experiment_template.split('SPLIT')
        foil_s = ''
        for f in self.foils_irradiated:
            foil_s += '\\include{{plot/{}}}\n'.format(self.foil_library[f].plotname + '_' + experimentname)
        plot_s = 'plot/amalgamated_' + self.experimentname
        experiment_tex = experiment_temp[0] + plot_s + experiment_temp[1] + foil_s + experiment_temp[2]
        experiment_filename = self.experimentname + '.tex'
        with open(experiment_filename, 'w+') as F:
            F.write(experiment_tex)

        # dump data
        with open(self.dataname, 'wb') as F:
            pickle.dump(self.data, F)

        subprocess.run(['pdflatex', experiment_filename])
        subprocess.run(['rm', '{}.log'.format(self.experimentname)])
        subprocess.run(['rm', '{}.aux'.format(self.experimentname)])
        subprocess.run(['rm', '{}.out'.format(self.experimentname)])


if __name__ == '__main__':
    assert len(sys.argv) > 1, 'Must include cases to run.'
    cases = sys.argv[1:]
    for case in cases:
        assert case in ['wisconsin', 'library', 'bismuth'], '{} not an option.'.format(case)

    for case in cases:
        if case == 'wisconsin':
            irradiation_data = []
            irradiation_data.append(['indium', 'In', np.array([2.1, 2.1, 2.1, 2.2]), 60, 3600*5, 120, 100])
            irradiation_data.append(['indium', 'InCd', np.array([2.2, 2.2, 2.0, 2.2]), 60, 3600*4, 120, 100])
            irradiation_data.append(['gold', 'Au', np.array([4.1, 4.1, 4.3, 4.2]), 60, 3600*2, 120, 100])
            irradiation_data.append(['gold', 'AuCd', np.array([4.1, 4.3, 4.1, 4.1]), 60, 3600*2, 120, 100])
            irradiation_data.append(['molybdenum', 'Mo', np.array([2.9, 2.9, 2.9, 2.8]), 3600, 3600*24*1, 600, 100])
            irradiation_data.append(['molybdenum', 'MoCd', np.array([2.7, 2.5, 2.7, 2.7]), 3600, 3600*24*1, 600, 100])
            irradiation_data.append(['rhodium', 'Rh', np.array([5.5, 5.1, 5.1, 5.3]), 3600, 3600*1, 1800, 100])
            irradiation_data.append(['aluminum', 'Al', np.array([3.5, 3.5, 3.5, 3.5]), 3600, 3600*24*1, 3600, 100])
            irradiation_data.append(['titanium', 'Ti1', np.array([5.5, 5.1, 5.1, 5.3]), 3600, 3600*1, 1800, 100])
            irradiation_data.append(['titanium', 'Ti2', np.array([5.5, 5.1, 5.1, 5.3]), 3600, 3600*1, 1800, 100])
            Theoretical('wisconsin', irradiation_data, 'ksu', 'trigaC')

        if case == 'library':
            irradiation_data = []
            irradiation_data.append(['indium', 'In', False, 'estimate', 30, 3600*2, 60, 100])
            irradiation_data.append(['molybdenum', 'Mo', False, 'estimate', 3600, 3600*24*2, 600, 100])
            irradiation_data.append(['zinc', 'Zn', False, 'estimate', 3600*2, 3600*24*2, 3600, 100])
            irradiation_data.append(['copper', 'Cu', False, 'estimate', 3600, 3600*24*2, 3600, 100])
            irradiation_data.append(['magnesium', 'Mg', False, 'estimate', 60, 3600*2, 3600, 100])
            irradiation_data.append(['aluminum', 'Al', False, 'estimate', 3600, 3600*24*2, 3600, 100])
            irradiation_data.append(['calcium', 'Ca', False, 'estimate', 3600, 3600*24*2, 3600, 100])
            irradiation_data.append(['scandium', 'Sc', False, 'estimate', 3600, 3600*24*2, 3600, 100])
            irradiation_data.append(['manganese', 'Mn', False, 'estimate', 60, 3600*2, 600, 100])
            irradiation_data.append(['iron', 'Fe', False, 'estimate', 60, 3600*2, 3600, 100])
            irradiation_data.append(['nickel', 'Ni', False, 'estimate', 3600, 3600*24*2, 3600, 100])
            irradiation_data.append(['zirconium', 'Zr', False, 'estimate', 3600, 3600*24*2, 3600, 100])
            irradiation_data.append(['rhodium', 'Rh', False, 0.6, 60, 3600*2, 3600, 100])
            irradiation_data.append(['dysprosium', 'Dy', False, 'estimate', 60, 3600*2, 120, 100])
            irradiation_data.append(['iridium', 'Ir', False, 'estimate', 60, 3600*2, 600, 100])
            irradiation_data.append(['gold', 'Au', False, 'estimate', 60, 3600*2, 60, 100])
            Theoretical('library', irradiation_data)

        if case == 'bismuth':
            irradiation_data = []
            irradiation_data.append(['bismuth', 'Bi', np.array([1.0, 1.0, 1.0, 1.0])*70, 3600, 1, 3600*24*20, 250])
            Theoretical('bismuth', irradiation_data, 'perfect', 'trigaC')
