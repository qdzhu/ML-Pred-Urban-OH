#!/Users/Josh/Documents/Python/VirtualEnvs/PECANS/bin/python
import os
import numpy as np

from pecans.ensembles.api import EnsembleRunner
from pecans.utilities.config import load_config_file


def name_output_files(index, **config_opts):
    lifetime_hours = config_opts['CHEMISTRY/mechanism_opts/lifetime_seconds'] / 3600
    emissions_width_km = config_opts['EMISSIONS/emission_opts/width_x'] / 1000
    return 'pecans_ens_tau-{}h_emwidth-{}km'.format(lifetime_hours, emissions_width_km)


_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
config_file = os.path.join(_mydir, 'pecans_config.cfg')

# We want lifetimes that vary from 1-9 hours. This covers about the most extreme values we'd expect for summer NOx
# lifetime
taus = np.arange(3600, 9*3600+1, 3600)

# We also want to test what happens when emissions widths are similar or greater than lifetimes. So we'll calculate
# emissions widths equal to each expected lifetime
config = load_config_file(config_file)
winds = config.get('TRANSPORT', 'wind_speeds')
x_wind = winds['x']
widths = taus * x_wind
widths = np.concatenate(([3000], widths))  # add a smaller width as an extra test

ens = EnsembleRunner(config_file,
                     ensemble_variables={'CHEMISTRY/mechanism_opts/lifetime_seconds': taus,
                                         'EMISSIONS/emission_opts/width_x': widths},
                     ensemble_mode='combinations',
                     save_in_individual_dirs=False,
                     save_final_output_only=True,
                     member_naming_fxn=name_output_files,
                     root_output_dir=os.path.join(_mydir, '..', 'Workspaces', 'PECANS', 'lifetime-ensemble'))

ens.run()
