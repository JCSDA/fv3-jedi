# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.getBackground as generic


class getBackgroundUFS(generic.getBackground):

    def setup(self, config, fc, fix):

        # Get generic defaults
        generic.getBackground.setup(self, config, fc, fix)

        if 'hack_step_bg' in config and config['hack_step_bg'] == True:
            self.RUNTIME_YAML['hack_step_bg'] = True

        # Remove extra info in names of files so it is readable by FMS
        self.output['bg']['datapath'] = os.path.join(self.output['bg']['datapath'], 'INPUT')
        self.output['bg']['filename_core'] = 'fv_core.res.nc'
        self.output['bg']['filename_trcr'] = 'fv_tracer.res.nc'
        self.output['bg']['filename_sfcd'] = 'sfc_data.nc'
        self.output['bg']['filename_sfcw'] = 'fv_srf_wnd.res.nc'
        self.output['bg']['filename_cplr'] = 'coupler.res'

        # Needed to link static data
        self.RUNTIME_YAML['fcworkdir'] = self.workdir['wdir']
        self.RUNTIME_YAML['ufs_modeldir'] = fix['ufs_modeldir']
        self.RUNTIME_YAML['fv3repo'] = os.path.join(os.environ.get("JEDI_SRC"), "fv3-jedi")

        self.RUNTIME_YAML['fc_length'] = config['forecast_length']
        self.RUNTIME_YAML['fc_freq'] = config['forecast_output_frequency']

        # Set the model run dir the current workdir in main config
        config['MODEL']['ufs_run_directory'] = self.workdir['wdir']

        # Use GFS specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetBackgroundUFS.py")

        self.exec_cmd = ''   # Run on login node for S3 and R2D2 Database access
        self.include_header = ''
        self.login_node_limit = 'True'
