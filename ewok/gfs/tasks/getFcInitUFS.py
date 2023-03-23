# (C) Copyright 2020-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import ewok.tasks.getFcInit as generic
import yamltools

class getFcInitUFS(generic.getFcInit):

    def setup(self, config):

        # Get generic defaults
        generic.getFcInit.setup(self, config)

        if 'hack_step_bg' in config and config['hack_step_bg'] == True:
            self.RUNTIME_YAML['hack_step_bg'] = True

        # Remove extra info in names of files so it is readable by FMS
        self.output['fcout']['PT6H']['datapath'] = os.path.join(self.output['fcout']['PT6H']['datapath'], 'RESTART')
        self.output['fcout']['PT6H']['filename_core'] = '{{ufs_current_cycle}}.fv_core.res.nc'
        self.output['fcout']['PT6H']['filename_trcr'] = '{{ufs_current_cycle}}.fv_tracer.res.nc'
        self.output['fcout']['PT6H']['filename_sfcd'] = '{{ufs_current_cycle}}.sfc_data.nc'
        self.output['fcout']['PT6H']['filename_sfcw'] = '{{ufs_current_cycle}}.fv_srf_wnd.res.nc'
        self.output['fcout']['PT6H']['filename_cplr'] = '{{ufs_current_cycle}}.coupler.res'

        # Use specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetForecastUFS.py")

        self.exec_cmd = ''   # Run on login node for S3 and R2D2 Database access
        self.include_header = ''
        self.login_node_limit = 'True'
