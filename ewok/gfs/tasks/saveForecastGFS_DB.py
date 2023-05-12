# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.saveForecast as generic


class saveForecastGFS_DB(generic.saveForecast):

    def setup(self, config, fc):

        # Get generic defaults
        generic.saveForecast.setup(self, config, fc)

        # Get target values for `store`
        self.RUNTIME_YAML['resolution'] = config['GEOMETRY_OUT']['_resol_name']

        if 'target_DB' in config:
            self.RUNTIME_YAML['target_DB'] = config['target_DB']

        if 'target_expid' in config:
            self.RUNTIME_YAML['target_expid'] = config['target_expid']
        else:
            self.RUNTIME_YAML['target_expid'] = config['expid']

        self.walltime = '00:10:00'

        # Use GFS specific script
        self.command = os.path.join(config['model_path'], "tasks/runSaveForecastDB.py")

        self.exec_cmd = ''   # Run on login node for S3 and R2D2 Database access
        self.include_header = ''
        self.login_node_limit = 'True'
