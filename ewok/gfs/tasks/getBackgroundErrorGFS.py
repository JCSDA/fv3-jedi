# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.getBackgroundError as generic


class getBackgroundErrorGFS(generic.getBackgroundError):

    def setup(self, config, fc):

        # Get generic defaults
        generic.getBackgroundError.setup(self, config, fc)

        if 'hack_step_ens' in config and config['hack_step_ens'] == True:
            self.RUNTIME_YAML['hack_step_bg'] = True

        # Use GFS specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetForecast.py")

        self.exec_cmd = ''   # Run on login node for S3 and R2D2 Database access
        self.include_header = ''
        self.login_node_limit = 'True'
