# (C) Copyright 2020-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.getBackground as generic


class getBackgroundGEOS(generic.getBackground):

    def setup(self, config, fc):

        # Get generic defaults
        generic.getBackground.setup(self, config, fc)

        # Use GEOS specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetForecast.py")

        self.exec_cmd = ''   # Run on login node for S3 and R2D2 Database access
        self.include_header = ''
        self.login_node_limit = 'True'
