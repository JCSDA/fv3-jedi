# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import ewok.tasks.getExpInit as generic


class getExpInitGFS(generic.getExpInit):

    def setup(self, config, fix):

        # Get generic defaults
        generic.getExpInit.setup(self, config, fix)

        # Use specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetAnalysis.py")

        self.exec_cmd = ''   # Run on login node for S3 and R2D2 Database access
        self.include_header = ''
        self.login_node_limit = 'True'
