# (C) Copyright 2020-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import ewok.tasks.getExpInit as generic


class getExpInitGFS(generic.getExpInit):

    def setup(self, config):

        # Get generic defaults
        generic.getExpInit.setup(self, config)

        # Use specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetAnalysis.py")
