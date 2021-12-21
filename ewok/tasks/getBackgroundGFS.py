# (C) Copyright 2020-2021 UCAR
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.getBackground as generic

class getBackgroundGFS(generic.getBackground):

    def setup(self, config, fc):

        # Get generic defaults
        generic.getBackground.setup(self, config, fc)

        # Use GFS specific script
        self.command = os.path.join(config['model_path'],"tasks/getForecastGFS.py")

