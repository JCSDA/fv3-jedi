# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.saveForecast as generic


class saveForecastGFS(generic.saveForecast):

    def setup(self, config, fc):

        # Get generic defaults
        generic.saveForecast.setup(self, config, fc)

        self.walltime = '00:10:00'

        # Use GFS specific script
        self.command = os.path.join(config['model_path'], "tasks/runSaveForecast.py")
