# (C) Copyright 2020-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.getInitialConditions as generic


class getInitialConditionsGEOS(generic.getInitialConditions):

    def setup(self, config):

        # Get generic defaults
        generic.getInitialConditions.setup(self, config)

        self.walltime = '00:05:00'

        # Use specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetForecast.py")
