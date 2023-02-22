# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.saveAnalysis as generic


class saveAnalysisGFSAero(generic.saveAnalysis):

    def setup(self, config, an):

        # Get generic defaults
        generic.saveAnalysis.setup(self, config, an)

        # Use GFS specific script
        self.command = os.path.join(config['model_path'], "tasks/runSaveAnalysisAero.py")
