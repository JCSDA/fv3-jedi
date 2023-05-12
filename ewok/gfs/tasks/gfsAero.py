# (C) Copyright 2020-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import getFixFilesGFSAero
import saveAnalysisGFSAero
import saveForecastGFSAero
import saveForecastGFS_DB

import gfs

class ModelTasks(gfs.ModelTasks):

    def __init__(self):
        gfs.ModelTasks.__init__(self)

        self.getStaticModel = getFixFilesGFSAero.getFixFilesGFSAero
        self.saveAnalysis = saveAnalysisGFSAero.saveAnalysisGFSAero
        self.saveForecast = saveForecastGFSAero.saveForecastGFSAero
        self.saveForecastDB = saveForecastGFS_DB.saveForecastGFS_DB
