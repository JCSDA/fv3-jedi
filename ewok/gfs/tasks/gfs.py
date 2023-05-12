# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import ewok.tasks.GenericModel
import getBackgroundGFS
import getBackgroundErrorGFS
import getExpInitGFS
import getFcInitGFS
import getFixFilesGFS
import getInitialConditionsGFS
import saveAnalysisGFS
import saveForecastGFS
import saveForecastGFS_DB


class ModelTasks(ewok.tasks.GenericModel.ModelTasks):

    def __init__(self):
        ewok.tasks.GenericModel.ModelTasks.__init__(self)

        self.convertstate = ewok.convertstate
        self.enstats = ewok.enstats
        self.getBackground = getBackgroundGFS.getBackgroundGFS
        self.getBackgroundError = getBackgroundErrorGFS.getBackgroundErrorGFS
        self.getExpInit = getExpInitGFS.getExpInitGFS
        self.getFcInit = getFcInitGFS.getFcInitGFS
        self.getStaticModel = getFixFilesGFS.getFixFilesGFS
        self.getInitialConditions = getInitialConditionsGFS.getInitialConditionsGFS
        self.plotAnalysis = ewok.plotAnalysis
        self.plotEnsStats = ewok.plotEnsStats
        self.plotObsStats = ewok.plotObsStats
        self.publishAnPlots = ewok.publishAnPlots
        self.publishEnsPlots = ewok.publishEnsPlots
        self.publishObsPlots = ewok.publishObsPlots
        self.saveAnalysis = saveAnalysisGFS.saveAnalysisGFS
        self.saveForecast = saveForecastGFS.saveForecastGFS
        self.saveForecastDB = saveForecastGFS_DB.saveForecastGFS_DB
