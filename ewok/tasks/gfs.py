# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import ewok.tasks.GenericModel
import createPlotsFields
import getBackgroundGFS
import getExpInitGFS
import getFcInitGFS
import getFixFilesGFS
import getInitialConditionsGFS
import saveAnalysisGFS
import saveForecastGFS


class ModelTasks(ewok.tasks.GenericModel.ModelTasks):

    def __init__(self):
        ewok.tasks.GenericModel.ModelTasks.__init__(self)

        self.getBackground = getBackgroundGFS.getBackgroundGFS
        self.getExpInit = getExpInitGFS.getExpInitGFS
        self.getFcInit = getFcInitGFS.getFcInitGFS
        self.getStaticModel = getFixFilesGFS.getFixFilesGFS
        self.getInitialConditions = getInitialConditionsGFS.getInitialConditionsGFS
        self.saveAnalysis = saveAnalysisGFS.saveAnalysisGFS
        self.saveForecast = saveForecastGFS.saveForecastGFS
        self.createPlots = ewok.createPlots
        self.createPlotsFields = createPlotsFields.createPlotsFields
        self.savePlots = ewok.savePlots
        self.savePlotsFields = ewok.savePlotsFields
