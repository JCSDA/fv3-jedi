# (C) Copyright 2020-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import ewok.tasks.GenericModel
import getBackgroundGEOS
import getExpInitGEOS
import getFcInitGEOS
import getFixFilesGEOS
import getInitialConditionsGEOS
import getStaticB
import saveAnalysisGEOS
import saveForecastGEOS


class ModelTasks(ewok.tasks.GenericModel.ModelTasks):

    def __init__(self):
        ewok.tasks.GenericModel.ModelTasks.__init__(self)

        self.getBackground = getBackgroundGEOS.getBackgroundGEOS
        self.getExpInit = getExpInitGEOS.getExpInitGEOS
        self.getFcInit = getFcInitGEOS.getFcInitGEOS
        self.getStaticB = getStaticB.getStaticB
        self.getStaticModel = getFixFilesGEOS.getFixFilesGEOS
        self.getInitialConditions = getInitialConditionsGEOS.getInitialConditionsGEOS
        self.saveAnalysis = saveAnalysisGEOS.saveAnalysisGEOS
        self.saveForecast = saveForecastGEOS.saveForecastGEOS
        self.createPlots = ewok.createPlots
        self.savePlots = ewok.savePlots
