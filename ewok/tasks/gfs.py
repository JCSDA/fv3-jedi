# (C) Copyright 2020-2021 UCAR
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import ewok.tasks.GenericModel
import getBackgroundGFS
import getFixFilesGFS
import getInitialConditionsGFS

class ModelTasks(ewok.tasks.GenericModel.GenericModelTasks):

    def __init__(self):
        ewok.tasks.GenericModel.GenericModelTasks.__init__(self)

        self.getBackground = getBackgroundGFS.getBackgroundGFS
        self.getFixFiles = getFixFilesGFS.getFixFilesGFS
        self.getInitialConditions = getInitialConditionsGFS.getInitialConditionsGFS

