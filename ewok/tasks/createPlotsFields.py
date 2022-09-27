# (C) Copyright 2020-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
from ewok import Task
import yamltools


class createPlotsFields(Task):

    def setup(self, config, statefile):
        fieldfile = statefile['an']

        plotsconf = yamltools.parse_config(os.path.join(config['model_path'], 'defaults/plotFields.yaml'))

        self.RUNTIME_YAML['levels'] = plotsconf['levels']
        self.RUNTIME_YAML['filepath'] = os.path.join(fieldfile['datapath'],
                                                    fieldfile['filename_core'])
        self.RUNTIME_YAML['gridfiledir'] = os.path.join(config['model_path'],
                                                    'tasks/plots/fv3grid')

        imgname = config['expid']

        self.RUNTIME_YAML['output'] = os.path.join(fieldfile['datapath'], imgname)
        self.RUNTIME_YAML['variables'] = plotsconf['variables']

        self.command = os.path.join(config['model_path'],
                                "tasks/plots/plot_gfs.py")


        self.walltime = '00:20:00'

        self.output['datapath'] = fieldfile['datapath']
        self.output['levels'] = plotsconf['levels']
        self.output['variables'] = plotsconf['variables']
