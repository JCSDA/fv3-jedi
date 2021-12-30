# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
from ewok import Task


class getFixFilesGFS(Task):

    def setup(self, config, **inputs):

        localconf = {}
        localconf['fixdir'] = self.config['currentdir']
        localconf['resol'] = config['GEOMETRY']['_resol_name']
        localconf['nlevs'] = str(config['GEOMETRY']['npz'])

        tmplfile = os.path.join(config['model_path'], "templates/fixgeom.yaml")
        geomtmpl = yamltools.parse_config(tmplfile)
        geomconf = yamltools.substitute_template_variables(geomtmpl, localconf)
        self.output['FIXGEOM'] = geomconf

        tmplfile = os.path.join(config['model_path'], "templates/fixmodel.yaml")
        modeltmpl = yamltools.parse_config(tmplfile)
        modelconf = yamltools.substitute_template_variables(modeltmpl, localconf)
        self.output['FIXMODEL'] = modelconf

        self.config['ENV']['FV3REPO'] = os.path.join(os.environ.get("JEDI_SRC"), "fv3-jedi")
        self.config['ENV']['FIXDIR'] = self.config['currentdir']
        self.config['ENV']['RESOL'] = config['GEOMETRY']['_resol_name']
        self.command = os.path.join(config['model_path'], "tasks/getfixfiles.sh")
