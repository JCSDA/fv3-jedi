# (C) Copyright 2023 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
from ewok import Task
class getFixFilesUFS(Task):
    def setup(self, config, **inputs):
        localconf = {}
        localconf['fixdir'] = self.workdir['wdir']
        localconf['resol'] = config['GEOMETRY']['_resol_name']
        localconf['nlevs'] = str(config['GEOMETRY']['npz'])
        ilayout = config['GEOMETRY']['layout']
        layout = str(ilayout[0]) + "x" + str(ilayout[1])

        tmplfile = os.path.join(config['model_path'], "templates/fixgeomUFS.yaml")
        geomtmpl = yamltools.parse_config(tmplfile)
        geomconf = yamltools.substitute_template_variables(geomtmpl, localconf)
        self.output['FIXGEOM'] = geomconf

        tmplfile = os.path.join(config['model_path'], "templates/fixmodelUFS.yaml")
        modeltmpl = yamltools.parse_config(tmplfile)
        modelconf = yamltools.substitute_template_variables(modeltmpl, localconf)
        self.output['FIXMODEL'] = modelconf

        list_filename = 'defaults/'+config['MODEL']['_list_static_files']
        tmplfile = os.path.join(config['model_path'], list_filename)
        listfiles = yamltools.parse_config(tmplfile)

#       Fill RUNTIME_YAML conf with necessary paths
        self.RUNTIME_YAML['fv3repo'] = os.path.join(os.environ.get("JEDI_SRC"), "fv3-jedi")
        self.RUNTIME_YAML['listfiles'] = listfiles
        self.RUNTIME_YAML['static_data'] = config['static_data']

        self.RUNTIME_YAML['fixdir'] = self.workdir['wdir']
        self.RUNTIME_YAML['nlevs'] = str(config['GEOMETRY']['npz'])
        self.RUNTIME_YAML['layout'] = layout
        self.RUNTIME_YAML['resol'] = config['GEOMETRY']['_resol_name']


        self.output['ufs_modeldir'] = os.path.join(self.workdir['wdir'], 'UFSFixFiles')

        self.command = os.path.join(config['model_path'], "tasks/runGetFixFilesUFS.py")
