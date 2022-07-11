# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
from ewok import Task


class getStaticB(Task):

    def setup(self, config, build):

        if ('_balfile' in config['BACKGROUND_ERROR']):
            balfile = config['BACKGROUND_ERROR']['_balfile']
            bump2d = config['BACKGROUND_ERROR']['_prefix2d']
            bump3d = config['BACKGROUND_ERROR']['_prefix3d']
    
            fv3_data_dir = os.path.join(build['builddir'], 'fv3-jedi/test/Data')
            static_b_dir = os.path.join(self.workdir['wdir'], 'staticB')
    
            self.output['STATICB'] = {}
            self.output['STATICB']['datadir'] = static_b_dir
    
            self.RUNTIME_ENV['FV3DATA'] = fv3_data_dir
            self.RUNTIME_ENV['STATBDIR'] = static_b_dir
            self.RUNTIME_ENV['BALFILE'] = balfile
            self.RUNTIME_ENV['PREFIX2D'] = bump2d
            self.RUNTIME_ENV['PREFIX3D'] = bump3d
            self.command = os.path.join(config['model_path'], "tasks/getstaticb.sh")
        else:
            self.notask = True
