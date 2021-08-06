# (C) Copyright 2017-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os

# Loop over test log files
run_files = os.listdir()
for run_file in run_files:

    # Check if the file is a log file
    if run_file[-9:] == '.test.out':

        # Output file
        ref_file = run_file[0:-8]+'ref'
        out_file = open(ref_file, 'w')
        read_file = open(run_file, 'r')

        content = read_file.readlines()
        type(content)

        for ii in range(0, len(content)):
            if(content[ii] != "\n"):
                out_file.write(content[ii])
            else:
                pass

        read_file.close()
        out_file.close()
