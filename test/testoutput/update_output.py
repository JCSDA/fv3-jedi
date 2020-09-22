# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import re

test_string = 'Test     :'

# Loop over log files
runfiles = os.listdir()
for runfile in runfiles:

    # Check if the file is a log file
    if runfile[-4:] == '.run':

        # Output file
        runreffile = runfile[0:-3]+'ref'
        outfile = open(runreffile, "w")

        # Open file
        for line in open(runfile, 'r'):

            if re.search(test_string, line):
                outfile.write(line)
                if line == None:
                    print('no matches found')

        outfile.close()

            # Grep for test string
            #print(line)
