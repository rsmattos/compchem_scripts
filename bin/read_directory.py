#! /usr/bin/env python3

import os
import fnmatch

###################      READING FILES     ###################
# searches in the folder or passed arguments and creates a list
# of paths to the calculation outputs
def find_outputs(args):
    paths=[]
    # output_file=''
    print("The output files being used are:")

    for inp in args.input:
        if os.path.isfile(inp):
            paths.append(inp)

        elif os.path.isdir(inp):
            # if the path given is a directory, searches inside for output files
            for path,dir,file in os.walk(inp):
                for file_name in file:
                    if fnmatch.fnmatch(file_name, '*'+args.extension):
                        print(path+'/'+file_name)
                        paths.append(path+'/'+file_name)
                        # output_file = file_name

        else:
            # output_file = inp
            for path,dir,file in os.walk('.'):
                for file_name in file:
                    if fnmatch.fnmatch(file_name, inp):
                        print(path+'/'+file_name)
                        paths.append(path+'/'+file_name)

    return paths
