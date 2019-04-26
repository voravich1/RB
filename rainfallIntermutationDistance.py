#!/Users/worawich/miniconda3/envs/Thalassemia_CDSS/bin/python

# In order to use matplotlib on macOs on conda env
# You must conda install python.app
# then change the interpreter from python to pythonw
# After that use can import any matplotlib library as common

import sys
import RB_Utils as rb

inputFile=sys.argv[1]
saveFilePath=sys.argv[2]
saveFilePrefix=sys.argv[3]



variant_pos_list_dict, mutationType_list_dict = rb.readVCF(inputFile)

rb.plotRainFallMutationType(variant_pos_list_dict, mutationType_list_dict, saveFilePath, saveFilePrefix)