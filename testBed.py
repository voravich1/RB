#!/Users/worawich/miniconda3/envs/Thalassemia_CDSS/bin/python

# In order to use matplotlib on macOs on conda env
# You must conda install python.app
# then change the interpreter from python to pythonw
# After that use can import any matplotlib library as common

import sys
import RB_Utils as rb


variant_pos_list_dict, mutationType_list_dict = rb.readVCF("/Users/worawich/Download_dataset/Ratina_cancer/Mutect_dnabrick_result/set1_dnabrick/604_somatic.vcf")

rb.plotRainFallMutationType(variant_pos_list_dict, mutationType_list_dict)
