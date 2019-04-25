#!/Users/worawich/miniconda3/envs/Thalassemia_CDSS/bin/python

import sys
import RB_Utils as rb


variant_pos_list_dict, mutationType_list_dict = rb.readVCF("/Users/worawich/Download_dataset/Ratina_cancer/Mutect_dnabrick_result/set1_dnabrick/604_somatic.vcf")

rb.plotRainFallMutationType(variant_pos_list_dict, mutationType_list_dict)
