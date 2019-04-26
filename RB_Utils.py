import sys
import os
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf
import math

def readVCF(inVCF):

    array_Column_Header = ()
    array_vcfInfo = []
    chr = ""
    pos = ""
    id = ""
    ref = ""
    alt = ""
    qual = ""
    filter = ""
    info = ""
    info_dict = {}
    format = ""
    format_dict = {}
    sample_format_dict = {}

    variantNumber = 0
    mutation_type = ""
    prior_chr = ""
    prior_pos = 0
    first_line_data_flag = True


    # variable require for create rain fall plot
    variant_pos_list = []        # list of SNP position
    variant_pos_list_dict = {}   # key is chromosome | value is list of SNP position
    mutationType_list = []       # list of mutation type
    mutationType_list_dict = {}  # key is chromosome | value is list of mutation type

    with open(inVCF) as vcfFile:

        for line in vcfFile:

            if line[0:2] == "##":  #skip header
                continue
            elif line[0:2] == "#C":
                array_Column_Header = line.split()
                continue

            array_vcfInfo = line.split()

            chr = array_vcfInfo[0]
            pos = array_vcfInfo[1]
            id = array_vcfInfo[2]
            ref = array_vcfInfo[3]
            alt = array_vcfInfo[4]
            qual = array_vcfInfo[5]
            filter = array_vcfInfo[6]
            info = array_vcfInfo[7]
            format = array_vcfInfo[8]

            for idx in list(range(len(array_Column_Header))):
                format_dict = {}
                if idx > 8:
                    dummyFormatHead = format.split(":")
                    dummyFormatInfo = array_vcfInfo[idx].split(":")
                    dummySample = array_Column_Header[idx]

                    for colNum in list(range(len(dummyFormatHead))):
                        key = dummyFormatHead[colNum]
                        if colNum >= len(dummyFormatInfo):
                            value = "."
                        else:
                            value = dummyFormatInfo[colNum]
                        format_dict.update({key: value})

                    sample_format_dict.update({dummySample: format_dict})
                    continue

            # Start part analyze data for rain fall plot
            #print("555")


            mutation_type = classifyMutationType(ref , alt)

            if mutation_type == "null":
                continue    # Skip other variant that is not SNP

            if(first_line_data_flag == True):
                prior_pos = pos
                prior_chr = chr
                first_line_data_flag = False

            if prior_chr == chr:
                variant_pos_list.append(pos)
                mutationType_list.append(mutation_type)
            elif prior_chr != chr:
                variant_pos_list_dict.update({prior_chr: variant_pos_list})
                mutationType_list_dict.update({prior_chr: mutationType_list})
                variant_pos_list = []
                mutationType_list = []
                variant_pos_list.append(pos)
                mutationType_list.append(mutation_type)
                prior_chr = chr
                prior_pos = pos

        variant_pos_list_dict.update({prior_chr: variant_pos_list})
        mutationType_list_dict.update({prior_chr: mutationType_list})

    return variant_pos_list_dict, mutationType_list_dict



def classifyMutationType(ref, alt):

    mutationType = str(ref + alt)
    mutationType.upper()

    if len(mutationType) > 2:
        return "null"

    if mutationType == "CA" or mutationType == "GT":
        return "CA"

    if mutationType == "CG" or mutationType == "GC":
        return "CG"

    if mutationType == "CT" or mutationType == "GA":
        return "CT"

    if mutationType == "TA" or mutationType == "AT":
        return "TA"

    if mutationType == "TC" or mutationType == "AG":
        return "TC"

    if mutationType == "TG" or mutationType == "AC":
        return "TG"

def plotRainFallMutationType(variant_pos_list_dict: dict, mutationType_list_dict: dict, savePath, saveFilePrefix):
    saveFile = os.path.join(savePath, saveFilePrefix + '_plotTestest.pdf')
    pdf = matplotlib.backends.backend_pdf.PdfPages(saveFile)
    column_max = 4
    row_max = 6
    count = 0
    for chrName in variant_pos_list_dict.keys():
        if chrName == "chrM" or chrName == "MT":
            continue
        count+=1
        #plt.figure()
        fig, ax = plt.subplots(figsize=(20,10))

        mutation_list = mutationType_list_dict.get(chrName)
        variant_list = variant_pos_list_dict.get(chrName)

        n = len(mutation_list)

        axisPlot = []
        listCA_x = []
        listCA_y = []
        listCG_x = []
        listCG_y = []
        listCT_x = []
        listCT_y = []
        listTA_x = []
        listTA_y = []
        listTC_x = []
        listTC_y = []
        listTG_x = []
        listTG_y = []

        for idx in range(1,len(variant_list)):
            query_idx = idx-1

            y = math.log10(int(variant_list[query_idx+1]) - int(variant_list[query_idx]))
            x = int(variant_list[query_idx+1])
            mutationType = mutation_list[query_idx+1]

            if mutationType == "CA":
                listCA_x.append(x)
                listCA_y.append(y)
            elif mutationType == "CG":
                listCG_x.append(x)
                listCG_y.append(y)
            elif mutationType == "CT":
                listCT_x.append(x)
                listCT_y.append(y)
            elif mutationType == "TA":
                listTA_x.append(x)
                listTA_y.append(y)
            elif mutationType == "TC":
                listTC_x.append(x)
                listTC_y.append(y)
            elif mutationType == "TG":
                listTG_x.append(x)
                listTG_y.append(y)

        for idx in range(6):

            if idx == 0:
                color = "blue"
                mutationType="C>A"
                ax.scatter(listCA_x, listCA_y, c=color, label=mutationType, alpha=0.75, edgecolors='none')
            elif idx == 1:
                color = "black"
                mutationType = "C>G"
                ax.scatter(listCG_x, listCG_y, c=color, label=mutationType, alpha=0.75, edgecolors='none')
            elif idx == 2:
                color = "red"
                mutationType = "C>T"
                ax.scatter(listCT_x, listCT_y, c=color, label=mutationType, alpha=0.75, edgecolors='none')
            elif idx == 3:
                color = "pink"
                mutationType = "T>A"
                ax.scatter(listTA_x, listTA_y, c=color, label=mutationType, alpha=0.75, edgecolors='none')
            elif idx == 4:
                color = "yellow"
                mutationType = "T>C"
                ax.scatter(listTC_x, listTC_y, c=color, label=mutationType, alpha=0.75, edgecolors='none')
            elif idx == 5:
                color = "green"
                mutationType = "T>G"
                ax.scatter(listTG_x, listTG_y, c=color, label=mutationType, alpha=0.75, edgecolors='none')


        saveFile = os.path.join(savePath,saveFilePrefix+'_plot_'+chrName+'.pdf')
        ax.legend()
        ax.grid(True)

        #plt.subplot(row_max, column_max, count)
        #plt.figure(count)
        plt.title("Rainfall plot of "+chrName)
        plt.xlabel("Genomic position")
        plt.ylabel("Intermutation distance (log(bp))")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        #plt.savefig(saveFile, bbox_inches='tight')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    pdf.close()
    #plt.show()

#def plotMultipleRainFallMutationType()