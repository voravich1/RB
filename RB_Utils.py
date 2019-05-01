import ntpath
import sys
import os
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf
import math
import numpy as np

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

def plotHeatMapTrinucleotide(Trinucleotide_raw_File):
    x_3prime = []
    y_5prime = []
    freq_data = []
    fileName = ntpath.basename(Trinucleotide_raw_File)
    sampleName = fileName.split("_")[0]
    firstLineFlag = True
    firstSecondLineFlag = True
    data_row_count = 0

    mutation_type_dict = {}

    with open(Trinucleotide_raw_File) as rawDataFile:
        for line in rawDataFile:
            line = line.strip('\n')
            if firstLineFlag==True:
                firstLineFlag=False
                continue
            elif line[0] == "M":
                if firstSecondLineFlag == False:
                    mutation_type_dict.update({mutation_type: freq_data_array})
                info = line.split()
                mutation_type=info[3]
                firstSecondLineFlag = False
            elif line[0] == "B":
                info = line.split(",")
                freq_data = []
                for idx in range(1,len(info)):
                    x_3prime.append(str(info[idx]))
            else:
                info = line.split(",")
                for idx in range(len(info)):
                    if idx == 0:
                        y_5prime.append(str(info[idx]))
                        value_list = []
                    else:
                        value_list.append(int(info[idx]))
                freq_data.append(value_list)
                freq_data_array = np.asarray(freq_data)
        mutation_type_dict.update({mutation_type: freq_data_array})


    # Continue Here
    # start plot heatMap
    fig, ax = plt.subplots()

    for key, value in mutation_type_dict.items():
        mutation_type = key
        freq_data_array = value


        im, cbar = heatmap(freq_data_array,y_5prime,x_3prime,ax=ax,cmap="YlGn",cbarlabel="InterMutation (log10)")

        fig.tight_layout()
        plt.show()

    print()






def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                     This should either use the string format method, e.g.
                     "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                     used for values below a threshold, the second for those
                     above.
        threshold  : Value in data units according to which the colors from
                     textcolors are applied. If None (the default) uses the
                     middle of the colormap as separation.

    Further arguments are passed on to the created text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[im.norm(data[i, j]) > threshold])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts