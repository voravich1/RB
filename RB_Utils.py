import sys

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
    prior_pos = ""

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

                if idx > 8:
                    dummyFormatHead = format.split(":")
                    dummyFormatInfo = array_vcfInfo[idx].split()
                    dummySample = array_Column_Header[idx]

                    for colNum in list(range(len(dummyFormatHead))):
                        key = dummyFormatHead[colNum]
                        value = dummyFormatInfo[colNum]
                        format_dict.update({key:value})

                    sample_format_dict.update({dummySample:format_dict})
                    continue

            # Start part analyze data for rain fall plot


