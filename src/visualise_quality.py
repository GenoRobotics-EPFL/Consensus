"""
This script visualises the distribution of quality of bases in reads of a single fastq file.
The fastq file comes from the sequensor, so it has multiple reads (though the script also works with one read only).
The output is a png image in the current directory with the name of the fastq file provided as input.

(works on linux and windows)

usage:

python3 visualise_quality.py path/To/FastQFile
"""

import sys
import os
import matplotlib.pyplot as plt


def getContentFile(pathFastqFile: str):

    # open file and get all the reads
    allReadsAsString = []

    with open(pathFastqFile) as fastq:
        allReadsAsString = fastq.readlines()

    return allReadsAsString


def extractReadQuality(fileContent):

    # the fileContent must contain multiple reads, where each read is 4 lines:
    if not (len(fileContent) % 4 == 0):
        raise ValueError("fastq file must have reads composed of 4 lines")

    qualityOfReads = []

    # the quality of each read is always the 4th element
    for i in range(3, len(fileContent), 4):

        line = fileContent[i]

        # remove the newline character at the end
        line = line[:-1]

        qualityOfReads.append(line)

    return qualityOfReads


def convertQuality(allReadsQuality):

    # note that here we can't distinguish the reads from each other anymore
    # they are all in one list
    convertedQualities = []

    for readQuality in allReadsQuality:

        for rawQuality in readQuality:

            # transform the character into the int quality
            score = ord(rawQuality) - 33

            convertedQualities.append(score)

    return convertedQualities


def visualiseQualities(pathFastqFile, readQualityConverted):

    plt.figure(figsize=(10, 6))

    plt.hist(readQualityConverted, bins=range(0, 95))

    plt.ylabel("Frequency")

    plt.xlabel("Base quality score 0-93 (higher is better)")
    plt.xlim(0, 95)

    nameFile = os.path.basename(pathFastqFile)
    plt.title(nameFile)

    nameOutput = nameFile + ".png"
    plt.savefig(nameOutput)  # saves in current directory


def main():

    # parse arguments
    if len(sys.argv) != 2:
        raise ValueError("You provide exactly one fastq file")

    pathFastqFile = sys.argv[1]

    if not pathFastqFile.endswith(".fastq"):
        raise ValueError("must provide a fastq file")

    # open file and get content
    fileContent = getContentFile(pathFastqFile)

    # extract quality of reads from file
    allReadsQuality = extractReadQuality(fileContent)

    # convert the quality score to score between 0 and 64
    readQualityConverted = convertQuality(allReadsQuality)

    # visualise histogram
    visualiseQualities(pathFastqFile, readQualityConverted)


if __name__ == "__main__":
    main()
