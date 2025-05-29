import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import matplotlib
import matplotlib.image as mplimg
import argparse

#initialize the argument parser and style sheet
parser = argparse.ArgumentParser()
parser.add_argument("--outFile", "-o" ,type=str,action="store",help="output file")
parser.add_argument("--psl1", "-f1" ,type=str,action='store',help="psl1 file")
parser.add_argument("--psl2", "-f2" ,type=str,action='store',help="pls2 file")
parser.add_argument("--gtf", "-g" ,type=str,action='store',help="gtf file")
parser.add_argument("--chr", "-c" ,type=str,action='store',help="chromosome region")
args = parser.parse_args()
outFile=args.outFile
psl1File=args.psl1
psl2File=args.psl2
gtfFile=args.gtf
chrRegion=args.chr

#Delete beofre sending#############################################################
from matplotlib import rcParams                                                   #
rcParams['font.family'] = 'Arial'                                                 #
rcParams['font.sans-serif'] = ['Arial']                                           #
matplotlib.font_manager.fontManager.addfont("/hb/home/bwong44/bme163/ARIAL.TTF")  #
###################################################################################

plt.style.use("BME163")
#Set up figure
figureWidth=5
figureHeight=6
plt.figure(figsize=(figureWidth,figureHeight))

# Parse the chromosome region flag
chrom, coords = chrRegion.split(":")
areaStart, areaEnd = map(int, coords.split("-"))

#colors
iOrange=(230/255,87/255,43/255)
iBlue=(88/255,85/255,120/255)
grey="Grey"

def makePanel(xpos,ypos,panelWidth, panelHeight, figureWidth, figureHeight):
    """Creates a panel in the figure with specified parameters"""
    relativePanelWidthCenter=panelWidth/figureWidth
    relativePanelHeightCenter=panelHeight/figureHeight
    panel=plt.axes([xpos/figureWidth,ypos/figureHeight,relativePanelWidthCenter,relativePanelHeightCenter],frameon=True)
    panel.set_xlim(0, 8)
    panel.set_ylim(0, 1262)
    panel.set_xticks([])
    panel.set_yticks([])
    return panel

def parsePSL(inputFile): 
    """Parse the input psl files and return a list of lists, reads and their metadata."""
    reads=[]
    with open(inputFile) as f:
        for line in f:
            stripLine=line.strip().split("\t")
            chr, start, end=stripLine[13], int(stripLine[15]), int(stripLine[16])
            blockwidths=np.array(stripLine[18].split(",")[:-1], dtype=int)
            blockstarts=np.array(stripLine[20].split(",")[:-1], dtype=int)
            keep=False
            if chr==chrom:
                if areaStart<start<areaEnd:
                    keep=True
                if areaStart<end<areaEnd:
                    keep=True
                if start<areaStart and end>areaEnd:
                    keep=True
            if keep:
                reads.append([chr,start,end,blockwidths,blockstarts])
    return reads

def parseGTF(gtfFile):
    """Parse the GTF file to extract gene information, including feature type."""
    transcripts = {}
    with open(gtfFile) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chromosome = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in ("exon", "CDS"):
                continue
            transcript_id = fields[8].split('transcript_id "')[1].split('"')[0]
            if transcript_id not in transcripts:
                transcripts[transcript_id] = []
            transcripts[transcript_id].append((chromosome, start, end, feature))

    reads = []
    for transcript_id, features in transcripts.items():
        chromosome = features[0][0]
        starts = []
        ends = []
        blockstarts = []
        blockwidths = []
        blockfeatures = []
        for chromosome, featureStart, featureEnd, feature in features:
            starts.append(featureStart)
            ends.append(featureEnd)
            blockstarts.append(featureStart)
            blockwidths.append(featureEnd - featureStart)
            blockfeatures.append(feature)  # Store the feature type for each block
        start = min(starts)
        end = max(ends)
        keep = False
        if chromosome == chrom:
            if areaStart < start < areaEnd:
                keep = True
            if areaStart < end < areaEnd:
                keep = True
            if start < areaStart and end > areaEnd:
                keep = True
        if keep:
            # Add blockfeatures to the read entry
            reads.append([chromosome, start, end, blockwidths, blockstarts, blockfeatures])
    return reads

def sortReads(reads,sortBy):
    """Sort the reads by either start or end position."""
    if sortBy == "start":
        reads.sort(key=lambda x: x[1])
    elif sortBy == "end":
        reads.sort(key=lambda x: x[2])
    return reads

def stackReads(reads):
    """Stack the reads efficiently, modifed from copilot"""
    yPosList=[None] * len(reads) #stores a list to track the given y value for each read
    assigned=[False] * len(reads) #stores a list to track if a read has been assigned a y value

    for y in range(1,len(reads)): #iterates through y values starting from 1
        lastEnd=-1 #tracks the end of last read, starts at -1 to allow first read to be placed
        for i in range(len(reads)):
            if not assigned[i]: #checks assignment list at index i
                start=reads[i][1]
                end=reads[i][2]
                if start > lastEnd:
                    yPosList[i]=y #stores the y values for each read
                    assigned[i]=True
                    lastEnd = end
        if all(assigned):
            break
    #add the assigned y value to each read
    stackedReads=[reads[i]+[yPosList[i]] for i in range(len(reads))]
    return stackedReads

def plotReads(reads, panel, color, edgewidth=0.2):
    """Plot the reads on the panels. Handles both PSL and GTF reads. Used Copilot to modify to handle GTF."""
    ymax = 0
    ymin = 0
    for read in reads:
        # Unpack read, handling both PSL and GTF (with or without blockfeatures)
        if len(read) == 7:
            chr, start, end, blockwidths, blockstarts, blockfeatures, y = read
        else:
            chr, start, end, blockwidths, blockstarts, y = read
            blockfeatures = None

        ymax = max(ymax, y)
        ymin = min(ymin, y)
        # Draw the full read line (for PSL, or for GTF if you want)
        rectangle = mplpatches.Rectangle(
            (start, y-0.025), end-start, 0.05,
            facecolor=color, linewidth=edgewidth, edgecolor="black", fill=True
        )
        panel.add_patch(rectangle)
        # Draw blocks
        for i in range(len(blockstarts)):
            blockstart = blockstarts[i]
            blockwidth = blockwidths[i]
            if blockfeatures:
                feature = blockfeatures[i]
                height = 0.5 if feature == "CDS" else 0.25
            else:
                height = 0.5  # Default for PSL
            rectangle = mplpatches.Rectangle(
                (blockstart, y - height/2), blockwidth, height,
                facecolor=color, linewidth=edgewidth, edgecolor="black", fill=True
            )
            panel.add_patch(rectangle)
    panel.set_xlim(areaStart, areaEnd)
    panel.set_ylim(ymin, ymax * 1.1)

def main():
    """Main function to run the script"""
    panel0=makePanel(0.1, 0.1, 4, 0.4, figureWidth, figureHeight) #bottom panel not used for histogram
    panel1=makePanel(0.1, 0.5, 4, 1.5, figureWidth, figureHeight) #bottom panel iBlue
    panel2=makePanel(0.1, 2.2, 4, 1.5, figureWidth, figureHeight) #middle panel iOrange
    panel3=makePanel(0.1, 3.9, 4, 1.5, figureWidth, figureHeight) #top panel grey

    readPSL2=parsePSL(psl2File)
    readPSL1=parsePSL(psl1File)
    readGTF=parseGTF(gtfFile)

    sortedPSL2=sortReads(readPSL2, "start")
    sortedPSL1=sortReads(readPSL1, "end")
    sortedGTF=sortReads(readGTF, "start")

    stackedPSL2Reads=stackReads(sortedPSL2)
    stackedPSL1Reads=stackReads(sortedPSL1)
    stackedGTFReads=stackReads(sortedGTF)
    
    plotReads(stackedPSL2Reads, panel1, iBlue, edgewidth=0)
    plotReads(stackedPSL1Reads, panel2, iOrange, edgewidth=0)
    plotReads(stackedGTFReads, panel3, grey, edgewidth=0.2)
main()

# save figure
plt.savefig(outFile, dpi=2400)

#python Wong_Brandon_BME163_Assignment_Final.py -f1 BME163_Input_Data_Final_1.psl -f2 BME163_Input_Data_Final_2.psl -g gencode.vM12.annotation.gtf -c chr7:45232000-45241000 -o test.png