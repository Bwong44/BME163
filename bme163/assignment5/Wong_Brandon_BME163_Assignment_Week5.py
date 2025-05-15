import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import matplotlib
import matplotlib.image as mplimg
import argparse

#initialize the argument parser and style sheet
parser = argparse.ArgumentParser()
parser.add_argument("--outFile", "-o" ,type=str,action="store",help="output file")
parser.add_argument("--spliceSequences", "-s" ,type=str,action='store',help="spliceSequences")
parser.add_argument("--Apng", "-A" ,type=str,action='store',help="A.png file")
parser.add_argument("--Tpng", "-T" ,type=str,action='store',help="T.png file")
parser.add_argument("--Gpng", "-G" ,type=str,action='store',help="G.png file")
parser.add_argument("--Cpng", "-C" ,type=str,action='store',help="C.png file")
args = parser.parse_args()
outFile=args.outFile
spliceSequences=args.spliceSequences
Apng=args.Apng
Tpng=args.Tpng
Gpng=args.Gpng
Cpng=args.Cpng

#Delete beofre sending#############################################################
from matplotlib import rcParams                                                   #
rcParams['font.family'] = 'Arial'                                                 #
rcParams['font.sans-serif'] = ['Arial']                                           #
matplotlib.font_manager.fontManager.addfont("/hb/home/bwong44/bme163/ARIAL.TTF")  #
###################################################################################
plt.style.use("BME163")

#Set up figure
figureWidth=5
figureHeight=2
plt.figure(figsize=(figureWidth,figureHeight))

def makePanel(xpos,ypos,panelWidth, panelHeight, figureWidth, figureHeight, ylabel="", hideTicks=False, title=""):
    """Creates a panel in the figure with specified parameters"""
    relativePanelWidthCenter=panelWidth/figureWidth
    relativePanelHeightCenter=panelHeight/figureHeight
    panel=plt.axes([xpos/figureWidth,ypos,relativePanelWidthCenter,relativePanelHeightCenter],frameon=True)
    panel.set_xlim(-10, 10)
    panel.set_ylim(0, 2)
    if hideTicks:
        panel.set_yticks([])
    panel.set_xlabel("Distance to\nSplice Site") #x-axis labels
    panel.set_ylabel(ylabel) #y-axis labels
    panel.set_title(title)
    panel.plot([0, 0], [0, 2], color="black", linestyle="-", linewidth=0.5)
    return panel

def fastaReader(file, search_char=None):
    """Reads a FASTA file and returns all sequences from copilot."""
    sequences = []

    with open(file, "r") as f:
        current_sequence = ""
        current_header = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if current_sequence:  # Save the previous sequence if it matches the search criteria
                    if search_char is None or search_char in current_header:
                        sequences.append(current_sequence)
                    current_sequence = ""
                current_header = line  # Update the current header
            else:
                current_sequence += line  # Append sequence lines

        # Add the last sequence if it matches the search criteria
        if current_sequence and (search_char is None or search_char in current_header):
            sequences.append(current_sequence)

    return sequences

def calculateLogo(sequences):
    """Calculates the heights of nucleotides for a sequence logo."""
    totalSeq=len(sequences)
    posList=[]
    freqList=[]
    icList=[]
    heightDict={}

    #loop to create a list of dictionaries for each position
    for i in range(20): 
        nucDict={"A":0,"T":0,"G":0,"C":0}
        posList.append(nucDict)

    #loop through each sequence and each nuc to find nuc counts
    for seq in sequences: 
        for index, nuc in enumerate(seq):
            if nuc in posList[index]:
                posList[index][nuc]+=1

    #loops to find freq of each nucleotide at each position
    for pos in posList:
        freq={nuc:count/totalSeq for nuc, count in pos.items()}
        freqList.append(freq)     

    #calculates small-sample correction
    en=(1/np.log(2))*((4-1)/(2*totalSeq)) 
    
    #calcualtes ic for each position
    for i,freq in enumerate(freqList): 
        entropy=0
        for nuc in freq:
            if freq[nuc]>0: #avoids log 0 error
                entropy-=freq[nuc]*np.log2(freq[nuc])
        ic=2-(entropy+en) #information content
        icList.append(ic)

    #calcualtes height for each nucleotide at each position
    for index, (ic, freq) in enumerate(zip(icList, freqList)):
        heightDict[index]={}
        for nucleotide, frequency in freq.items():
            height=frequency*ic
            heightDict[index][nucleotide] = height

    return heightDict

def plotPanel(panel, heightDict):
    """Plots the sequence logo using heightDict, a lot from copilot."""
    # Load nucleotide images
    A = mplimg.imread(Apng)
    T = mplimg.imread(Tpng)
    G = mplimg.imread(Gpng)
    C = mplimg.imread(Cpng)

    nucleotideImages = {"A": A, "T": T, "G": G, "C": C} #Map nucleotides to their respective images
 
    #Loop through each position in heightDict
    for position, heights in heightDict.items():
        xStart = position - 10  #Map position index to x-axis (-10 to 10)
        yStart = 0  #Start stacking from y=0

        # Sort nucleotides by height (most frequent on top)
        sortedBases=sorted(heights.items(), key=lambda x: x[1])

        # Plot each base at the current position
        for nucleotide, height in sortedBases:
            if height>0:  # Only plot if height is greater than 0
                panel.imshow(
                    nucleotideImages[nucleotide],
                    extent=(xStart, xStart+1, yStart, yStart+height),
                    aspect="auto",
                    origin="upper"
                )
                yStart += height  # Stack the next base on top

def main():
    """Main function to run the script"""
    panel1=makePanel(0.5, 0.3, 1.5, 0.5, 5, 2, "Bits", title="5'SS")
    panel2=makePanel(2.2, 0.3, 1.5, 0.5, 5, 2, hideTicks=True, title="3'SS")
    bin5=fastaReader(spliceSequences,">5")
    bin3=fastaReader(spliceSequences,">3")
    heightDict5=calculateLogo(bin5)
    heightDict3=calculateLogo(bin3)
    plotPanel(panel1, heightDict5)
    plotPanel(panel2, heightDict3)
main()

# save figure
plt.savefig(outFile, dpi=600)