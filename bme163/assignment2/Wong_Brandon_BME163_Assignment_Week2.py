import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import matplotlib
import argparse

#initialize the argument parser and style sheet
parser = argparse.ArgumentParser()
parser.add_argument("--outFile", "-o" ,type=str,action="store",help="output file")
parser.add_argument("--inFile", "-i" ,type=str,action='store',help="input file")
args = parser.parse_args()
outFile=args.outFile
input=args.inFile
plt.style.use("BME163")

#colors
iBlue=(88/255,85/255,120/255)
Grey="grey"
iGreen=(120/255,172/255,145/255)

#Set up figure
figureWidth=3
figureHeight=3
plt.figure(figsize=(figureWidth,figureHeight))

#set up panels

#main panel
panelLeftMain=0.7
panelBottomMain=0.3
panelWidthMain=1.5
panelHeightMain=1.5
panelMain=plt.axes([panelLeftMain/figureWidth,panelBottomMain/figureHeight,panelWidthMain/figureWidth,panelHeightMain/figureHeight])
panelMain.set_xlim(0,15)
panelMain.set_ylim(0,15)
panelMain.set_xticks([0,15,5,10])
panelMain.set_yticks([])

#top panel
panelLeftTop=0.7
panelBottomTop=1.87
panelWidthTop=1.5
panelHeightTop=0.25
panelTop=plt.axes([panelLeftTop/figureWidth,panelBottomTop/figureHeight,panelWidthTop/figureWidth,panelHeightTop/figureHeight])
panelTop.set_xlim(0,15)
panelTop.set_ylim(0,20)
panelTop.set_xticks([])

#left panel
panelLeftLeft=0.38
panelBottomLeft=0.3
panelWidthLeft=0.25
panelHeightLeft=1.5
panelLeft=plt.axes([panelLeftLeft/figureWidth,panelBottomLeft/figureHeight,panelWidthLeft/figureWidth,panelHeightLeft/figureHeight])
panelLeft.set_xlim(0,20)
panelLeft.set_ylim(0,15)
panelLeft.set_xticks([20,0])
panelLeft.set_yticks([0,15,5,10])
panelLeft.set_xticklabels([0,20])

#extract x and y positions from text file
xList = [] 
yList = []
for line in open(input):
    splitLine = line.strip().split("\t") #cleans each line
    xList.append(float(splitLine[1])) #pulls x values
    yList.append(float(splitLine[2])) #pulls y values

#convert to log2
xArray=np.array(xList)
yArray=np.array(yList)
logxArray=np.log2(xArray+1)
logyArray=np.log2(yArray+1)

#plotting main panel
panelMain.plot(logxArray, logyArray, marker="o", linestyle="", markersize=3, color=iBlue, alpha=0.1, markeredgewidth=0)

#plotting top panel
bins=np.arange(0, 15, 0.5)
topHisto, bins=np.histogram(logxArray,bins)
for i in range(0,len(topHisto),1):
    left=bins[i]
    bottom=0
    width=bins[i+1]-left
    height=np.log2(topHisto[i]+1)

    rectangle1=mplpatches.Rectangle((left,bottom),width,height,facecolor=iGreen,linewidth=0.3,edgecolor="black")
    panelTop.add_patch(rectangle1)

#plotting left panel
"""Created using copilot, I promted it by asking for it to
create a histogram by changing my top panel code for the left panel."""
leftHisto,bins = np.histogram(logyArray, bins)
for i in range(len(leftHisto)):
    bottom=bins[i] 
    height=bins[i+1]-bottom 
    width=np.log2(leftHisto[i]+1)
    left=20-width #start position flipped to be on the right side of panel

    # Create a rectangle for each bin
    rectangle = mplpatches.Rectangle((left,bottom),width,height,facecolor=Grey, linewidth=0.3, edgecolor="black")
    panelLeft.add_patch(rectangle)

# save figure
plt.savefig(outFile, dpi=600)