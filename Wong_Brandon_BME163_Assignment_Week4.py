import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import matplotlib
import argparse

#initialize the argument parser and style sheet
parser = argparse.ArgumentParser()
parser.add_argument("--outFile", "-o" ,type=str,action="store",help="output file")
parser.add_argument("--identityFile", "-i" ,type=str,action='store',help="identity file")
parser.add_argument("--coverageFile", "-c" ,type=str,action='store',help="coverage file")
args = parser.parse_args()
outFile=args.outFile
identity=args.identityFile
coverage=args.coverageFile

#Delete beofre sending#############################################################
from matplotlib import rcParams                                                   #
rcParams['font.family'] = 'Arial'                                                 #
rcParams['font.sans-serif'] = ['Arial']                                           #
matplotlib.font_manager.fontManager.addfont("/hb/home/bwong44/bme163/ARIAL.TTF")  #
from customPlots import boxPlot #deelete after                                    #
###################################################################################
plt.style.use("BME163")

#Set up figure
figureWidth=6
figureHeight=2.5
plt.figure(figsize=(figureWidth,figureHeight))

#colors
iBlue=(44/255,86/255,134/255)
iGreen=(32/255,100/255,113/255)
iYellow=(248/255,174/255,51/255)
iOrange=(230/255,87/255,43/255)

#set up panel
plt.figure(figsize=(figureWidth,figureHeight))
panelWidthCenter=4.5
panelHeightCenter=1.5
relativePanelWidthCenter=panelWidthCenter/figureWidth
relativePanelHeightCenter=panelHeightCenter/figureHeight
panelCenter=plt.axes([0.5/figureWidth,0.15,relativePanelWidthCenter,relativePanelHeightCenter])
panelCenter.set_xlim(0,300)
panelCenter.set_ylim(75,100)
panelCenter.set_xticks([37.5,112.5,187.5,262.5])
panelCenter.set_xticklabels(["1-3","4-6","7-9",">=10"])
panelCenter.set_yticks([75,80,85,90,95,100])
panelCenter.set_xlabel("Subread Coverage") #x-axis labels
panelCenter.set_ylabel("Identity (%)") #y-axis labels

#read in coverage file
coverageDict={}
for line in open(coverage):
    splitLine=line.strip().split()
    coverageDict[splitLine[0]]=int(splitLine[1])

#read in identity file and bin data to 500 subsample per
bin1=[] #1-3
bin2=[] #4-6
bin3=[] #7-9
bin4=[] #>=10
for line in open(identity): #interages through identity file
    splitLine=line.strip().split() #cleans each line
    if splitLine[0] in coverageDict: #correlates the files by their name
        if int(coverageDict[splitLine[0]]) <= 3: #puts in 1-3 bin
            bin1.append(float(splitLine[1]))
        elif 4 <= int(coverageDict[splitLine[0]]) <= 6: #puts in 4-6 bin
            bin2.append(float(splitLine[1]))
        elif 7 <= int(coverageDict[splitLine[0]]) <= 9: #puts in 7-9 bin
            bin3.append(float(splitLine[1]))
        elif int(coverageDict[splitLine[0]]) >= 10: #puts in >=10 bin
            bin4.append(float(splitLine[1]))

#plot swarm plot
def plotSwarm(panel,yvalues,position,panelWidth,pandelHeight,xmin,xmax,ymin,ymax,ms,color,label): #ms is marker size
    """Plotting function for swarm plot"""
    #defining the parameters and variables
    xrange=xmax-xmin
    yrange=ymax-ymin
    markerSize=ms
    minDistance=markerSize/72
    increment=((minDistance/10)*xrange)/panelWidth
    span=(0.51*xrange)/panelWidth

    #find all possible positions
    possiblePositions=[]
    for shift in np.arange(0,span,increment):
        possiblePositions.append(position+shift)
        possiblePositions.append(position-shift)

    #plot subsample of points
    plottedPoints=[]
    for y1 in yvalues[:500]:
        if len(plottedPoints)==0: #adds in the first point
            plottedPoints.append((position,y1))
        else: #else if plottedPoints list isn't empty
            pointPlaced=False #true false statements from copilot
            for x1 in possiblePositions:
                distList=[]
                for x2,y2 in plottedPoints: #iterates through the points already plotted
                    xdist=((x2-x1)/xrange)*panelWidth
                    ydist=((y2-y1)/yrange)*pandelHeight
                    distance=((xdist**2)+(ydist**2))**0.5
                    distList.append(distance)
                if min(distList) > minDistance:
                    plottedPoints.append((x1,y1))
                    pointPlaced=True
                    break
            if not pointPlaced:  # If no valid position is found, stop plotting
                print(f"{500-len(plottedPoints)} points could not be plotted at x position {label}.")
                break
                

    for x1,y1 in plottedPoints: #plots points once loop has completed using the plottedPoints list
        panel.plot(x1,y1,marker="o",ms=markerSize,mew=0,lw=0,color=color)    

plotSwarm(panelCenter,bin1,37.5,4.5,1.5,0,300,75,100,2,iBlue,"1-3") #1-3
plotSwarm(panelCenter,bin2,112.5,4.5,1.5,0,300,75,100,2,iGreen,"4-6") #4-6
plotSwarm(panelCenter,bin3,187.5,4.5,1.5,0,300,75,100,2,iYellow,"7-9") #7-9
plotSwarm(panelCenter,bin4,262.5,4.5,1.5,0,300,75,100,2,iOrange,">=10") #>=10

# save figure
plt.savefig(outFile, dpi=600)