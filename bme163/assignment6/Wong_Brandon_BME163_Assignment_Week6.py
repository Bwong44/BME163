import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import matplotlib
import matplotlib.image as mplimg
import argparse

#initialize the argument parser and style sheet
parser = argparse.ArgumentParser()
parser.add_argument("--outFile", "-o" ,type=str,action="store",help="output file")
parser.add_argument("--phase", "-p" ,type=str,action='store',help="phase file")
parser.add_argument("--exp", "-e" ,type=str,action='store',help="exp file")
parser.add_argument("--gene", "-g" ,type=str,action='store',help="gene list")
args = parser.parse_args()
outFile=args.outFile
phase=args.phase
exp=args.exp
gene=args.gene

#Delete beofre sending#############################################################
from matplotlib import rcParams                                                   #
rcParams['font.family'] = 'Arial'                                                 #
rcParams['font.sans-serif'] = ['Arial']                                           #
matplotlib.font_manager.fontManager.addfont("/hb/home/bwong44/bme163/ARIAL.TTF")  #
###################################################################################

plt.style.use("BME163")
#Set up figure
figureWidth=5
figureHeight=3
plt.figure(figsize=(figureWidth,figureHeight))

def colormap():
    """Creates a colormap using the viridis color scheme"""
    viridis5 = (253/255, 231/255, 37/255)
    viridis4 = (94/255, 201/255, 98/255)
    viridis3 = (33/255, 145/255, 140/255)
    viridis2 = (59/255, 82/255, 139/255)
    viridis1 = (68/255, 1/255, 84/255)
    R1=np.linspace(viridis1[0],viridis2[0],26)
    G1=np.linspace(viridis1[1],viridis2[1],26)
    B1=np.linspace(viridis1[2],viridis2[2],26)

    R2=np.linspace(viridis2[0],viridis3[0],26)
    G2=np.linspace(viridis2[1],viridis3[1],26)
    B2=np.linspace(viridis2[2],viridis3[2],26)

    R3=np.linspace(viridis3[0],viridis4[0],26)
    G3=np.linspace(viridis3[1],viridis4[1],26)
    B3=np.linspace(viridis3[2],viridis4[2],26)

    R4=np.linspace(viridis4[0],viridis5[0],26)
    G4=np.linspace(viridis4[1],viridis5[1],26)
    B4=np.linspace(viridis4[2],viridis5[2],26)

    R=np.concatenate((R1[:-1],R2[:-1],R3[:-1],R4),axis=None)
    G=np.concatenate((G1[:-1],G2[:-1],G3[:-1],G4),axis=None)
    B=np.concatenate((B1[:-1],B2[:-1],B3[:-1],B4),axis=None)

    #for lop for creating the full color map
    colormap=[]
    for i in range(0,101,1):
        colormap.append((R[i],G[i],B[i]))
    return colormap

def makePanel(xpos,ypos,panelWidth, panelHeight, figureWidth, figureHeight, xlabel="", ylabel="", hideTicks=False, title=""):
    """Creates a panel in the figure with specified parameters"""
    relativePanelWidthCenter=panelWidth/figureWidth
    relativePanelHeightCenter=panelHeight/figureHeight
    panel=plt.axes([xpos/figureWidth,ypos/figureHeight,relativePanelWidthCenter,relativePanelHeightCenter],frameon=True)
    panel.set_xlim(0, 8)
    panel.set_ylim(0, 1262) #number of genes from file
    panel.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
    panel.set_xticklabels(["0", "", "6", "", "12", "", "18", ""])
    panel.set_yticks([])
    panel.set_xlabel("CT")
    
    return panel

def makeColormapPanel(xpos,ypos,panelWidth, panelHeight, figureWidth, figureHeight, xlabel="", ylabel="", hideTicks=False):
    """Creates a panel in the figure with specified parameters"""
    relativePanelWidthCenter=panelWidth/figureWidth
    relativePanelHeightCenter=panelHeight/figureHeight
    panel=plt.axes([xpos/figureWidth,ypos/figureHeight,relativePanelWidthCenter,relativePanelHeightCenter],frameon=True)
    panel.set_xlim(0,10)
    panel.set_ylim(0, 101)
    panel.set_yticks([0, 101])
    panel.set_yticklabels(["Min", "Max"])
    if hideTicks:
        panel.set_xticks([])
    panel.set_xlabel(xlabel) #x-axis labels
    panel.set_ylabel(ylabel) #y-axis labels
    panel.tick_params(labelright=True, right=True, labelleft=False, left=False)
    
    #create gradient plot
    colormapValues=colormap()
    for i in range(0,101,1):
        rectangle1=mplpatches.Rectangle((0,i),20,1,
                                    facecolor=colormapValues[i],
                                    linewidth=0,
                                    edgecolor="black"
        )
        panel.add_patch(rectangle1)

    return panel

def readInputData(inputFile, phaseFile):
    """Parse the input exp file, normalize data, and sort by phase"""
    geneDict = {}
    with open(inputFile) as f:
        next(f)  # Skip the header line
        for line in f:
            lineSplit=line.strip().split("\t")
            # Convert values to int
            values=[int(x) for x in lineSplit[4:12]]
            normalizedValues=[int(((v-min(values))/(max(values)-min(values)))*100) for v in values]
            geneDict[lineSplit[1]]=(lineSplit[0], normalizedValues) 
    
    #opens phase file and creates a dict
    phaseDict={}
    with open(phaseFile) as f:
        next(f)  # Skip the header line
        for line in f:
            lineSplit=line.strip().split("\t")
            ensemblID=lineSplit[0]
            phaseDict[ensemblID]=float(lineSplit[1])

    #sort geneDict by phase
    sortedGeneList=sorted(geneDict.items(), key=lambda x: phaseDict[x[0]], reverse=True)

    return sortedGeneList

def plotGeneData(geneList, panel):
    """Plot the gene data on the panel"""
    if gene is None: #handles if no -g is given
        pass
    else:
        geneNames=gene.split(",")
    yTicks=[]
    yLabels=[]
    #panel.set_yticks([])
    for index, (ensembleID, (geneName, values)) in enumerate(geneList): #iterates through the gene list (y-axis)
        for i in range(len(values)): #iterates throught the normalized values (x-axis)
            value=values[i]    
            rectangle=mplpatches.Rectangle([i,index],1,1, #i=x position, index=y position
                                           facecolor=colormap()[value],
                                           edgecolor="black",
                                           linewidth=0,
                                           fill=True
                                           ) 
            panel.add_patch(rectangle)
        if  gene is not None and geneName in geneNames: #checks if the gene name is in the gene list and stores its position
            yTicks.append(index)
            yLabels.append(geneName)
    panel.set_yticks(yTicks)
    panel.set_yticklabels(yLabels)

def main():
    """Main function to run the script"""
    panel1=makePanel(0.7, 0.3, 0.75, 2.5, figureWidth, figureHeight, "CT")
    makeColormapPanel(1.5, 1.45, 0.1, 0.2, figureWidth, figureHeight,"","",hideTicks=True)
    geneList=readInputData(exp, phase)
    plotGeneData(geneList, panel1)
main()

# save figure
plt.savefig(outFile, dpi=600)

