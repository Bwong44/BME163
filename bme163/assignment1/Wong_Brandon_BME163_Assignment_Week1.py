import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import matplotlib
import argparse

#initialize the argument parser and style sheet
parser = argparse.ArgumentParser()
parser.add_argument("--outFile","-o" ,type=str,action="store",help="output file")
args = parser.parse_args()
outFile=args.outFile
plt.style.use("BME163")


#set up figure
figureWidth=5
figureHeight=2
plt.figure(figsize=(figureWidth,figureHeight))

#set up panels
panelLeft1=0.20
panelBottom=0.20
panelWidth=2
panelHeight=1
panel1=plt.axes([panelLeft1/figureWidth,panelBottom/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
panel1.set_xticks([])
panel1.set_yticks([])
panel1.set_xlim(0,101)
panel1.set_ylim(0,200)

#plotting data
#heatmap1:
viridis5 = (253/255, 231/255, 37/255)
viridis4 = (94/255, 201/255, 98/255)
viridis3 = (33/255, 145/255, 140/255)
viridis2 = (59/255, 82/255, 139/255)
viridis1 = (68/255, 1/255, 84/255)

#heatmap2:
plasma5 = (237/255, 252/255, 27/255)
plasma4 = (245/255, 135/255, 48/255)
plasma3 = (190/255, 48/255, 101/255)
plasma2 = (87/255, 0/255, 151/255)
plasma1 = (15/255, 0/255, 118/255)

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

#colormap2
R1Plasma=np.linspace(plasma1[0],plasma2[0],26)
G1Plasma=np.linspace(plasma1[1],plasma2[1],26)
B1Plasma=np.linspace(plasma1[2],plasma2[2],26)

R2Plasma=np.linspace(plasma2[0],plasma3[0],26)
G2Plasma=np.linspace(plasma2[1],plasma3[1],26)
B2Plasma=np.linspace(plasma2[2],plasma3[2],26)

R3Plasma=np.linspace(plasma3[0],plasma4[0],26)
G3Plasma=np.linspace(plasma3[1],plasma4[1],26)
B3Plasma=np.linspace(plasma3[2],plasma4[2],26)

R4Plasma=np.linspace(plasma4[0],plasma5[0],26)
G4Plasma=np.linspace(plasma4[1],plasma5[1],26)
B4Plasma=np.linspace(plasma4[2],plasma5[2],26)

RPlasma=np.concatenate((R1Plasma[:-1],R2Plasma[:-1],R3Plasma[:-1],R4Plasma),axis=None)
GPlasma=np.concatenate((G1Plasma[:-1],G2Plasma[:-1],G3Plasma[:-1],G4Plasma),axis=None)
BPlasma=np.concatenate((B1Plasma[:-1],B2Plasma[:-1],B3Plasma[:-1],B4Plasma),axis=None)

#for loop for creating color map 2
colormap2=[]
for i in range(0,101,1):
    colormap2.append((RPlasma[i],GPlasma[i],BPlasma[i]))

#for loop for plotting the rectangles in steps of 1, correlating with the color map
for i in range(0,101,1):
    rectangle1=mplpatches.Rectangle((i,100),1,200,
                                facecolor=colormap[i],
                                linewidth=0,
                                edgecolor="black"

    )
    panel1.add_patch(rectangle1)

    panel1.plot(i,i, #draws the diagonal line from the bottom
         marker="o",
         markeredgewidth=0,
         markerfacecolor=colormap2[i],
         color="black",
         linewidth=1,
         linestyle="--",
         markersize=4,
         zorder=0
    )

#plot a line down the middle from copilot
panel1.plot([0, 101], [100, 100], color='black', linewidth=1)

# save figure
plt.savefig(outFile,dpi=600)

