import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np

def boxPlot(panel,yvalues,position,width,faceColor="white"): #default color is white
    median=np.median(yvalues) #median
    percentile25=np.percentile(yvalues,25) #25th percentile
    percentile75=np.percentile(yvalues,75) #75th percentile
    IQR=percentile75-percentile25
    lowBound=percentile25-(1.5*IQR) #lower bound
    highBound=percentile75+(1.5*IQR) #upper bound
    bottom=min(yvalues) #minimum
    top=max(yvalues) #maximum
    previous=-np.inf
    outliers=[] #outliers
    for element in sorted(yvalues):
        if element < lowBound or element > highBound:
            outliers.append(element)
        if element >= lowBound and previous < lowBound:
            bottom=element
        if element > highBound and previous <=highBound:
            top=previous
        previous=element
    panel.plot([position]*len(outliers),outliers,marker="o",ms=4,mew=0,mfc="black",lw=0,alpha=0.5) #outliers plot
    rectangle=mplpatches.Rectangle([position-(width/2), percentile25],
                                   width,
                                   IQR,
                                   facecolor=faceColor,
                                   edgecolor="black",
                                   lw=0.75,
                                   fill=True)
    panel.add_patch(rectangle)
    panel.plot([position-width/2,position+width/2],[median,median],color="red") #median lin
    panel.plot([position-width/4,position+width/4],[top,top],color="black",zorder=0) #top cap
    panel.plot([position-width/4,position+width/4],[bottom,bottom],color="black",zorder=0) #bottom cap
    panel.plot([position,position],[percentile75,top],color="black",zorder=0) #top line
    panel.plot([position,position],[percentile25,bottom],color="black",zorder=0) #bottom line