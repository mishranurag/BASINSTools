# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 18:30:58 2019

@author: Anurag
This program is used to calculate FTABLE for storm drains and natural drainages

"""
import os
import pandas as pd
import numpy as np
from io import StringIO

datafile='C:\BASINS45\modelout\LowWill_APGIS\Hydraulics\RCHRES_Details_PVT_Drainages.csv'
#datafile='C:\BASINS45\modelout\LowWill_APGIS\Hydraulics\RCHRES_Details_SW_CSO.csv'
filelocation,filename=os.path.split(datafile)
OutputFile=os.path.splitext(filename)[0]+'_FTABLE.txt'
Hydraulicdata=pd.read_csv(datafile)        


testdata='SUBBASINS,TAREAMI2,Type,Pipe_Dia_ft,Manning n,FloodPlain n,SLO2_mean\n'\
'101,0.05,Sewer,2,0.024,,0.01\n'\
'102,0.1,Natural,,0.1,0.2,0.02\n'
Hydraulicdata=pd.read_csv(StringIO(testdata))


def SidePiece(aDepth, aDepthSlopeChange, aSideSlopeUpperFP, \
              aDepthChannel, aSideSlopeLowerFP, aSlopeSide):
    if aDepth>aDepthSlopeChange:
        lSidePiece=(aDepth - aDepthSlopeChange) / aSideSlopeUpperFP
    elif(aDepth > aDepthChannel):
        lSidePiece = (aDepth - aDepthChannel) / aSideSlopeLowerFP
    else:
        lSidePiece = aDepth / aSlopeSide
    return lSidePiece
ftableAll=''
for index,row in Hydraulicdata.iterrows():
    ftableheader='  FTABLE    ' + str(row['SUBBASINS']) +'\n'
    ftable=''
    if row['Type']=='Sewer':            
        DT=row['Pipe_Dia_ft']
        NP=row['Manning n']
        SP=row['SLO2_mean']/100.0
        L=row['LEN2']*3.28
        if DT<=1:
            dblDepthIncrement=DT/10.0
        else:
            dblDepthIncrement = 0.1;
        dblDepth=0.0
        nrows=0
        ftableheader+='***Sewer Pipe With Diameter={:0.1f} and Manning n={:0.3f}\n'.format(DT,NP)
        ftableDict={}
        while dblDepth<=DT:
            #Flow area calculations
            nrows+=1
            b1 = DT * DT / 8.0
            b2 = 1.0 - (2.0 * dblDepth / DT)
            b3 = 2.0 * (2.0 * np.arctan(1.0) - np.arctan(b2 / np.sqrt(1.0 - b2 * b2)))
            b4 = np.sin(b3)
            b5 = 2.0 * (2.0 * np.arctan(1.0) - np.arctan(b2 / np.sqrt(1.0 - b2 * b2)))
            R=DT/2
            dblArea = b1 * (b5 - b4)
            
            #top width calculation
            
            t1 = ((2.0 * dblDepth - DT) / DT)
            t2 = 2.0 * np.arctan(1.0) - np.arctan(t1 / np.sqrt(1.0 - t1 * t1))
            t = DT * np.sin(t2)
            wp1 = (1.0 - (2.0 * dblDepth / DT))
            wP = DT * (2.0 * np.arctan(1.0) - np.arctan(wp1 / np.sqrt(1.0 - wp1 * wp1)))
            CONS=1.486
            ACONS=43560
            dblOutFlow = CONS * dblArea * np.power((dblArea / wP), (2.0 / 3.0)) * np.sqrt(SP) / NP;
            #Volume calculations
            V1 = (R - dblDepth) * np.sqrt(2.0 * R * dblDepth - dblDepth * dblDepth);
            V2 = (R - dblDepth) / R;
            V3 = 2.0 * np.arctan(1.0) - np.arctan(V2 / np.sqrt(1.0 - V2 * V2));
            dblVolume = (L * ((0.5 * DT) * (0.5 * DT) * V3 - V1)) / ACONS;
            #Velocity Calculatios
            dblVelocity = CONS * np.power((dblArea / wP), (2.0 / 3.0)) * np.sqrt(SP) / NP;
            #Area Calculations
            dblArea = (L * t) / ACONS;
            if (dblDepth == 0.0): 
                dblOutFlow = 0.0
                dblArea=0.0
            dblOutFlow2=0.0
            PreviousDepth=round(dblDepth-dblDepthIncrement,3)
            if len(ftableDict) > 2 and dblOutFlow<ftableDict[PreviousDepth][2]:
                dblOutFlow2=ftableDict[PreviousDepth][2]-dblOutFlow
                dblOutFlow=ftableDict[PreviousDepth][2]*1.1
            
            
            
            ftableDict[dblDepth]=[dblArea, dblVolume, dblOutFlow,dblOutFlow2]
            ftable+=" {:9.2g} {:9.4g} {:9.4g} {:9.4g} {:9.4g}". \
                    format(dblDepth, dblArea, dblVolume, dblOutFlow,dblOutFlow2)+'\n'
            if dblDepth!=0.0: 
                dblDepthIncrement=DT/10.0
            if np.isclose(dblDepth,DT):
                print('I am here for Subbasin:' + str(row['SUBBASINS']))
                break
            dblDepth+=dblDepthIncrement
            dblDepth=round(dblDepth,3)
            #print(str(row['SUBBASINS']) + ' : ' + str(dblDepthIncrement))
        
        SArea=0.0
        for key in ftableDict.keys():
            if SArea<ftableDict[key][1]:SArea=ftableDict[key][1]
        dblArea=SArea*10    
        
        
        ftable+=" {:9.2g} {:9.4g} {:9.4g} {:9.4g} {:9.4g}".format(dblDepth*1.5, 
                 dblArea, dblVolume*10, dblOutFlow,dblOutFlow*10)+'\n'    
                 
        ftableheader+=' rows cols                               ***' + '\n'
        ftableheader+="{:5.0f}{:5.0f}".format(nrows+1, 5)+'\n'
        ftableheader+='     depth      area    volume  outflow1   outflow2***\n'
        ftable=ftableheader+ftable
        ftable+='  END FTABLE' + str(row['SUBBASINS'])+'\n'
        ftable+='\n'
        
    elif row['Type']=='Natural':
        D_A_sq_km=row['TAREAMI2']*1.6*1.6
        Wid_m=1.29*D_A_sq_km**0.6
        Dep_m=0.13*D_A_sq_km**0.4
        Len_m=row['LEN2']
        SL_PCT=row['SLO2_mean']
        Stream_n=row['Manning n']
        FP_n=row['FloodPlain n']
        slopeSideUpperFPLeft=0.5
        slopeSideUpperFPRight=0.5
        WidthZeroSlopeLeft=Wid_m*3.28
        slopeSideLeft=1
        slopeSideRight=1
        WidthZeroSlopeRight=Wid_m*3.28
        slopeSideLowerFPLeft=0.5
        slopeSideLowerFPRight=0.5
        ChannelDep=Dep_m*3.28*1.25
        MaxDep=Dep_m*3.28*62.5
        SlopeChange=Dep_m*3.28*1.875
        ftableheader+='***Natural Channel. Avg Width={:0.2f}m, Avg Depth={:0.2f}m, \n'.format(Wid_m,Dep_m)
        ftableheader+='***Stream Manning n={:0.3f}, and Floodplain n={:0.3f} \n'.format(Stream_n, FP_n)
        
        
        if Dep_m*3.28 < ChannelDep:
            dWb = Wid_m*3.28 - (Dep_m*3.28/slopeSideLeft) \
                -(Dep_m*3.28/slopeSideRight)
        if Dep_m*3.28 > ChannelDep and Dep_m*3.28 < SlopeChange:
            dWb = Wid_m*3.28 - WidthZeroSlopeLeft - WidthZeroSlopeRight \
            - ((Dep_m*3.28 - ChannelDep) / slopeSideLowerFPLeft) - \
            ((Dep_m*3.28 - ChannelDep) / slopeSideLowerFPRight) - \
            (ChannelDep/slopeSideLeft) - (ChannelDep / slopeSideRight)
        if Dep_m*3.28 > SlopeChange and Dep_m*3.28 < MaxDep:
            dWb = Wid_m*3.28 - ((Dep_m*3.28 - SlopeChange)\
                                / slopeSideUpperFPLeft) - ((Dep_m*3.28 - \
                                 SlopeChange) / slopeSideUpperFPRight) - \
                                WidthZeroSlopeLeft - WidthZeroSlopeRight - \
                                ((SlopeChange - ChannelDep) / \
                                 slopeSideLowerFPLeft) - ((SlopeChange - \
                                ChannelDep) / slopeSideLowerFPRight) - \
                                 (ChannelDep / slopeSideLeft) - \
                                 (ChannelDep / slopeSideRight)
        if dWb < 0: dWb = 0.0001
                
        dWc=dWb+(ChannelDep/slopeSideLeft)+(ChannelDep/slopeSideRight)
        dWt1=dWc+WidthZeroSlopeLeft+WidthZeroSlopeRight+ \
                    ((SlopeChange-ChannelDep)/slopeSideLowerFPLeft+ \
                         (SlopeChange-ChannelDep)/slopeSideLowerFPRight)
        dWt2=dWt1+((MaxDep-SlopeChange)/slopeSideUpperFPLeft + \
                (MaxDep-SlopeChange)/slopeSideUpperFPRight)
            
        lDepth=[]
        lDepth.append(0)
        lDepth.append(Dep_m*3.28/10)
        lDepth.append(Dep_m*3.28)
        lDepth.append(ChannelDep)
        lDepth.append((ChannelDep+SlopeChange)/2)
        lDepth.append(SlopeChange)
        lDepth.append((SlopeChange+MaxDep)/2)
        lDepth.append(MaxDep)
        
        nrows=8
        
        for depth in lDepth:
            
            if depth>SlopeChange:
                lNearestBase=dWt1
            elif depth>ChannelDep:
                lNearestBase=dWc+WidthZeroSlopeLeft+WidthZeroSlopeRight
            else:
                lNearestBase=dWb
            
            if ChannelDep>depth:
                lDepthD=depth
            else:
                lDepthD=ChannelDep
                
            lLeftPiece=SidePiece(lDepthD, SlopeChange, \
                                 slopeSideUpperFPLeft, ChannelDep,\
                                 slopeSideLowerFPLeft, slopeSideLeft)
            lRightPiece=SidePiece(lDepthD, SlopeChange, \
                                 slopeSideUpperFPRight, ChannelDep,\
                                 slopeSideLowerFPRight, slopeSideRight)
            lCrossSectionArea = lDepthD * \
                    (dWb + (lLeftPiece * 0.5) + (lRightPiece * 0.5))
            if depth>ChannelDep:
                if SlopeChange>depth:lDepthD=depth
                else: lDepthD=SlopeChange
                lLeftPiece = SidePiece(lDepthD, SlopeChange, \
                                   slopeSideUpperFPRight, ChannelDep,\
                                   slopeSideLowerFPRight, slopeSideRight)
                lRightPiece = SidePiece(lDepthD, SlopeChange, \
                            slopeSideUpperFPLeft, ChannelDep, \
                            slopeSideLowerFPLeft, slopeSideLeft)
                lCrossSectionArea += (lDepthD - ChannelDep) * (dWc + \
                               WidthZeroSlopeLeft + WidthZeroSlopeRight + \
                               (lLeftPiece * 0.5) + (lRightPiece * 0.5))
                               
            if depth>SlopeChange:
                if MaxDep>depth:lDepthD=depth
                else: lDepthD=MaxDep
                
                lLeftPiece = SidePiece(lDepthD, SlopeChange, \
                                       slopeSideUpperFPRight, ChannelDep, \
                                       slopeSideLowerFPRight, slopeSideRight)
                lRightPiece = SidePiece(lDepthD, SlopeChange,\
                                        slopeSideUpperFPLeft, ChannelDep, \
                                        slopeSideLowerFPLeft, slopeSideLeft)
                lCrossSectionArea += (lDepthD - SlopeChange) * \
                                        (dWt1 + (lLeftPiece * 0.5) + \
                                         (lRightPiece * 0.5))
            #calculating hydraulic Radius
            lDenominator = dWb
            if ChannelDep > depth: lDepthD = depth
            else: lDepthD = ChannelDep
            lDenominator += lDepthD * (np.sqrt(1.0 + \
                                1.0 / (slopeSideLeft * slopeSideLeft)) + \
                                np.sqrt(1.0 + 1.0 / slopeSideRight**2))
            if depth > ChannelDep:
                if SlopeChange > depth: lDepthD = depth
                else: lDepthD = SlopeChange
                lDenominator += WidthZeroSlopeLeft + WidthZeroSlopeRight + \
                    (lDepthD - ChannelDep) * (np.sqrt(1.0 + \
                    1.0 / slopeSideLowerFPLeft**2) + np.sqrt(1.0 + \
                          1.0 / slopeSideLowerFPRight**2))
                        
            if depth > SlopeChange:
                if MaxDep > depth:lDepthD = depth
                else: lDepthD = MaxDep
                lDenominator += (lDepthD - SlopeChange) * (np.sqrt(1.0 + \
                                1.0 / slopeSideUpperFPLeft**2) + \
                                np.sqrt(1.0 + 1.0 / slopeSideUpperFPRight**2))
                        
            lHydraulicRadius = lCrossSectionArea / lDenominator

            lLeftPiece = SidePiece(depth, SlopeChange, slopeSideUpperFPRight, \
                                   ChannelDep, slopeSideLowerFPRight, slopeSideRight)
            lRightPiece = SidePiece(depth, SlopeChange, \
                                    slopeSideUpperFPLeft, ChannelDep, \
                                    slopeSideLowerFPLeft, slopeSideLeft)

            lSurfaceArea = Len_m*3.28 * (lNearestBase + lLeftPiece \
                                        + lRightPiece) / 43560.0
            lVolume = Len_m*3.28 * lCrossSectionArea / 43560.0
            if depth>SlopeChange:ManningN=FP_n
            else:ManningN=Stream_n
            lDischarge = 1.49 / ManningN * (lHydraulicRadius ** (2/3))\
                        * np.sqrt(SL_PCT/100) * lCrossSectionArea
            if depth==0:
                lSurfaceArea=0
                lVolume=0
                lDischarge=0
            ftable+=" {:9.2f} {:9.4g} {:9.4g} {:9.4g}".format(depth, lSurfaceArea, lVolume, lDischarge)+'\n'
        ftableheader+=' rows cols                               ***' + '\n'
        ftableheader+="{:5.0f}{:5.0f}".format(nrows, 4)+'\n'
        ftableheader+='     depth      area    volume  outflow1 ***\n'
        ftable=ftableheader+ftable
        ftable+='  END FTABLE' + str(row['SUBBASINS'])+'\n'
        ftable+='\n'
        #print(ftable)  
    ftableAll+=ftable
f=open(os.path.join(filelocation, OutputFile),'w+')
f.write(ftableAll)
f.close()
                
            
            
                    
            
            
            

          
        
        
        
        
