# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 09:59:31 2019

@author: Anurag
This code reads the WDM file and outputs TAUCD and TAUCS values 
for each reach in the model.
"""

from wdmtoolbox import wdmtoolbox as wdm
import pandas as pd

"""Locate the wdm file"""
wdmfile='C:\BASINS45\modelout\LowWill_APGIS\LowerWil.wdm'
DataList=wdm.listdsns(wdmfile)

RCHRES_List=['R:101','R:103','R:105','R:107','R:109','R:111','R:113','R:115','R:117',
'R:119','R:121','R:123','R:125','R:127','R:129','R:131','R:133','R:135','R:137',
'R:139','R:141','R:143','R:201','R:203','R:205','R:207','R:209','R:211','R:213',
'R:215','R:217','R:219','R:221','R:223','R:225','R:227','R:229','R:231','R:233',
'R:235','R:237','R:301','R:303','R:305','R:307','R:309','R:401','R:403','R:405',
'R:407','R:409','R:411','R:413','R:415','R:417']

TauAnalysisText=''
for key in DataList.keys():
    data=DataList[key]
    Location=data['location'].decode('utf-8')
    Constituent=data['constituent'].decode('utf-8')
    if Constituent=='TAU':
        dataSeries=wdm.extract(wdmfile,key)
        TauAnalysisText+= Location +','+ '{:.3E}'.format(dataSeries.quantile(0.01)[0]) + \
                    ','+ '{:.3E}'.format(dataSeries.quantile(0.02)[0]) + \
                    ','+ '{:.3E}'.format(dataSeries.quantile(0.98)[0]) + \
                    ','+ '{:.3E}'.format(dataSeries.quantile(0.99)[0]) + '\n'
f=open('C:\BASINS45\modelout\LowWill_APGIS\LowWillTau.txt','w+')
f.write(TauAnalysisText)
f.close()


