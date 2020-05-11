# -*- coding: utf-8 -*-
"""
Created on Mon May 11 08:35:47 2020

@author: Anurag
"""


import numpy as np


def PanEvaporationValueComputedByHamon(aTAVC, aMonth, aDay, 
                                       aLatDeg, aDegF):
    '''
    This function is adapted from BASINS code. It takes a value or a timeseries
    of daily average temperature as input, and calculates daily potential 
    evapotranspiration using Hamon Method.
    To understand the mathematical basis for this calculation, please
    refer to BASINS documentation.
    
    :param float aTAVC: Average air temperature value or series in degF or degC
    :param integer aMonth: Month of the year
    :param integer aDay: Day of the month
    :param float aLatDeg: Latitude of the loation in Degrees
    :param float aDegF: True if temperature is in degF, False if in degC
    '''
    
    MetComputeLatitudeMin=-66.5
    MetComputeLatitudeMax=66.5
    #The PEVT and disaggregation calculations are valid for this range 
    #of Latitude and Longitude only.
    if aLatDeg < MetComputeLatitudeMin or aLatDeg > MetComputeLatitudeMax: #invalid latitude 
        print('Latitude value is not within range!')
        return -999
        
    else: #latitude ok,convert to radians
        
        JulDay = 30.5 * (aMonth - 1) + aDay
        LatRdn = np.radians(aLatDeg)
        Phi = LatRdn
        AD = 0.40928 * np.cos(0.0172141 * (172.0 - JulDay))
        SS = np.sin(Phi) * np.sin(AD)
        CS = np.cos(Phi) * np.cos(AD)
        X2 = -SS / CS
        Delt = 7.6394 * (1.5708 - np.arctan(X2 / np.sqrt(1.0 - np.power(X2,2))))
        SunR = 12.0 - Delt / 2.0
        SUNS = 12.0 + Delt / 2.0
        DYL  = (SUNS - SunR) / 12

        #convert temperature to degC, if necessary
        if aDegF: aTAVC = (aTAVC - 32.0) * (5.0 / 9.0)

        #Hamon equation
        VPSAT = 6.108 * np.exp(17.26939 * aTAVC / (aTAVC + 237.3))
        VDSAT = 216.7 * VPSAT / (aTAVC + 273.3)
        
        aCTS=np.where(aMonth>0,0.0055,0)
        #make aCTS as a series based on aMonth
        #The BASINS applications I have seen all use the value of 0.0055
        #for all the months
        
        lPanEvap = aCTS * DYL * DYL * VDSAT
        #when the estimated pan evaporation is negative
        #the value is set to zero
        lPanEvap.where(lPanEvap<0,0)
        
        return lPanEvap


def PETDST(aDayPet, aLatDeg, aMonth, aDay):
    
    
    ''' This function is adapted from BASINS code. It takes a value of daily 
    potential evapotranspiration and outputs a list of 24 values for the 
    day.
    To understand the mathematical basis for this calculation, please
    refer to BASINS documentation.
    
    :param float aDayPet: Average PET in inches
    :param float aLatDeg: Latitude of the loation in Degrees
    :param integer aMonth: Month of the year
    :param integer aDay: Day of the month 
    '''
    MetComputeLatitudeMin=-66.5
    MetComputeLatitudeMax=66.5
    #The PEVT and disaggregation calculations are valid for this range 
    #of Latitude and Longitude only.


    #julian date
    JulDay = 30.5 * (aMonth - 1) + aDay
    
    #check latitude
    if aLatDeg < MetComputeLatitudeMin or \
        aLatDeg > MetComputeLatitudeMax: #invalid latitude, return
        aRetCod = -1
        print(aRetCod)
    else: #latitude ok
        #convert to radians
        LatRdn = np.radians(aLatDeg)

        Phi = LatRdn
        AD = 0.40928 * np.cos(0.0172141 * (172.0 - JulDay))
        SS = np.sin(Phi) * np.sin(AD)
        CS = np.cos(Phi) * np.cos(AD)
        X2 = -SS / CS
        Delt = 7.6394 * (1.5708 - np.arctan(X2 / np.sqrt(1.0 - np.power(X2,2))))
        SunR = 12.0 - Delt / 2.0

        #develop hourly distribution given sunrise,
        #sunset and length of day (DELT)
        DTR2 = Delt / 2.0
        DTR4 = Delt / 4.0
        CRAD = 0.66666667 / DTR2
        SL = CRAD / DTR4
        TRise = SunR
        TR2 = TRise + DTR4
        TR3 = TR2 + DTR2
        TR4 = TR3 + DTR4
        aHrPET=24*[0]
        CURVE=24*[0]
        #calculate hourly distribution curve
        for IK in range(24):
            #print(IK)
            RK = IK
            if RK > TRise:
                if RK > TR2:
                    if RK > TR3:
                        if RK > TR4:
                            CURVE[IK] = 0.0
                            aHrPET[IK]=CURVE[IK]
                        else:
                            CURVE[IK] = (CRAD - (RK - TR3) * SL)
                            aHrPET[IK] = CURVE[IK] * aDayPet
                    else:
                        CURVE[IK] = CRAD
                        aHrPET[IK] = CURVE[IK] * aDayPet
                else:
                    CURVE[IK] = (RK - TRise) * SL
                    aHrPET[IK] = CURVE[IK] * aDayPet
            else:
                CURVE[IK] = 0.0
                aHrPET[IK] = CURVE[IK]
            #print('aHRPET Value', IK, aHrPET[IK])
            #print('CURVE Value', IK, CURVE[IK])
            if aHrPET[IK] > 40:
                print("Bad Hourly Value ", aHrPET[IK])
            

    return aHrPET
