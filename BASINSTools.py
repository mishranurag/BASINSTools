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


def ExtractHydrologicStateVariables(aBinaryData, aGroup):
    '''
    This function outputs the state variable table based on the 
    data passed and the group 
    
    :param Series aBinaryData: Values of all group members at the date when 
    the state variables are needed to be extracted
    :param string aGroup: Name of the group
    
    aBinaryData should look like this
    PERLND_21_AGWET_4    0.000000
    PERLND_21_AGWI_4     1.883073
    PERLND_21_AGWO_4     1.006952
    PERLND_21_AGWS_4     1.468427
    PERLND_21_BASET_4    0.128918
    
    aGroup should be 'PWATER' or 'IWATER' or 'HYDR'
    '''
    
    ListOfOperations=[]
    for OPN in aBinaryData.index:
        OPNNumber=OPN.split('_')[1]
        if OPNNumber in ListOfOperations:continue
        ListOfOperations.append(OPNNumber)
    tStep=aBinaryData.index[0][-2:]
    ListOfStateVariables=[]
    if aGroup=='PWATER':
        ListOfStateVariables=['CEPS','SURS','UZS','IFWS','LZS','AGWS','GWVS']
        output='  PWAT-STATE1\n'
        output+='*** < PLS>  PWATER state variables (in)\n'
        output+='*** x  - x      CEPS      SURS       UZS      IFWS       LZS      AGWS      GWVS\n'

        for OPN in ListOfOperations:    
            output+='{:>5}'.format(OPN)
            output+='     '
            for StateVariable in ListOfStateVariables:
                ColumnName='PERLND_'+ str(OPN) + '_' + StateVariable + tStep
                output+='{:10.4f}'.format(aBinaryData[ColumnName])
            output+='\n'
        output+='  END PWAT-STATE1\n'

    elif aGroup=='IWATER':
        ListOfStateVariables=['RETS','SURS']
        output+='  IWAT-STATE1\n'
        output+='*** <ILS >  IWATER state variables (inches)\n'
        output+='*** x -  x      RETS      SURS\n'
        
        for OPN in ListOfOperations:
            output+='{:>5}'.format(OPN)
            output+='     '
            for StateVariable in ListOfStateVariables:
                ColumnName='IMPLND_'+ str(OPN) + '_' + StateVariable + tStep
                output+='{:10.4f}'.format(aBinaryData[ColumnName])
            output+='\n'
        output+='  END IWAT-STATE1\n'

    elif aGroup=='HYDR':
        ListOfStateVariables=['VOL']
        output+='  HYDR-INIT\n'
        output+='***         Initial conditions for HYDR section\n'
        output+='***RC HRES       VOL  CAT Initial value  of COLIND     initial  value  of OUTDGT\n'
        output+='*** x  - x     ac-ft      for each possible   exit  for each possible exit,ft3\n'
        
        
        for OPN in ListOfOperations:
            output+='{:>5}'.format(OPN)
            output+='     '
            for StateVariable in ListOfStateVariables:
                ColumnName='RCHRES_'+ str(OPN) + '_' + StateVariable + tStep
                output+='{:10.4f}'.format(aBinaryData[ColumnName])
                output+='       4.2  4.5  4.5  4.5  4.2       2.1  1.2  0.5  1.2  1.8'
            output+='\n'
        
        output+='  END HYDR-INIT\n'
    else:
        output='Not Implemented'
    return output

    

    
        