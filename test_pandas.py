# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:38:15 2016

@author: smudd
"""

import numpy as np
import pandas as pd


def load_MIDAS_data():
    #fname = "station_data-196101010000-196112312359.csv"
    fname = "station_data-201001010000-201611161701.csv"
    dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M')    
    
    
    yo  = pd.read_csv(fname,parse_dates=True,date_parser=dateparse)

    a = list(yo.columns.values)
    
    #print a  


    yo.columns = yo.columns.str.strip()

    b = list(yo.columns.values)
    
    #print b  
    
    
    #print yo
    
    
    #print yo['ob_date']
    
    yo['ob_date'] =  pd.to_datetime(yo['ob_date'], format='%Y-%m-%d %H:%M')  
    yo = yo.drop_duplicates()
    yo = yo[yo.ob_day_cnt == 1]
    #print yo['ob_date']
    
    # nested organisation
    for station, df_station in yo.groupby('src_id'):
        print(station)
        print(df_station)

    
    
    
if __name__ == "__main__":
    #compare_linear_to_loop()
    #test_FoS()
    load_MIDAS_data()

