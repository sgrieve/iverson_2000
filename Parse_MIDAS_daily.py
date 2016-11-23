# -*- coding: utf-8 -*-
"""
This script parses data derived from MIDAS weather station data
Created on Tue Nov 22 16:38:15 2016

@author: smudd
"""

import pandas as pd


def load_MIDAS_data():
    #fname = "station_data-196101010000-196112312359.csv"
    #fname = "station_data-201001010000-201611161701.csv"
    fname = "new_small.csv"
    dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M')    
    
    # Read in the data
    MIDAS_df  = pd.read_csv(fname,parse_dates=True,date_parser=dateparse)
    
    # get rid of stupid whitespace in the header
    MIDAS_df.columns = MIDAS_df.columns.str.strip()

    # Make sure that whitespace has been stripped
    a = list(MIDAS_df.columns.values)
    print "Headers are:"    
    print a
    
    # Get rid of some useless columns
    MIDAS_df.drop('id', axis=1, inplace=True)
    MIDAS_df.drop('id_type', axis=1, inplace=True)
    MIDAS_df.drop('met_domain_name', axis=1, inplace=True)
    MIDAS_df.drop('ob_end_ctime', axis=1, inplace=True)
    MIDAS_df.drop('version_num', axis=1, inplace=True)
    MIDAS_df.drop('ob_day_cnt_q', axis=1, inplace=True)
    MIDAS_df.drop('meto_stmp_time', axis=1, inplace=True)
    MIDAS_df.drop('midas_stmp_etime', axis=1, inplace=True)
    MIDAS_df.drop('prcp_amt_j', axis=1, inplace=True)
    
    # Parse dates
    MIDAS_df['ob_date'] =  pd.to_datetime(MIDAS_df['ob_date'], format='%Y-%m-%d %H:%M')  
    
    # get rid of exactly duplicate lines    
    MIDAS_df = MIDAS_df.drop_duplicates()
    

    # Make sure we are only dealing with 1 day records (not monthly totals)
    MIDAS_df = MIDAS_df[MIDAS_df.ob_day_cnt == 1]
    #print yo['ob_date']

    # Make a timestamp for the year 1900
    yr_1900 = pd.Timestamp('1900-01-01')

    # Now we are going to have to group by station
    # print to file, organised by station. 
    for station, df_station in MIDAS_df.groupby('src_id'):
        fname = "MidasNEW_"+str(station)+".csv"

        df_station = df_station.drop_duplicates(['ob_date'], keep="last")
        new_MIDAS = df_station.copy()
        
        # get a column that has the days since 1900
        new_MIDAS['days_since_1900'] = (new_MIDAS['ob_date'] - yr_1900).dt.days        
        
        new_MIDAS.to_csv(fname)
    
    
    
if __name__ == "__main__":
    #compare_linear_to_loop()
    #test_FoS()
    load_MIDAS_data()

