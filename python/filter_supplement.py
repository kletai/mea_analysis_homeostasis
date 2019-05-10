import pandas as pd
import itertools as it
import seaborn as sns
import numpy as np
from pymea import matlab_compatibility as mc
from matplotlib import pyplot as plt
from matplotlib import mlab as mlab
import random
from datetime import datetime, timedelta

def filter_neurons_homeostasis(cat_table, baseline_table, stim_table, ind_filter = True, var=10, minHz = 0.001, maxHz = 100000, foldMin = 0.001, filter_wells = False, data_col = 'spike_freq'):
    '''
    Returns a cat_table only including neurons that pass the filters for min/maxHz, baseline "var", staying alive
    throughout the experiment, and responding to drug.
    '''
    c_filter = pd.DataFrame()
    b_filter = pd.DataFrame()
    count_real = 0
    count_live = 0
    count_final = 0
    last_time = max(cat_table['time'])
    #filter individual neurons based on baseline firing, whether they stay alive, and whether they respond to stim
    for cond in cat_table['condition'].unique():
        c = cat_table.query('condition == "%s"'%cond)
        b = baseline_table.query('condition == "%s"'%cond)
        s = stim_table.query('condition == "%s"'%cond)
        c_filter_cond = pd.DataFrame()
        b_filter_cond = pd.DataFrame()
        s_filter_cond = pd.DataFrame()
        for unit_name in c['unit_name'].unique():
            unit = c.query('unit_name == @unit_name')
            unit_b = b.query('unit_name == @unit_name')
            unit_s = s.query('unit_name == @unit_name')
            meanOfBaseline = np.mean(unit_b[data_col])
            varOfBaseline = (unit_b.loc[:,data_col].max() - unit_b.loc[:,data_col].min())/meanOfBaseline
            meanAfterDrug = np.mean(unit_s[6:42][data_col]) #mean firing in first 3hr of stim, skipping first 30min
            folds = unit[data_col]/meanOfBaseline
            folds_b = unit_b[data_col]/meanOfBaseline
            folds_s = unit_s[data_col]/meanOfBaseline
            if meanOfBaseline > minHz and meanOfBaseline < maxHz and varOfBaseline < var:
                count_real = count_real+1
                if (sum((folds) < foldMin)<50):
                    count_live = count_live+1
                    if (ind_filter == False) or (meanOfBaseline < 1*meanAfterDrug): 
                        count_final = count_final+1
                        well_num = mc.get_well_number(unit_name)
                        unit.loc[:,'well'] = well_num
                        unit.loc[:,'folds'] = folds
                        unit_b.loc[:,'well'] = well_num
                        unit_b.loc[:,'folds'] = folds_b
                        unit_s.loc[:,'well'] = well_num
                        unit_s.loc[:,'folds'] = folds_s
                        c_filter_cond = c_filter_cond.append(unit, ignore_index=True)
                        b_filter_cond = b_filter_cond.append(unit_b, ignore_index=True)
                        s_filter_cond = s_filter_cond.append(unit_s, ignore_index=True)
            else:
                continue
                
        if c_filter_cond.empty:
            return (c_filter, b_filter, count_real, count_live, count_final)
    
        #filter out wells that have fewer than two neurons, or whose median doesn't increase by at least 1.3
        for w in c_filter_cond['well'].unique():
            s_well = s_filter_cond.query('well == @w')
            num_units = len(s_well['unit_name'].unique())
            if num_units <= 2:
                c_filter_cond = c_filter_cond[c_filter_cond.well != w]
                b_filter_cond = b_filter_cond[b_filter_cond.well != w]
                count_real = count_real - num_units
                count_live = count_real - num_units
                count_final = count_final - num_units
                print('Fewer than 3 neurons. Omitting well ' + repr(w))
            elif ind_filter == True: 
                well_med = s_well.groupby('time')['folds'].median()
                fold_med = np.mean(well_med[6:42])
                if fold_med < 1.3:
                    c_filter_cond = c_filter_cond[c_filter_cond.well != w]
                    b_filter_cond = b_filter_cond[b_filter_cond.well != w]
                    count_final = count_final - num_units
                    print('Weak median FR response. Omitting well ' + repr(w))
    return (c_filter, b_filter, count_real, count_live, count_final)
