import pandas as pd
import itertools as it
import seaborn as sns
import numpy as np
from pymea_use import matlab_compatibility as mc
from matplotlib import pyplot as plt
from matplotlib import mlab as mlab
import random
from datetime import datetime, timedelta

def select_neurons(cat_table, min_freq = 0, max_freq = 100000000, data_col='spike_freq'):
    '''
    Returns a version of cat_table only containing units whose mean frequency is between min_freq and max_freq
    '''
    unit_freq_mean = cat_table.groupby(('unit_name'))[data_col].mean()
    unit_freq_mean = unit_freq_mean.rename(data_col).reset_index()
    filt = (unit_freq_mean[data_col] >= min_freq) & (unit_freq_mean[data_col] < max_freq)
    selected_units = unit_freq_mean.loc[filt,'unit_name']
    return (cat_table.loc[cat_table['unit_name'].isin(selected_units)], selected_units)

def filter_neurons_homeostasis(cat_table, baseline_table, stim_table, ind_filter = True, var=10, minHz = 0.001, maxHz = 100000, foldMin = 0.001, filter_wells = False, data_col = 'spike_freq', n_units=2):
    '''
    Returns a cat_table only including neurons that pass the filters for min/maxHz, baseline var, staying alive
    throughout the experiment, responding to drug, and whose wells are behaving similarly to others.
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
                if (sum((folds) < foldMin)<50):# and (sum(folds > 500)<100):
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
            if num_units <= n_units:
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
        # filter out wells that behave differently from others
        if filter_wells == True:
            if c_filter_cond.empty:
                return (c_filter_cond, b_filter_cond, count_real, count_live, count_final)
            c_filter_cond, b_filter_cond = drop_outlier_wells(c_filter_cond, b_filter_cond, data_col = data_col)
        
        c_filter = c_filter.append(c_filter_cond, ignore_index=True) #concatenate filtered condition tables
        b_filter = b_filter.append(b_filter_cond, ignore_index=True) #concatenate filtered condition tables
    return (c_filter, b_filter, count_real, count_live, count_final)

def drop_outlier_wells(c_filter, b_filter, data_col = 'spike_freq'):
    '''
    Crops wells whose correlation coefficient with both the mean/median frequency of the entire condition is below 0.6
    '''
    mean_freq_traces = c_filter.groupby(('condition', 'time'))[data_col].mean()
    mean_freq_traces = mean_freq_traces.rename(data_col).reset_index()
    mean_freq_traces_b = b_filter.groupby(('condition', 'time'))[data_col].mean()
    mean_freq_traces_b = mean_freq_traces_b.rename(data_col).reset_index()
    median_freq_traces = c_filter.groupby(('condition', 'time'))[data_col].median()
    median_freq_traces = median_freq_traces.rename(data_col).reset_index()
    median_freq_traces_b = b_filter.groupby(('condition', 'time'))[data_col].median()
    median_freq_traces_b = median_freq_traces_b.rename(data_col).reset_index()
    meanOfMean = np.mean(mean_freq_traces_b[data_col])
    meanOfMedian = np.mean(median_freq_traces_b[data_col])
    mean_folds = mean_freq_traces[data_col]/meanOfMean
    median_folds = median_freq_traces[data_col]/meanOfMedian
    
    mean_by_well = c_filter.groupby(('well', 'time'))[data_col].mean()
    mean_by_well = mean_by_well.rename(data_col).reset_index()
    median_by_well = c_filter.groupby(('well', 'time'))[data_col].median()
    median_by_well = median_by_well.rename(data_col).reset_index()
    mean_by_well_b = b_filter.groupby(('well', 'time'))[data_col].mean()
    mean_by_well_b = mean_by_well_b.rename(data_col).reset_index()
    median_by_well_b = b_filter.groupby(('well', 'time'))[data_col].median()
    median_by_well_b = median_by_well_b.rename(data_col).reset_index()

    for well in mean_by_well['well'].unique():
        this_well_mean = mean_by_well.query('well == @well')
        this_well_mean_b = mean_by_well_b.query('well == @well')
        this_well_median = median_by_well.query('well == @well')
        this_well_median_b = median_by_well_b.query('well == @well')
        meanOfMean_well = np.mean(this_well_mean_b[data_col])
        meanOfMedian_well = np.mean(this_well_median_b[data_col])
        
        mean_folds_well = this_well_mean[data_col]/meanOfMean_well
        median_folds_well = this_well_median[data_col]/meanOfMedian_well
        corr_mean = np.corrcoef(mean_folds_well, mean_folds)
        corr_median = np.corrcoef(median_folds_well, median_folds)
        #print(corr_mean)
        #print(corr_median)
        if corr_mean[0,1] < 0.6 and corr_median[0,1] < 0.6:
            c_filter = c_filter[c_filter.well != well]
            b_filter = b_filter[b_filter.well != well]
            print('Omitting well ' + repr(well))
    return (c_filter, b_filter)

def cdf(data):
    '''
    returns a sorted version of data (small to big) and an array of their proportions for a cdf plot
    '''
    sorted_data = np.sort(data)
    # make array of proportions from 0:1, the length of the data
    p = 1. * np.arange(len(data)) / (len(data) - 1)
    return (sorted_data, p)

def select_homeo_units(plot_group, c_filter, b_filter):
    ''' 
    Returns the data table of neurons whose end firing is below (plot_group=1) or above (plot_group=2) 
    double the baseline firing. Used to show plots with only the neurons that homeostased or didn't.
    '''
    low_units = []
    mid_units = []
    high_units = []
    baseline_stop = max(b_filter['time'])
    end_time = max(c_filter['time']) - timedelta(hours = 1)
    lasthr_drug = end_time - timedelta(hours = 1)
    lasthr_baseline = baseline_stop - timedelta(hours = 1)
        
    for unit in c_filter['unit_name'].unique():
        unit_end = c_filter.query('unit_name == @unit and time >= @lasthr_drug and time <= @end_time')
        unit_start = c_filter.query('unit_name == @unit and time >= @lasthr_baseline and time < @baseline_stop')
        mean_end = unit_end['spike_freq'].mean()
        mean_start = unit_start['spike_freq'].mean()
        homeo_val = mean_end/mean_start
        if homeo_val < 0.5:
            low_units.append(unit)
        elif homeo_val <= 1.5:
            mid_units.append(unit)
        else:
            high_units.append(unit)
    if plot_group == 1:
        # select neurons that died
        c_filter = c_filter.loc[c_filter['unit_name'].isin(low_units)]
        b_filter = b_filter.loc[b_filter['unit_name'].isin(low_units)]
    elif plot_group == 2:
        # select neurons that worked
        c_filter = c_filter.loc[c_filter['unit_name'].isin(mid_units)]
        b_filter = b_filter.loc[b_filter['unit_name'].isin(mid_units)]
    elif plot_group == 3:
        # select neurons that never decreased their firing
        c_filter = c_filter.loc[c_filter['unit_name'].isin(high_units)]
        b_filter = b_filter.loc[b_filter['unit_name'].isin(high_units)]
    return (c_filter, b_filter)
    
def calc_end_vs_start(units, baseline_stop, cat_table, end_time = 0):
    '''
    Calculates the ratio of FR at the end of the experiment to pre-drug FR.
    '''
    if end_time == 0:
        end_time = max(cat_table['time'])
    lasthr_drug = end_time - timedelta(hours = 1)
    lasthr_baseline = baseline_stop - timedelta(hours = 1)
    
    units_table = cat_table.loc[cat_table['unit_name'].isin(units)]
    units_end = units_table.query('time >= @lasthr_drug')
    units_start = units_table.query('time >= @lasthr_baseline and time < @baseline_stop')
    mean_end = units_end.groupby(('unit_name'))['spike_freq'].mean()
    mean_end = mean_end.rename('spike_freq').reset_index()
    mean_start = units_start.groupby(('unit_name'))['spike_freq'].mean()
    mean_start = mean_start.rename('spike_freq').reset_index()
    ratio = mean_end['spike_freq']/mean_start['spike_freq']
    return pd.DataFrame(data = {'unit_name': mean_end['unit_name'], 'ratio': ratio})