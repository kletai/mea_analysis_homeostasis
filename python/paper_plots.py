import pandas as pd
import itertools as it
import seaborn as sns
import numpy as np
from pymea import matlab_compatibility as mc
from matplotlib import pyplot as plt
from matplotlib import mlab as mlab
import random
from datetime import datetime, timedelta
from pymea import supplement_to_plotting as psupp
import math
from scipy import stats

def plot_exp_traces_plus_means(all_table, yscale = 'linear', data_col = 'spike_freq', x_label = 'Time (days)', title = 'Experiment Traces plus Mean', ymax=10, end_time = 1.375, norm_dmso = False, dmso_table = pd.DataFrame([]), c = 'multi', **plot_kwargs):
    """
    Plots spike frequency traces for each experiment in the provided all_table DataFrame, along with 
    the mean trace (in black)
    """
    
    time_vector = all_table['time']
    cat_table_norm = pd.DataFrame([])
    
    for exp in all_table['exp'].unique():
        exp_table = all_table.query('exp == @exp')
        time_vector = exp_table['time']
        if norm_dmso == True:
            exp_dmso_table = dmso_table.query('exp == @exp')
            if exp_dmso_table.empty:
                norm_exp = exp_table[data_col]
            else:
                norm_exp = np.divide(exp_table[data_col], exp_dmso_table[data_col])
            if c == 'multi':
                plt.plot(time_vector, norm_exp, alpha=0.4, **plot_kwargs)
            else:
                plt.plot(time_vector, norm_exp, alpha=0.2, color=c, label='_nolegend_', **plot_kwargs)
            cat_table_norm = pd.concat([cat_table_norm, (pd.DataFrame(data = {'spike_freq': norm_exp, 'time': time_vector, 'exp': exp}))])
            tt_drug, tt_end = exp_stats(cat_table_norm, end_time)
        else:
            if c == 'multi':
                plt.plot(time_vector, exp_table[data_col], alpha=0.4, **plot_kwargs)
            else:
                plt.plot(time_vector, exp_table[data_col], alpha=0.2, color = c, label='_nolegend_', **plot_kwargs)
            tt_drug, tt_end = exp_stats(all_table, end_time)
    
    mean_freq_traces = all_table.groupby(('time'))[data_col].mean()
    mean_freq_traces = mean_freq_traces.rename(data_col).reset_index() # Convert the multiindexed series back to a dataframe
    
    if norm_dmso == True:
        mean_freq_traces_dmso = dmso_table.groupby(('time'))[data_col].mean()
        mean_freq_traces_dmso = mean_freq_traces_dmso.rename(data_col).reset_index()
        if c == 'multi':
            plt.plot(mean_freq_traces['time'], np.divide(mean_freq_traces[data_col], mean_freq_traces_dmso[data_col]), 'k', **plot_kwargs)
        else:
            plt.plot(mean_freq_traces['time'], np.divide(mean_freq_traces[data_col], mean_freq_traces_dmso[data_col]), c, **plot_kwargs)
    else:
        if c == 'multi':
            plt.plot(mean_freq_traces['time'], mean_freq_traces[data_col], 'k', **plot_kwargs)
        else:
            plt.plot(mean_freq_traces['time'], mean_freq_traces[data_col], c, **plot_kwargs)
    
    plt.axhline(1, color='k', label='_nolegend_')
    plt.ylim([0,ymax])
    
    
    print('Drug Stats: ')
    print(tt_drug)
    print('End Stats: ')
    print(tt_end)
    
    plt.yscale(yscale)
    plt.xlabel(x_label)
    plt.ylabel('Fold Induction')
    plt.title(title)
    if c == 'multi':
        plt.legend(all_table['exp'].unique())
    
    return cat_table_norm
    
def exp_stats(all_table, end_time):
    """
    Performs t-test on data comparing fold induction to baseline after drug addition and at the end of the experiment
    """
    baseline = all_table.query('time < 0').groupby('exp').mean()
    drug = all_table.query('time > 0 and time < 0.125').groupby('exp').mean()
    end = all_table.query('time > @end_time').groupby('exp').mean()
    
    stats_table = pd.DataFrame(data = {'baseline': baseline.spike_freq, 'drug': drug.spike_freq, 'end': end.spike_freq})
    
    tt_drug = stats.ttest_rel(stats_table['baseline'], stats_table['drug'])
    tt_end = stats.ttest_rel(stats_table['baseline'], stats_table['end'])
    
    return(tt_drug, tt_end)

def plot_unit_rate_change(cat_table, title='Unit FR Change', norm_dmso = False, dmso_table = pd.DataFrame([]), **plot_kwargs):
    """
    Plots end mean FR vs baseline mean FR for individual units
    """
    unit_freq_mean_base = cat_table.query('time < 0 and time > -0.125').groupby(('exp','unit_name'))['spike_freq'].mean()   
    unit_freq_mean_end = cat_table.query('time < 1.25 and time > 1.125').groupby(('exp','unit_name'))['spike_freq'].mean()
    
    

    for e in cat_table['exp'].unique():
        for unit in cat_table.query('exp == @e')['unit_name'].unique(): 
           # if unit_freq_mean_base.loc[e,unit] > 0 and unit_freq_mean_end.loc[e,unit] > 0:
                plt.plot(unit_freq_mean_base.loc[e,unit], unit_freq_mean_end.loc[e,unit], '.', **plot_kwargs)
        
    max_base = max(unit_freq_mean_base)
    max_end = max(unit_freq_mean_end)
        
    plt.plot([0.001,np.ceil(max(max_base,max_end))],[0.001,np.ceil(max(max_base,max_end))],'k')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([0.001,np.ceil(max(max_base,max_end))])
    plt.ylim([0.001,np.ceil(max(max_base,max_end))])
    plt.axis('equal')
    
    
    plt.ylabel('End mean firing rate (Hz)')
    plt.xlabel('Baseline mean firing rate (Hz)')
    plt.title(title)
    return

def collate_exp(exp_list, file_string, base_start = -0.25, end_time = 1.25):
    all_table = pd.DataFrame()
    u_table = pd.DataFrame()
    for e in exp_list:
        data_path_a = '/home/sean/mea data 2/medians/'+e+'_'+file_string+'_med.csv'
        data_path_u = '/home/sean/mea data 2/units/'+e+'_'+file_string+'_u.csv'
        all_table = all_table.append(pd.read_csv(data_path_a), ignore_index=True)
        u_table = u_table.append(pd.read_csv(data_path_u), ignore_index=True)
    all_table.drop(labels='Unnamed: 0', axis=1,inplace=True)
    u_table.drop(labels='Unnamed: 0', axis=1,inplace=True)

    all_table = all_table.query('time <= @end_time and time >= @base_start')
    all_table.time = all_table.time.round(decimals=2)
    bin_table = all_table.groupby(('exp', 'time'))['spike_freq'].mean()
    bin_table = bin_table.rename('spike_freq').reset_index()

    u_table = u_table.query('time <= @end_time and time >= @base_start')
    u_table.time = u_table.time.round(decimals=2)
    
    return(bin_table,u_table)