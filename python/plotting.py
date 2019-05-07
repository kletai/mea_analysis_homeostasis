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

def plot_units_from_spike_table(spike_table):
    time_vector = spike_table['time'].map(mc.datetime_str_to_datetime)
    unit_table = spike_table.copy()
    del unit_table['time']
    num_units = len(unit_table.columns)
    #plt.figure(figsize=(10, 0.1 * num_units))
    for i, unit_name in enumerate(unit_table.columns):
        #plt.subplot(num_units, 1, i + 1)
        plt.figure()
        plot_unit(time_vector, unit_table[unit_name])
        plt.xlabel(unit_name)

def smooth(A, kernel_size=5, mode='same'):
    """
    Computes the moving average of A using a kernel_size kernel.
    """
    kernel = np.ones(kernel_size)/kernel_size
    return np.convolve(A, kernel, mode=mode)

def plot_unit_traces(category_dataframe, yscale = 'linear', **plot_kwargs):
    """
    Plots spike frequency unit traces for each neural unit in the provided category dataframe
    """
    for unit in category_dataframe['unit_name'].unique():
        unit_table = category_dataframe.query('unit_name == @unit')
        plt.plot(unit_table['time'], unit_table['spike_freq'], **plot_kwargs)
        plt.yscale(yscale)

def plot_unit_traces_plus_means(category_dataframe, yscale = 'linear', repeated = False, data_col = 'spike_freq', alt_x = False, x_label = 'Time (days)', title = 'Unit Traces and Mean', **plot_kwargs):
    """
    Plots spike frequency unit traces for each neural unit in the provided category dataframe, along with 
    the mean trace (in black)
    """
    time_days = (category_dataframe['time']-category_dataframe['time'].iloc[0]).map(lambda x: x.days)
    time_seconds = (category_dataframe['time']-category_dataframe['time'].iloc[0]).map(lambda x: x.seconds)
    time_vector = (time_days + (time_seconds/3600/24)).unique()
    
    for unit in category_dataframe['unit_name'].unique():
        unit_table = category_dataframe.query('unit_name == @unit')
        if repeated == True:
            time_vector = unit_table['time']
        plt.plot(time_vector, unit_table[data_col], alpha=0.4, **plot_kwargs)
    
    mean_freq_traces = category_dataframe.groupby(('condition', 'time'))[data_col].mean()
    mean_freq_traces = mean_freq_traces.rename(data_col).reset_index() # Convert the multiindexed series back to a dataframe
    for condition in mean_freq_traces['condition'].unique():
        condition_trace = mean_freq_traces.query('condition == @condition')
        if repeated == True:
            time_vector = condition_trace['time']
        plt.plot(time_vector, condition_trace[data_col], 'k')

    plt.yscale(yscale)
    plt.xlabel(x_label)
    plt.ylabel(data_col)
    plt.title(title)
    plt.legend(mean_freq_traces['condition'].unique()) 

def plot_unit_traces_plus_medians(category_dataframe, yscale = 'linear', data_col = 'spike_freq', alt_x = False, x_label = 'time', title = 'Spike Frequency Traces', **plot_kwargs):
    """
    Plots spike frequency unit traces for each neural unit in the provided category dataframe, along with 
    the mean trace (in black)
    """
    time_days = (category_dataframe['time']-category_dataframe['time'].iloc[0]).map(lambda x: x.days)
    time_seconds = (category_dataframe['time']-category_dataframe['time'].iloc[0]).map(lambda x: x.seconds)
    time_vector = (time_days + (time_seconds/3600/24)).unique()
    
    for unit in category_dataframe['unit_name'].unique():
        unit_table = category_dataframe.query('unit_name == @unit')
        plt.plot(time_vector, unit_table[data_col], **plot_kwargs)
    
    mean_freq_traces = category_dataframe.groupby(('condition', 'time'))[data_col].median()
    mean_freq_traces = mean_freq_traces.rename(data_col).reset_index() # Convert the multiindexed series back to a dataframe
    for condition in mean_freq_traces['condition'].unique():
        condition_trace = mean_freq_traces.query('condition == @condition')
        plt.plot(time_vector, condition_trace[data_col], 'k')

    plt.yscale(yscale)
    plt.xlabel(x_label)
    plt.ylabel('spike frequency')
    plt.title(title)
    plt.legend(mean_freq_traces['condition'].unique()) 
    
def plot_unit_points_plus_means(category_dataframe, title, divide_fn, **plot_kwargs):
    """
    Plots spike frequency points for each neural unit in the provided category dataframe, in each section returned by 
    divide_fn. Mean of all units of the same category for each section is also shown.
    """
    color_map = plt.cm.get_cmap('viridis', category_dataframe['condition'].unique().size)
    color_index = 0
    for cond in category_dataframe['condition'].unique():
        cond_table = category_dataframe.query('condition == @cond')
        cond_table = cond_table.reset_index()
        # Get spks/pulse for each neuron in each time period
        unit_table = cond_table.groupby(('unit_name', lambda x: divide_fn(cond_table, x, 'time')))['spike_freq'].mean() 
        unit_table = unit_table.rename('spike frequency').reset_index() # Convert the multiindexed series back to a dataframe
         # Get average spks/pulse for all neurons in each time period
        mean_table = unit_table.groupby('level_1')['spike frequency'].mean()
        mean_table = mean_table.reset_index() # Convert the multiindexed series back to a dataframe
        
        plt.plot(unit_table['level_1'].astype('int')*3, unit_table['spike frequency'], 'o', color = color_map(color_index), markerfacecolor = 'none', alpha = 0.25, label = '_nolegend_')
        plt.plot(mean_table['level_1'].astype('int')*3, mean_table['spike frequency'], 'o', color = color_map(color_index))
        color_index += 1 # Move to next color for next condition
        
    plt.axhline(y=1, xmin=0, xmax=1, hold=None, color='black')
    plt.legend(category_dataframe['condition'].unique())
    plt.xlabel('hours since start of 1st recording')
    plt.ylabel('spikes/pulse')
    plt.title(title)

def average_timecourse_plot(category_dataframe, **kwargs):
    """
    Generates an average timecourse with error bars for each category in category_dataframe
    see construct_categorized_dataframe for details on generateing the category_dataframe
    """
    sns.pointplot(x='time', y='spike_freq', hue='condition', data=category_dataframe, **kwargs)

def avg_timecourse_plot_2(category_dataframe, **kwargs):
    mean_freqs = category_dataframe.groupby(('condition', 'time'))['spike_freq'].mean()
    std_freqs = category_dataframe.groupby(('condition', 'time'))['spike_freq'].std()
    plt.errorbar()

def plot_unit_frequency_distributions(category_dataframe, **kwargs):
    """
    Plots the distribution of mean frequencies for units in each condition
    """
    mean_freqs_by_condition = category_dataframe.groupby(('condition', 'unit_name'))['spike_freq'].mean()
    mean_freqs_by_condition = mean_freqs_by_condition.rename('mean_freq').reset_index()
    for condition in mean_freqs_by_condition['condition']:
        sns.distplot(mean_freqs_by_condition.query('condition == @condition')['mean_freq'].map(np.log), bins=100)

def plot_mean_frequency_traces(category_dataframe, data_col = 'spike_freq', alt_x = False, x_label = 'time', title = 'Mean Traces', Ret = False, **kwargs):
    """
    Plots the mean frequency trace for each condition in category_dataframe
    """
    mean_freq_traces = category_dataframe.groupby(('condition', 'time'))[data_col].mean()
    mean_freq_traces = mean_freq_traces.rename(data_col).reset_index() # Convert the multiindexed series back to a dataframe
    for condition in mean_freq_traces['condition'].unique():
        condition_trace = mean_freq_traces.query('condition == @condition')
        if alt_x == False:
            plt.plot(condition_trace['time'], condition_trace[data_col])
        else:
            plt.plot(alt_x, condition_trace[data_col])

    plt.xlabel(x_label)
    plt.ylabel(data_col)
    plt.title(title)
    plt.legend(mean_freq_traces['condition'].unique())
    if Ret == True:
        return mean_freq_traces       
    
def plot_median_frequency_traces(category_dataframe, yscale = 'linear', quartiles = True, data_col = 'spike_freq', alt_x = False, x_label = 'time', **kwargs):
    """
    Plots the median frequency trace for each condition in category_dataframe
    """
    median_freq_traces = category_dataframe.groupby(('condition', 'time'))[data_col].median()
    Q1 = category_dataframe.groupby(('condition', 'time'))[data_col].quantile(.25)
    Q3 = category_dataframe.groupby(('condition', 'time'))[data_col].quantile(.75)
    median_freq_traces = median_freq_traces.rename(data_col).reset_index() # Convert the multiindexed series back to a dataframe
    Q1 = Q1.rename(data_col).reset_index()
    Q3 = Q3.rename(data_col).reset_index()
    
    for condition in median_freq_traces['condition'].unique():
        condition_trace = median_freq_traces.query('condition == @condition')
        if alt_x == False:
            ct = plt.plot(condition_trace['time'], condition_trace[data_col])
        else:
            ct = plt.plot(alt_x, condition_trace[data_col])
        if quartiles == True:
            Q1_trace = Q1.query('condition == @condition')
            Q3_trace = Q3.query('condition == @condition')
            if alt_x == False:
                plt.plot(Q1_trace['time'], Q1_trace[data_col], '--', color = ct[0].get_color())
                plt.plot(Q3_trace['time'], Q3_trace[data_col], '--', color = ct[0].get_color())
            else:
                plt.plot(alt_x, Q1_trace[data_col], '--', color = ct[0].get_color())
                plt.plot(alt_x, Q3_trace[data_col], '--', color = ct[0].get_color())

    plt.yscale(yscale)
    plt.xlabel(x_label)
    plt.ylabel('spike frequency')
    plt.title('Median Spike Frequency Traces')
    plt.legend(median_freq_traces['condition'].unique())
    
def get_median_unit_traces(category_dataframe):
    """
    Finds the unit traces in category_dataframe with the median average firing rate
    """
    overall_mean_freq = category_dataframe.groupby(('unit_name', 'condition'))['spike_freq'].mean()
    overall_mean_freq = overall_mean_freq.rename('spike_freq').reset_index() # Convert the multiindexed series back to a dataframe
    median_traces = pd.DataFrame()
    for condition in overall_mean_freq['condition'].unique():
        condition_trace = overall_mean_freq.query('condition == @condition')
        n = len(condition_trace['spike_freq'])
        if n%2 == 0:
            sorted_freq = sorted(condition_trace['spike_freq'])
            median_freq = sorted_freq[n//2 - 1]
        else:
            median_freq = np.median(condition_trace['spike_freq'])
        median_unit = condition_trace[condition_trace.spike_freq == median_freq]['unit_name']
        median_unit.reset_index(drop=True, inplace = True)
        median_unit = median_unit.iloc[0]
        median_traces = pd.concat([median_traces, category_dataframe.query('unit_name == @median_unit')])
    return median_traces
    
def plot_median_unit_frequency_traces(category_dataframe, rec_starts, rec_ends, yscale = 'linear', **kwargs):
    """
    Plots the frequency trace of the unit with the median avg spike freq for each condition in category_dataframe
    """
    median_traces = get_median_unit_traces(category_dataframe)
    for condition in median_traces:
        plot_unit_means_per_rec(median_traces.query('condition == @cond'), rec_starts, rec_ends, num_rec, yscale)
        
    plt.ylabel('spike frequency')
    plt.title('Median Unit Spike Frequency Traces')
    plt.legend(overall_mean_freq['condition'].unique())

def plot_unit_means_per_rec(category_dataframe, rec_starts, rec_ends, num_rec, yscale = 'linear', **plot_kwargs):
    """
    Plots the mean firing of each unit per recording session
    """
    mean_unit_freq = pd.DataFrame()
    for index in range(0,num_rec):
        start1 = rec_starts[index]
        end1 = rec_ends[index]
        rec_table = category_dataframe.query('time >= @start1 and time <= @end1')
        rec_mean_unit_freq = rec_table.groupby('unit_name')['spike_freq'].mean()
        num_units = rec_mean_unit_freq.count()
        start_dt = datetime.strptime(start1, "%Y-%m-%d %H:%M:%S").date()
        start_times = pd.Series([start1]*num_units, index = rec_mean_unit_freq.index)
        rec_data = pd.DataFrame({"mean_freq": rec_mean_unit_freq, "start_time": start_times})
        del rec_data.index.name
        rec_data.reset_index()
        rec_data['unit_name'] = rec_mean_unit_freq.index
        mean_unit_freq = pd.concat([mean_unit_freq, rec_data])
        
    for unit in mean_unit_freq['unit_name'].unique():
        date_table = mean_unit_freq.query('unit_name == @unit')
        plt.plot_date(date_table['start_time'], date_table['mean_freq'], '-o')
        
    plt.yscale(yscale)
    plt.xlabel('time')
    plt.ylabel('mean spike frequency')
    plt.title('Mean Spike Frequency Per Recording')
        
def plot_means_per_rec(category_dataframe, rec_starts, rec_ends, num_rec, yscale = 'linear', **plot_kwargs):
    """
    Plots the mean firing of each condition per recording session
    """
    mean_unit_freq = pd.DataFrame()
    for index in range(0,num_rec):
        start1 = rec_starts[index]
        end1 = rec_ends[index]
        rec_table = category_dataframe.query('time >= @start1 and time <= @end1')
        rec_mean_unit_freq = rec_table.groupby('condition')['spike_freq'].mean()
        num_units = rec_mean_unit_freq.count()
        start_dt = datetime.strptime(start1, "%Y-%m-%d %H:%M:%S").date()
        start_times = pd.Series([start1]*num_units, index = rec_mean_unit_freq.index)
        rec_data = pd.DataFrame({"mean_freq": rec_mean_unit_freq, "start_time": start_times})
        del rec_data.index.name
        rec_data.reset_index()
        rec_data['condition'] = rec_mean_unit_freq.index
        mean_unit_freq = pd.concat([mean_unit_freq, rec_data])
        
    for cond in mean_unit_freq['condition'].unique():
        date_table = mean_unit_freq.query('condition == @cond')
        plt.plot_date(date_table['start_time'], date_table['mean_freq'], '-o')
        
    plt.yscale(yscale)
    plt.xlabel('time')
    plt.ylabel('mean spike frequency')
    plt.title('Mean Spike Frequency Per Recording')
    plt.legend(mean_unit_freq['condition'].unique())
        
def plot_medians_per_rec(category_dataframe, rec_starts, rec_ends, num_rec, yscale='linear', **plot_kwargs):
    """
    Plots the median firing of each condition per recording session
    """
    median_unit_freq = pd.DataFrame()
    for index in range(0,num_rec):
        start1 = rec_starts[index]
        end1 = rec_ends[index]
        rec_table = category_dataframe.query('time >= @start1 and time <= @end1')
        rec_median_unit_freq = rec_table.groupby('condition')['spike_freq'].median()
        num_units = rec_median_unit_freq.count()
        start_dt = datetime.strptime(start1, "%Y-%m-%d %H:%M:%S").date()
        start_times = pd.Series([start1]*num_units, index = rec_median_unit_freq.index)
        rec_data = pd.DataFrame({"median_freq": rec_median_unit_freq, "start_time": start_times})
        del rec_data.index.name
        rec_data.reset_index()
        rec_data['condition'] = rec_median_unit_freq.index
        median_unit_freq = pd.concat([median_unit_freq, rec_data])
        
    for cond in median_unit_freq['condition'].unique():
        date_table = median_unit_freq.query('condition == @cond')
        plt.plot_date(date_table['start_time'], date_table['median_freq'], '-o')
        
    plt.yscale(yscale)
    plt.xlabel('time')
    plt.ylabel('median spike frequency')
    plt.title('Median Spike Frequency Per Recording')
    plt.legend(median_unit_freq['condition'].unique())
    
def construct_categorized_dataframe(data_table, filter_dict, var_name = 'spike_freq'):
    """
    Takes the data from the matlab csv generated by preprocessing and applies filters to column names
    allowing for the categorization of data

    data_table - pandas DataFrame - should be populated from the .csv file generated by the 
        "generate_frequency_table.m" matlab script
    filter_dict - dictionary of the form {'condition_name': condition_filter}, where 
        condition_name is a string used to identify an experimental condition, and condition filter
        is a function that returns True for the unit_names corresponding to the desired condition
    """
    time_vector = data_table['time'].map(mc.datetime_str_to_datetime)
    unit_table = data_table.drop('time', axis=1)
    condition_dicts = (
        {
            'time': time_vector,
            'condition': condition_name,
            var_name: condition_column,
            'unit_name': condition_column.name,
            'well': mc.get_well_number(condition_column.name)
        } for condition_name, condition_filter in filter_dict.iteritems()
            for condition_column in filter_unit_columns(condition_filter, unit_table)
    )
    condition_tables = it.imap(pd.DataFrame, condition_dicts)
    return pd.concat(condition_tables)

def construct_categorized_dataframe_burst(data_table, filter_dict):
    """
    Takes the data from the matlab csv generated by preprocessing and applies filters to column names
    allowing for the categorization of data

    data_table - pandas DataFrame - should be populated from the .csv file generated by the 
        "generate_frequency_table.m" matlab script
    filter_dict - dictionary of the form {'condition_name': condition_filter}, where 
        condition_name is a string used to identify an experimental condition, and condition filter
        is a function that returns True for the unit_names corresponding to the desired condition
    get_power - function that returns the power of optical stimulation by mapping the unit name to the 
        well map
    get_width - function that returns the width of each pulse of optical stimulation by mapping the unit 
        name to the well map
    """
    condition_table = pd.DataFrame()
    for condition_name, condition_filter in filter_dict.iteritems():#iterates through each condition
        filtered_table = filter_unit_rows(condition_filter, data_table)
        #print(filtered_table)
        filtered_table['condition'] = condition_name
        condition_table = condition_table.append(filtered_table, ignore_index=True)
        
    return condition_table

def filter_unit_columns(predicate, unit_table):
    """
    Generates columns from unit_table whose names satisfy the condition specified in predicate

    predicate - function that returns true for desired unit names
    unit_table - data_mat containing firing rates over time from each unit, with the time column ommited
    """
    unit_column_names = filter(predicate, unit_table.columns)
    for column_name in unit_column_names:
        yield unit_table[column_name]
        
def filter_unit_rows(predicate, data_table):
    """
    Generates rows from unit_table whose times satisfy the condition specified in predicate

    predicate - function that returns true for desired unit names
    data_table - data_mat containing data over time from each unit, with the time column included
    """    
    data_row_units = filter(predicate, data_table['unit_name'].unique())
    selected_unit_table = pd.DataFrame()
    for row_unit in data_row_units:
        selected_unit_table = pd.concat([selected_unit_table, data_table[data_table['unit_name'] == row_unit]])
        
    return selected_unit_table

def smooth_categorized_dataframe_unit_traces(category_dataframe, kernel_size=5):
    cat_df_copy = category_dataframe.copy()
    for unit_name in cat_df_copy['unit_name'].unique():
        unit_table = cat_df_copy.query('unit_name == @unit_name')
        smooth_trace = smooth(unit_table['spike_freq'], kernel_size=kernel_size, mode='same')
        cat_df_copy.loc[cat_df_copy['unit_name'] == unit_name, 'spike_freq'] = smooth_trace

    return cat_df_copy


def makeTables(b_start, b_stop, s_start, e_start, cat_table):
    '''
    Makes tables of the baseline portion, stimulated portion and the end portion (i.e. the part of the time course that you deem to have adapted) from the table of the whole time course
    '''

    baseline_table = cat_table.query('time < "%s"'%b_stop).query('time > "%s"'%b_start)
    stim_table = cat_table.query('time > "%s"'%s_start)
    end_table = cat_table.query('time > "%s"'%e_start)
    return(baseline_table, stim_table, end_table)

def get_mean_med_traces(c_filter, data_col, b_filter, FR_gradient):
    mean_freq_traces = c_filter.groupby(('condition', 'time'))[data_col].mean()
    mean_freq_traces = mean_freq_traces.rename(data_col).reset_index() # Convert the multiindexed series back to a dataframe
    mean_freq_traces_b = b_filter.groupby(('condition', 'time'))[data_col].mean()
    mean_freq_traces_b = mean_freq_traces_b.rename(data_col).reset_index() # Convert the multiindexed series back to a dataframe
    
    median_freq_traces = c_filter.groupby(('condition', 'time'))[data_col].median()
    median_freq_traces = median_freq_traces.rename(data_col).reset_index() # Convert the multiindexed series back to a dataframe
    median_freq_traces_b = b_filter.groupby(('condition', 'time'))[data_col].median()
    median_freq_traces_b = median_freq_traces_b.rename(data_col).reset_index() # Convert the multiindexed series back to a dataframe
    
    if FR_gradient == True:
        b_mean_freq = b_filter.groupby(('unit_name'))['spike_freq'].mean()
        b_mean_freq = b_mean_freq.rename('spike_freq')#.reset_index()
        u_color = (np.log10(b_mean_freq)+3)/3
    else:
        b_mean_freq = b_filter.groupby(('unit_name'))['spike_freq'].mean()
        b_mean_freq = b_mean_freq.rename('spike_freq')#.reset_index()
        u_color = np.random.random_sample()*b_mean_freq
    return(mean_freq_traces, mean_freq_traces_b, median_freq_traces, median_freq_traces_b, u_color)


def make_fold_plot(c_filter, t_start, u_color, FR_gradient, plotFolds, norm_by_median, norm_by_mean, mean_freq_traces_b, median_freq_traces_b, mean_freq_traces, median_freq_traces, y_scale, data_col, data_col_mm, title, ymax):
    plt.xlabel('Time (days)')
    plt.ylim(0.00005,ymax)
    for unit_name in c_filter['unit_name'].unique():
        unit = c_filter.query('unit_name == @unit_name')
        u_time = unit['time']
        time_vector_u = u_time-t_start
        time_vector_u = time_vector_u.map(lambda x: x.total_seconds()/86400.0)
        this_color = u_color[unit_name]
        if FR_gradient == True:
            if plotFolds == True:
                if norm_by_median.empty == False:
                    plt.plot(time_vector_u, np.divide(unit['folds'], norm_by_mean), color=plt.cm.gnuplot2(this_color, .4))
                else:
                    plt.plot(time_vector_u, unit['folds'], color=plt.cm.gnuplot2(this_color, .4))
            else:
                plt.plot(time_vector_u, unit[data_col], color=plt.cm.gnuplot2(this_color, .4))
            #color_ind = color_ind+1
        else:
            if plotFolds == True:
                if norm_by_median.empty == False:
                    plt.plot(time_vector_u, np.divide(unit['folds'], norm_by_mean), color=(random.random(), random.random(), random.random(), .4))
                else:
                    plt.plot(time_vector_u, unit['folds'], color=(random.random(), random.random(), random.random(), .4))
            else:
                plt.plot(time_vector_u, unit[data_col], color=(random.random(), random.random(), random.random(), .4))
        
    #print(mean_freq_traces_b)
    meanOfMean = np.mean(mean_freq_traces_b[data_col_mm])
    meanOfMedian = np.mean(median_freq_traces_b[data_col_mm])
    m_time = mean_freq_traces['time']
    time_vector_m = m_time-t_start
    time_vector_m = time_vector_m.map(lambda x: x.total_seconds()/86400.0)
    if plotFolds == True:
        plt.axhline(y=1, xmin=0, xmax=1, hold=None, color='black')
        if norm_by_mean.empty == False:
            plt.plot(time_vector_m, np.divide(mean_freq_traces[data_col]/meanOfMean,norm_by_mean), color=(0,0,0))
            plt.plot(time_vector_m, np.divide(median_freq_traces[data_col]/meanOfMedian,norm_by_median), 'r')
        else:
            plt.plot(time_vector_m, mean_freq_traces[data_col_mm]/meanOfMean, color=(0,0,0))
            plt.plot(time_vector_m, median_freq_traces[data_col_mm]/meanOfMedian, 'r')
        plt.ylabel('Fold Induction of Spike Frequency (Hz)')
    else:
        plt.axhline(y=meanOfMean, xmin=0, xmax=1, hold=None, color='black')
        plt.plot(time_vector_m, mean_freq_traces[data_col], color=(0,0,0))
        plt.plot(time_vector_m, median_freq_traces[data_col], 'r')      
        plt.ylabel('Spike Frequency (Hz)')
    plt.yscale(y_scale)
    plt.title(title)
    plt.show()
    return(meanOfMean, meanOfMedian, time_vector_m)

def foldInductionPlusMean_stim(cat_table, baseline_table, stim_table, condition, title, var, minHz, maxHz, ymax, plotFolds, foldMin, y_scale, filter_wells, data_col, data_col_mm, plot_group, FR_gradient, norm_by_mean, norm_by_median, plot_wells):
    '''
    This function plots baseline-normalized plots for a given condition that include both all of the channels passing filters and the mean(black)+median(red) of those channels--use for stimulated samples b/c filters out things that don't change with stim
    '''
    c = cat_table.query('condition == "%s"'%condition)
    b = baseline_table.query('condition == "%s"'%condition)
    s = stim_table.query('condition == "%s"'%condition)
    t_start = min(s['time'])
    
    c_filter, b_filter, count_real, count_live, cf = psupp.filter_neurons_homeostasis(c, b, s, ind_filter=True, var=var, minHz=minHz, maxHz=maxHz, foldMin=foldMin, filter_wells=filter_wells, data_col=data_col)
    
    if c_filter.empty:
        print "No valid units for condition",condition
        print('respond to drug: 0')
        print('stay alive: ' + str(count_live))
        print('real: ' + str(count_real))
        print('condition: ' + str(len(c['unit_name'].unique())))
        return
    
    if plot_group != 0:
        c_filter, b_filter = psupp.select_homeo_units(plot_group, c_filter, b_filter)

    mean_freq_traces, mean_freq_traces_b, median_freq_traces, median_freq_traces_b, u_color = get_mean_med_traces(c_filter, data_col_mm, b_filter, FR_gradient)
        
    meanOfMean, meanOfMedian, time_vector_m = make_fold_plot(c_filter, t_start, u_color, FR_gradient, plotFolds, norm_by_median, norm_by_mean, mean_freq_traces_b, median_freq_traces_b, mean_freq_traces, median_freq_traces, y_scale, data_col, data_col_mm, title, ymax)
    
    #plot individual well plots
    if plot_wells == True:
        for w in c_filter['well'].unique():
            plt.figure()
            well_c = c_filter.query('well == @w')
            well_b = b_filter.query('well == @w')
            well_mft, well_mftb, well_mdft, well_mdftb, well_color = get_mean_med_traces(well_c, data_col_mm, well_b, FR_gradient)
            well_title = 'Well ' + str(w)
            make_fold_plot(well_c, t_start, well_color, FR_gradient, plotFolds, norm_by_median, norm_by_mean, well_mftb, well_mdftb, well_mft, well_mdft, y_scale, data_col, data_col_mm, well_title, ymax)          

    print('respond to drug: ' + str(len(c_filter['unit_name'].unique())))
    print('stay alive: ' + str(count_live))
    print('real: ' + str(count_real))
    print('condition: ' + str(len(c['unit_name'].unique())))    
    
    return (c_filter['unit_name'].unique(), mean_freq_traces[data_col_mm]/meanOfMean, median_freq_traces[data_col_mm]/meanOfMedian, time_vector_m)


def foldInductionPlusMean_ctrl(cat_table, baseline_table, stim_table, condition, title, var, minHz, maxHz, ymax, plotFolds, foldMin, y_scale, filter_wells, data_col, data_col_mm, plot_group, FR_gradient, norm_by_mean, norm_by_median, plot_wells):
    '''
    This function plots baseline-normalized plots for a given condition that include both all of the channels passing filters and the mean(black)+median(red) of those channels--use for unstim samples
    '''
    c = cat_table.query('condition == "%s"'%condition)
    b = baseline_table.query('condition == "%s"'%condition)
    s = stim_table.query('condition == "%s"'%condition)
    t_start = min(s['time'])

    c_filter, b_filter, count_real, count_live, count_final = psupp.filter_neurons_homeostasis(c, b, s, ind_filter=False, var=var, minHz=minHz, maxHz=maxHz, foldMin=foldMin, filter_wells=False, data_col = data_col)
    
    if c_filter.empty:
        print "No valid units for condition",condition
        print('stay alive: ' + str(count_live))
        print('real: ' + str(count_real))
        print('condition: ' + str(len(c['unit_name'].unique())))
        return (0,0,0)

    # select to show only neurons that do homeostase or don't
    if plot_group != 0:
        c_filter, b_filter = psupp.select_homeo_units(plot_group, c_filter, b_filter)
    
    mean_freq_traces, mean_freq_traces_b, median_freq_traces, median_freq_traces_b, u_color = get_mean_med_traces(c_filter, data_col_mm, b_filter, FR_gradient)
        
    meanOfMean, meanOfMedian, time_vector_m = make_fold_plot(c_filter, t_start, u_color, FR_gradient, plotFolds, norm_by_median, norm_by_mean, mean_freq_traces_b, median_freq_traces_b, mean_freq_traces, median_freq_traces, y_scale, data_col, data_col_mm, title, ymax)
    
    #plot individual well plots
    if plot_wells == True:
        for w in c_filter['well'].unique():
            plt.figure()
            well_c = c_filter.query('well == @w')
            well_b = b_filter.query('well == @w')
            well_mft, well_mftb, well_mdft, well_mdftb, well_color = get_mean_med_traces(well_c, data_col_mm, well_b, FR_gradient)
            well_title = 'Well ' + str(w)
            make_fold_plot(well_c, t_start, well_color, FR_gradient, plotFolds, norm_by_median, norm_by_mean, well_mftb, well_mdftb, well_mft, well_mdft, y_scale, data_col, data_col_mm, well_title, ymax)
    
    print('stay alive: ' + str(count_live))
    print('real: ' + str(count_real))
    print('condition: ' + str(len(c['unit_name'].unique())))
    
    plt.show()

    return (c_filter['unit_name'].unique(), mean_freq_traces[data_col_mm]/meanOfMean, median_freq_traces[data_col_mm]/meanOfMedian, time_vector_m)

def foldInductionPlusMean(cat_table, drug_time, condition, title, var=10, minHz = 0.001, maxHz = 100, ind_filter = True, ymax = 10, plotFolds = True, foldMin = 0.001, y_scale = 'linear', filter_wells = False, data_col ='spike_freq', data_col_mm = 'folds', plot_group = 0, FR_gradient = True, norm_by_mean = pd.Series([]), norm_by_median = pd.Series([]), plot_wells=True):
    '''
    Combine stim and ctrl fxns
    '''
    mean = False
    median = False
    baseline_table = cat_table.query('time < @drug_time')
    stim_table = cat_table.query('time >= @drug_time')
    if ind_filter:
        filtered_units, mean, median, time_vector = foldInductionPlusMean_stim(cat_table, baseline_table, stim_table, condition, title, var, minHz, maxHz, ymax, plotFolds, foldMin, y_scale, filter_wells, data_col, data_col_mm, plot_group, FR_gradient, norm_by_mean, norm_by_median, plot_wells)
    else:
        filtered_units, mean, median, time_vector = foldInductionPlusMean_ctrl(cat_table, baseline_table, stim_table, condition, title, var, minHz, maxHz, ymax, plotFolds, foldMin, y_scale, filter_wells, data_col, data_col_mm, plot_group, FR_gradient, norm_by_mean, norm_by_median, plot_wells)
        
    return filtered_units, mean, median, time_vector


def count_active_neurons(cat_table, baseline_table = 0, stim_table = 0, threshold = 0.001, folds = 0, kill_neurons = 0, return_value = 0):
    '''
    Count and plot the number of neurons firing above a threshold at each time point. If folds == 1, a neuron is deemed firing
    if its fold induction is greater than threshold. If kill_neurons == 1, a neuron is deemed dead for the rest of the
    experiment as soon as its fold induction goes below threshold.
    '''
    time_days = (cat_table['time']-cat_table['time'].iloc[0]).map(lambda x: x.days)
    time_seconds = (cat_table['time']-cat_table['time'].iloc[0]).map(lambda x: x.seconds)
    time_vector = (time_days + (time_seconds/3600/24)).unique()
    
    if folds == 0:
        count_table = cat_table
        
    elif folds == 1:
        meanOfBaseline = baseline_table.groupby('unit_name')['spike_freq'].mean()
        meanOfBaseline = meanOfBaseline.reset_index()
        count_table = pd.DataFrame()
        for unit in baseline_table['unit_name'].unique():
            unit_table = cat_table.query('unit_name == @unit')
            unit_mean_b = meanOfBaseline.query('unit_name == @unit')['spike_freq']
            unit_mean_b = unit_mean_b.reset_index()
            b = unit_mean_b.get_value(0,'spike_freq')
            unit_table.loc[:,'spike_freq'] = unit_table['spike_freq']/b
            below_fold_thresh = unit_table.query('spike_freq < @threshold')
            if not below_fold_thresh.empty:
                first_death = min(below_fold_thresh['time'])
                unit_table.set_value((unit_table.loc[unit_table['time'] > first_death]).index, 'spike_freq', 0)
                ccc=unit_table
            count_table = pd.concat([count_table, unit_table])
        
    above_threshold = count_table.query('spike_freq > @threshold')
    time_grouped_counts = above_threshold.groupby(('time'))['unit_name'].count()
    time_grouped_counts = time_grouped_counts.rename('count').reset_index() # Convert the multiindexed series back to a dataframe
    
    plt.plot(time_vector, time_grouped_counts['count'])
    plt.xlabel('time')
    plt.ylabel('Number of active units')
    plt.title(cat_table.get_value(0,'condition')[0])
    if return_value:
        return time_grouped_counts
    
def compare_active_per_recording(cat_table, threshold, rec_starts, rec_ends, num_rec):
    '''
    For each recording session, find the number of new neurons and the 
    number of neurons that have stopped firing
    '''
    above_threshold = cat_table.query('spike_freq > @threshold')
    only_1 = [0]*(num_rec-1);
    only_2 = [0]*(num_rec-1);
    for index in range(0,num_rec-1):
        start1 = rec_starts[index]
        end1 = rec_ends[index]
        start2 = rec_starts[index+1]
        end2 = rec_ends[index+1]
        group_1 = above_threshold.query('time >= @start1 and time <= @end1')
        group_2 = above_threshold.query('time >= @start2 and time <= @end2')
        units_1 = group_1['unit_name'].unique()
        units_2 = group_2['unit_name'].unique()
        both = list(set(units_1) | set(units_2))
        only_1[index] = len(both) - len(units_2) #Count the number of units in group1 but not group2
        only_2[index] = len(both) - len(units_1) 
   
    rec_starts_series = pd.Series(rec_starts)
    recs = rec_starts_series.map(mc.remapped_str_to_datetime)
    plt.plot_date(recs[1:num_rec], only_2, '-', label = "new")
    plt.plot_date(recs[1:num_rec], only_1, '-', label = "died")
    plt.legend()
    plt.xlabel('Recording session')
    plt.ylabel('Number of units')
    plt.title('Neuron turnover')
    
def compare_active_per_sec(cat_table, threshold):
    '''For each recording session, find the number of new neurons and the 
    number of neurons that have stopped firing
    '''
    above_threshold = cat_table.query('spike_freq > @threshold').reset_index()
    secs = above_threshold['time']
    num_sec = len(secs)
    only_1 = [0]*(num_sec-1);
    only_2 = [0]*(num_sec-1);
    for index in range(0,num_sec-1):
        start1 = secs.iloc[index]
        start2 = secs.iloc[index+1]
        group_1 = above_threshold.query('time == @start1')
        group_2 = above_threshold.query('time == @start2')
        units_1 = group_1['unit_name'].unique()
        units_2 = group_2['unit_name'].unique()
        both = list(set(units_1) | set(units_2))
        only_1[index] = len(both) - len(units_2) #Count the number of units in group1 but not group2
        only_2[index] = len(both) - len(units_1) 
   
    plt.plot_date(secs.iloc[0:num_sec-1], only_2, '-', label = "new")
    plt.plot_date(secs.iloc[0:num_sec-1], only_1, '-', label = "died")
    plt.legend()
    plt.xlabel('Recording session')
    plt.ylabel('Number of units')
    plt.title('Neuron turnover')

def unit_mean_freq_hist(category_dataframe, num_bins = 50, plot = 'linear', title = 'Mean Firing Rate Per Unit'):
    '''
    Plots histogram showing the distribution of mean firing rate of each unit in category_dataframe
    '''
    unit_freq_mean = category_dataframe.groupby(('unit_name'))['spike_freq'].mean()
    unit_freq_mean = unit_freq_mean.rename('spike_freq').reset_index() # Convert the multiindexed series back to a dataframe
    unit_freq_mean = unit_freq_mean.query('spike_freq > 0')
    sigma = unit_freq_mean['spike_freq'].std()
    mu = unit_freq_mean['spike_freq'].mean()
    if plot == 'linear':
        n, bins, patches = plt.hist(unit_freq_mean['spike_freq'], bins = num_bins)
    elif plot == 'log':
        n, bins, patches = plt.hist(np.log10(unit_freq_mean['spike_freq']), bins = num_bins)
        
    # add a 'best fit' line
   # y = mlab.normpdf(bins, mu, sigma)
    #plt.plot(bins, y, 'r--')
    #plt.axvline(mu, color ='r')
    plt.title(title)

def unit_mean_freq_hist_compare_cond(category_dataframe, num_bins = 50, plot = 'linear'):
    '''
    Plots histogram showing the distribution of mean firing rate of each unit in category_dataframe
    '''
    for cond in category_dataframe['condition'].unique():
        cond_table = category_dataframe.query('condition == @cond')
        unit_freq_mean = cond_table.groupby(('unit_name'))['spike_freq'].mean()
        unit_freq_mean = unit_freq_mean.rename('spike_freq').reset_index()
        unit_freq_mean = unit_freq_mean.query('spike_freq > 0')
        sigma = unit_freq_mean['spike_freq'].std()
        mu = unit_freq_mean['spike_freq'].mean()
        if plot == 'linear':
            n, bins, patches = plt.hist(unit_freq_mean['spike_freq'], bins = num_bins)
        elif plot == 'log':
            n, bins, patches = plt.hist(np.log10(unit_freq_mean['spike_freq']), bins = num_bins)
        
        # add a 'best fit' line
        y = mlab.normpdf(bins, mu, sigma)
        plt.plot(bins, y, 'r--')
        plt.axvline(mu, color ='r')
        plt.title('Mean Firing Rate per Unit')

def unit_mean_freq_bar(category_dataframe):
    '''
    Plots histogram showing the distribution of mean firing rate of each unit in category_dataframe
    '''
    num_cond = len(category_dataframe['condition'].unique())
    i=0
    for cond in category_dataframe['condition'].unique():
        cond_table = category_dataframe.query('condition == @cond')
        unit_freq_mean = cond_table.groupby(('unit_name'))['spike_freq'].mean()
        unit_freq_mean = unit_freq_mean.rename('spike_freq').reset_index() # Convert the multiindexed series back to a dataframe
        unit_freq_mean = unit_freq_mean.query('spike_freq > 0')
        sigma = unit_freq_mean['spike_freq'].std()
        mu = unit_freq_mean['spike_freq'].mean()
        plt.bar(i,mu, yerr=sigma)
        i=i+1

    plt.xticks([0,1,2,3], category_dataframe['condition'].unique())
    plt.title('Mean Firing Rate per Unit')
    plt.ylabel('Firing Rate (spk/s)')

def neurons_per_well(cat_table):
    '''
    Plots a bar plot of the number of active neurons per well
    '''
    units = cat_table['unit_name'].unique()
    wells = np.zeros(48)
    ind = range(1,49)
    for unit in units:
        well = mc.get_well_number(unit)
        wells[well-1] = wells[well-1]+1
    barlist = plt.bar(ind,wells)
    for i in range(0,6):
        barlist[i].set_color('r')
    for i in range(6,12):
        barlist[i].set_color('b')
    for i in range(12,18):
        barlist[i].set_color('g')
    for i in range(18,24):
        barlist[i].set_color('y')
    for i in range(24,30):
        barlist[i].set_color('c')
    for i in range(30,36):
        barlist[i].set_color('k')
    for i in range(36,42):
        barlist[i].set_color('m')
    plt.title('Neurons per Well')
    plt.xlabel('Well Number')
    plt.ylabel('Number of Neurons')

def neurons_per_electrode(cat_table):
    '''
    Plots a bar plot of the number of active neurons per well
    '''
    units = cat_table['unit_name'].unique()
    eles = np.zeros(768)
    ind = range(1,769)
    for unit in units:
        ele = mc.get_electrode_number(unit)
        eles[ele-1] = eles[ele-1]+1
    barlist = plt.bar(ind,eles)
    for i in range(0,96):
        barlist[i].set_color('r')
    for i in range(96,192):
        barlist[i].set_color('b')
    for i in range(192,288):
        barlist[i].set_color('g')
    for i in range(288,384):
        barlist[i].set_color('y')
    for i in range(384,480):
        barlist[i].set_color('c')
    for i in range(480,576):
        barlist[i].set_color('k')
    for i in range(576,672):
        barlist[i].set_color('m')
    plt.title('Neurons per Electrode')
    plt.xlabel('Electrode Number')
    plt.ylabel('Number of Neurons')

def heatmap_active_wells(unit_names):
    '''
    Makes a heatmap comparing the number of active neurons per well
    '''
    wells = np.zeros((6,8))
    for unit in unit_names:
        row, col = mc.get_row_col_number_tuple(unit)
        wells[row-1,col-1] = wells[row-1,col-1]+1
    plt.imshow(wells)
    plt.colorbar()

def heatmap_active_electrodes(unit_names):
    '''
    Makes a heatmap comparing the number of active neurons per electrode
    '''
    electrodes = np.zeros((24, 32))
    for unit in unit_names:
        well_row, well_col = mc.get_row_col_number_tuple(unit)
        ele = mc.get_electrode_number(unit)
        row_in_well, col_in_well = mc.get_ele_row_col_number_tuple(unit)
        ele_row = 4*(well_row-1) + row_in_well
        ele_col = 4*(well_col-1) + col_in_well
        ele_row = int(ele_row)
        electrodes[ele_row-1, ele_col-1] = electrodes[ele_row-1, ele_col-1]+1
    plt.imshow(electrodes)
    plt.colorbar()
    for xline in np.arange(-0.5,32,4):
        plt.axvline(x = xline)
    for yline in np.arange(-0.5,24,4):
        plt.axhline(y=yline)
        
def cdf_foldInduction(b_filter, s_filter, title = ""):
    '''
    Plots the cumulative distribution function of the fold induction post baseline at various 
    timepoints during a homeostasis experiment
    '''
    s_start = s_filter['time'].iloc[0]
    hours = np.array([0, 1, 3, 6, 12, 24, 36, 48,])# 72, 96, 120, 144, 168, 192])
    max_hours = (max(s_filter['time']) - s_start).days*24 + (max(s_filter['time']) - s_start).seconds/3600
    hours = hours[hours<= max_hours]
    
    color_idx = np.linspace(0.1, 1, 1+len(hours))
    stim_color_ind = 1
    
    baseline_fold = b_filter.groupby(('unit_name'))['folds'].mean()
    baseline_fold = baseline_fold.rename('folds').reset_index()
    sort, p = psupp.cdf(baseline_fold['folds'])
    plt.plot(sort, p, color='r')#plt.cm.gist_yarg(0.07))
    
    for timepoint in hours:
        period_start = s_start + timedelta(hours=timepoint)
        period_stop = period_start + timedelta(hours=1)
        s_period = s_filter.query('time > @period_start and time < @period_stop')
        period_fold = s_period.groupby('unit_name')['folds'].mean()
        period_fold = period_fold.rename('folds').reset_index()
        sort, p = psupp.cdf(period_fold['folds'])
        plt.plot(sort, p, color=plt.cm.gist_yarg(color_idx[stim_color_ind]))
        stim_color_ind = stim_color_ind+1
        
    legend_labels = np.append("Baseline", hours)
    plt.legend(legend_labels)
    plt.xlabel('FR/baseline')
    plt.ylabel('Fraction of Population')
    plt.xscale('log')
    plt.title('CDF ' + title)
    
    plt.xlim([0, 10])
    plt.axhline(y=0.5, linestyle='--', color = 'k')
    
def hist_end_vs_start(units, baseline_stop, cat_table, cond, end_time = 0, nbins=100):
    '''
    Plots a histogram of the firing rate (of each unit in units) during the last hour divided by the hour before drug was added.
    '''
    ratio_table = psupp.calc_end_vs_start(units, baseline_stop, cat_table, end_time)
    hist_val = ratio_table['ratio']
    
    plt.hist(hist_val, nbins)#, bins=np.logspace(np.log10(0.0001),np.log10(1000), nbins))
    plt.xscale('linear')
    plt.axvline(1, color='k')
    plt.ylabel('Frequency')
    plt.xlabel('FR end / FR start')
    plt.title(cond)

def scatter_homeo_vs_baseline(units, baseline_stop, cat_table, baseline_table, cond, end_time = 0, nbins=100):
    '''
    Makes a scatter plot of the end/start ratio on the y-axis, and the average baseline FR on the x-axis, for the units in 'units'.
    '''
    # Calculate end/start ratio
    ratio_table = psupp.calc_end_vs_start(units, baseline_stop, cat_table, end_time)
    ratio_table = ratio_table.sort_values(by = 'unit_name')
    # Calculate mean FR
    baseline_table = baseline_table.loc[baseline_table['unit_name'].isin(units)]
    unit_freq_mean = baseline_table.groupby(('unit_name'))['spike_freq'].mean()
    unit_freq_mean = unit_freq_mean.rename('spike_freq').reset_index() # Convert the multiindexed series back to a dataframe
    unit_freq_mean = unit_freq_mean.sort_values(by = 'unit_name')
    
    #joined = unit_freq_mean.set_index('unit_name').join(ratio_table.set_index('unit_name'), on = 'unit_name')
    #caller.join(other.set_index('key'), on='key')
    
    plt.scatter(unit_freq_mean['spike_freq'], ratio_table['ratio'])
    plt.xlabel('Mean Baseline FR (Hz)')
    plt.ylabel('End/Start')
    plt.title(cond)
    #return (ratio_table, unit_freq_mean)#, joined)