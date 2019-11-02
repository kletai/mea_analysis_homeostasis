# mea_analysis_homeostasis
Code for analyzing MEA data in the Gray Lab neuronal firing rate homeostasis assay. Python code written for python 2.7. Matlab pre-processing was run on O2 through Harvard Medical School Research Computing. Code written by Kate Letai, Sam Rendall, and Kelsey Tyssowski.

1.[Installing](#downloading-and-installing-the-code)
  * [Python_mea](#installing-the-python_mea-module)
  * [Matlab](#installing-the-matlab-code)
  
2.[Spike Sorting](#spike-sorting)

3.[Data Visualization](#visualizing-cluster-data)
  * [Exporting Matlab Data to CSV Format](#generating-spike-frequency-traces)
  * [Plotting Data With Python](#python_mea)
  * [Plotting Data with R](#r-functions-used-for-publication)

Downloading and Installing the Code
===================================

You can download this code using git. Open a terminal, navigate to a directory of your chosing and run the following command:
```
git clone https://github.com/kletai/mea_analysis_homeostasis
```
This will create a `mea_analysis_homeostasis` directory in your current directory that contains the latest version of the code. To update to the most recent version of the codebase, run the following from within the `mea_analysis_homeostasis` directory:
```
git pull origin master
```
Important Note: In order to run the MATLAB code in this pipeline, you will also need the "Axion MATLAB Files" package distributed by Axion Biosystems with their AxIS software.

Installing the python_mea module
---------------------------
Code for generating plots using python, pandas, matplotlib and seaborn is included in the pymea module. Using this module requires several dependencies. We recommend using the [Anaconda](https://www.continuum.io/downloads) python distribution, which includes a curated list of commonly used scientific python packages. If you don't want to use anaconda, you can download the rest of the dependencies with `pip`:
```
pip install matplotlib numpy scipy pandas seaborn
```

Once you have all of the dependencies installed, you have to add the python_mea directory to your `PYTHONPATH`. To do this, you can run the following command:
```
export PYTHONPATH="/path/to/mea_analysis_homeostasis:$PYTHONPATH"
```
This has to be done each time you open up a new terminal. However, you can have the command run automatically by adding it to your `bashrc` like this:
```
echo 'export PYTHONPATH="/path/to/mea_analysis_homeostasis:$PYTHONPATH"' >> ~/.bashrc (use ~/.bash_profile for Mac OSX)
```
Once you've done that, you should be able to import the python_mea module into python after opening a new terminal.

Installing the Matlab code
--------------------------
All of the matlab dependencies required by this package are included with the repo. All you need to do is add them to your path. The easiest way to do this is to use the buttons on the matlab console. Open matlab, and select the "set path" option.
Choose the "Add with Subfolders" option.
Navigate to and select the `matlab` folder within the `mea_analysis_homeostasis` directory.
Click "save" to save the changes to your Matlab path.
And you're done! You should be able to run the matlab code in this repo if you followed the instructions correctly!

Spike Sorting
=============
Spike sorting consists of two stages, an automated preprocessing step involving feature extraction and clustering, and a manual review step where features and clusters are visualized and the number of clusters is specified. The entire process is carried out in Matlab.

Preprocessing
-------------
Preprocessing is carried out by the `process_spk_files_parallel.m` matlab function. Our preprocessing was all performed on the O2 High Performance Compute Cluster at Harvard Medical School. The function can be called from matlab by specifying the amount of memory required for the job as the first argument, the amount of time required as the second argument, and finally a cell array containing the paths to the .spk files in a recording as the last argument (.spk files should be given in alphabetical order):
```
process_spk_files_parallel(num_GB, reqTime, {'/path/to/spk_file_1.spk', '/path/to/spk_file_2.spk'}}
```

It is usually inconvenient to do this from the matlab command line however, so we created a script - `processing_launcher_O2.sh` - that can be used to construct this matlab command, and launch the matlab engine with this command as input. This script simply takes in a bash array of spk files. It can be called like this:
```
/path/to/processing_launcher_O2.sh "#GB" "hr:mn:sc" /path/to/spk_file_1.spk /path/to/spk_file_2.spk
```

or if you have a large number of spike files in a directory, it is convenient to use a `find` command to automatically locate all of the .spk files in a directory like this:
```
O2~$ sbatch -p priority -t 5:00 "#GB" "hr:mn:sc" /path/to/processing_launcher_O2.sh $(find /path/to/spk_file_directory -name '*.spk' | sort)
```

Initiating preprocessing using `processing_launcher_O2.sh` will automatically name the output file after the first .spk file specified as an input, with the '.spk' suffix replaced with '.mat'.

How It Works
------------
The spikes on each electrode are clustered based on their first three principal components using Gaussian mixture models. For each electrode, data is attempted to be clustered using up to five clusters. If clustering fails before five clusterings are generated, then clustering is aborted. 

Manual Cluster Numbering
------------------------
After preprocessing completes, results can be visualized and the number of clusters can be determined using a matlab GUI: `spike_analyzer.m`. To use the gui, simply enter the command on the matlab command line:
```
spike_analyzer
```

A prompt should appear instructing you to select the .mat file generated by preprocessing. After you have selected the correct .mat file, a second prompt will appear asking you to select all of the .spk files. Hold shift or ctrl to select all of the files simultaneously, then click 'open'. A GUI will appear displaying the waveforms in 3D feature space, as well as the average waveform for each cluster. A list of commands is in the lower right-hand side of the user-interface. To specify the number of clusters, press a number 1-5. If it appears that the most accurate clustering would require merging two or more automatic clusters, press 'm' to merge clusters. Press 'x' to move to the next electrode. Individual traces can be visualized by selecting the 't' option. However, when in trace mode, data must be loaded from .spk files for each electrode, which is very slow. We reccomend staying in average trace mode ('v') and only using individual traces when absolutely necessary.

Visualizing Cluster Data
========================
Clustered data is located in the .mat file output by `process_spk_files_parallel.m`. The data contained there can be loaded into matlab with the load command:
```
load /path/to/processing_output_file.mat
```

This will create four variables in your matlab environment. `recording_start_time` is a matlab `datetime` array corresponding to the date and time of the recording start time for each .spk file. `final_spike_time` is a `datetime` object corresponding to the date and time of the final spike in the given .spk files. This is used as a proxy for the recording's end time, since that is not stored in any of Axis' output files. `stim_times` is a matlab `datetime` array of the stimulation tags recorded in the experiment, if any. `electrode_containers` is an [num_well_rows, num_well_cols, num_electrode_cols, num_electrode_rows] array of`ElectrodeContainer` objects that store preprocessing data. Each `ElectrodeContainer` stores data for a single maestro electrode. See `ElectrodeContainer.m` or type `help ElectrodeContainer` into the Matlab command prompt for a description of the data contained in each `ElectrodeContainer`.

Generating Spike Frequency Traces
---------------------------------
Although data can be loaded from .mat files and plotted, we have written matlab functions for converting the data into a more manageable form. `generate_spike_frequency_table.m` generates spike frequency timecourses from the data in a given .mat file, and outputs a .csv file containing the frequency timecourses and a time axis. You can use this function like this:
```
generate_spike_frequency_table('/path/to/input_file.mat', '/path/to/output_file.csv')
```
Spike frequencies are output in units of spikes/time bin. An optional `bin_size` parameter can be specified to set the size of time bins in seconds. By default, spikes are binned into 300 second (5 min) time bins. For example, you can set the time bins to be one minute long like this:
```
generate_spike_frequency_table('/path/to/input_file.mat', '/path/to/output_file.csv', 'bin_size', 60)
```
Creating the spike frequency table can take a long time and use a lot of memory, so we recommend using a computing cluster if available.

Memory and time needed will vary according to the size of the .mat file.

Python_mea
-----
We wrote a python module that can be used to quickly create some commonly used plots.\
\
Here is an index of what's in the module:
### matlab_compatibility.py
`datetime_str_to_datetime` - This converts the strings generated when matlab datetimes are written to a table to python datetime objects\
`well_row_letter_to_number` - Converts the letter corresponding to the Axis well row to the corresponding integer row number. Starts with A = 1\
`get_row_number` - Returns the row_number corresponding to the row_number of the unit specified by unit_name. Useful for filtering rows by condition.\
`get_col_number` - Returns the col_number tuple corresponding to the column of the unit specified by unit_name. Useful for filtering rows/columns by condition.\
`get_row_col_number_tuple` - Returns the (row_number, col_number) tuple corresponding to the row_number and column of the unit specified by unit_name. Useful for filtering rows/columns by condition.\
`map_classic_to_lumos` - Converts the electrode mapping of a category dataframe from CytoView Well (clear electrode bottom) to 48 Well Classic (opaque electrode bottom)\
`remapped_str_to_datetime` - Converts the strings generated when a remapped cat_table is written to a table to python datetime objects\
`create_stim_starts` - Adds the date and seconds columns of the stim_times table, and returns a datetime series of the stimulation tags\
`get_well_number` - Returns the well number corresponding to the unit specified by unit_name. "A1" is well 1, "B1" is well 2, "A2" is well 7, etc.\
`get_electrode_number` - Returns the electrode number corresponding to the unit specified by unit_name. "A111" is electrode 1, "A121" is electrode 2, "A112" is electrode 5, "B111" is electrode 17, "A211" is electrode 97, etc.\
`get_ele_row_col_number_tuple` - Returns the row and column of an electrode corresponding to a given unit.\

### plotting.py
`smooth` - Applies a moving average to an array. Useful for smoothing spike frequency traces\
`construct_categorized_dataframe` - Creates a dataframe that includes category information for unit spike timecourses. Uses `filter_unit_columns` as a supporting function.\
`foldInductionPlusMean` - This function plots baseline-normalized plots for a given condition that include both all of the channels passing a filters and all the mean(black)+median(red) of those channels. Uses `foldInductionPlusMean_stim`, `foldInductionPlusMean_ctrl`, `get_mean_med_traces`, and `make_fold_plot`as supporting functions
 - Supporting functions for `foldInductionPlusMean`

### supplement_to_plotting.py
`filter_neurons_homeostasis` - Returns a category dataframe only including neurons that pass the specified filters. Filters first based on baseline firing rate, choosing only neurons whose mean firing rate remains within 0.001-100,000 Hz (default) and who are not silent for more than 8 hours during the experiment. If "ind_filter" is true, neurons are eliminated if their firing rate does not increase upon PTX stimulation. Finally, entire wells are removed if there are fewer than three neurons, or if the median firing rate does not increase by at least 1.3 fold after PTX stimulation.

R functions used for publication
-----
`plotting_for_publication.R` - Contains R functions used to generate plots that are in the publication.
