from datetime import datetime
from datetime import timedelta
import pandas as pd

def datetime_str_to_datetime(datetime_str):
    """
    This converts the strings generated when matlab datetimes are written to a table to python datetime objects
    """
    if len(datetime_str) == 24: # Check for milliseconds
        return datetime.strptime(datetime_str, '%d-%b-%Y %H:%M:%S.%f')
    else:
        return datetime.strptime(datetime_str, '%d-%b-%Y %H:%M:%S')

def well_row_letter_to_number(letter):
    """
    Converts the letter corresponding to the Axis well row to the corresponding integer row number.
    Starts with A = 1
    """
    return ord(letter) - 64

def get_row_number(unit_name):
    """
    Returns the row_number corresponding to the row_number of the 
        unit specified by unit_name. Useful for filtering rows by condition.
    """
    row_num = well_row_letter_to_number(unit_name[0])
    return row_num

def get_col_number(unit_name):
    """
    Returns the col_number tuple corresponding to the column of the 
        unit specified by unit_name. Useful for filtering rows/columns by condition.
    """
    col_num = int(unit_name[1])
    return col_num

def get_row_col_number_tuple(unit_name):
    """
    Returns the (row_number, col_number) tuple corresponding to the row_number and column of the 
        unit specified by unit_name. Useful for filtering rows/columns by condition.
    """
    row_num = get_row_number(unit_name)
    col_num = get_col_number(unit_name)
    return (row_num, col_num)

def map_classic_to_lumos(cat_table, map_path, dest_path):
    """
    Converts Axion's electrode mapping system for opaque "Classic" plates to "Lumos" plates
    """
    mapping = pd.read_csv(map_path)
    for orig_str, new_str in zip(mapping['CytoView Well / Electrode'], mapping['48 Well Classic Well / Electrode']):
        cat_table['unit_name'] = cat_table['unit_name'].str.replace(orig_str, new_str)
    cat_table['unit_name'] = cat_table['unit_name'].str.replace("_", "")
    cat_table.to_csv(dest_path)
    return cat_table

def remapped_str_to_datetime(datetime_str):
    """
    This converts the strings generated when a remapped cat_table is written to a table to python datetime objects
    """
    return datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')

def create_stim_starts(stim_times):
    """
    Adds the date and seconds columns of the stim_times table, and returns a datetime series of the
    stimulation tags
    """
    stim_times['date_time'] = stim_times['date_time'].map(datetime_str_to_datetime)
    td=stim_times['seconds'].map(lambda x: timedelta(seconds = x))
    return stim_times['date_time']+td

def get_well_number(unit_name):
    """
    Returns the well number corresponding to the unit specified by unit_name. "A1" is well 1, "B1" is well 2, "A2" is well 7,
        etc.
    """
    (row, col) = get_row_col_number_tuple(unit_name)
    return ((col - 1)*6 + row)

def get_electrode_number(unit_name):
    """
    Returns the electrode number corresponding to the unit specified by unit_name. "A111" is electrode 1, "A121" is electrode 2, 
        "A112" is electrode 5, "B111" is electrode 17, "A211" is electrode 97, etc.
    """
    well = get_well_number(unit_name)
    return ((well-1)*16 + (int(unit_name[2])-1)*4 + int(unit_name[3]))

def get_ele_row_col_number_tuple(unit_name):
    """
    Returns the (row_number, col_number) tuple corresponding to the row_number and column of the 
        electrode specified by unit_name.
    """
    row = int(unit_name[3])
    col = int(unit_name[2])
    return (row, col)