function electrode_string = get_electrode_string(electrode_col, electrode_row)
%%  electrode_string = get_electrode_string(electrode_col, electrode_row)
%
% returns the string corresponding to the given electrode column and row.
%  NOTE that Axion decided to specify electrodes as [column, row]
    electrode_string(1) = char(electrode_col + 48); % Offsets 1 to '1'
    electrode_string(2) = char(electrode_row + 48);
