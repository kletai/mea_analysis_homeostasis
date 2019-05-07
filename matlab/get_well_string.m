function well_string = get_well_string(well_row, well_col)
%%  well_string = get_well_string(well_row, well_col)
%
% Returns the string corresponding to the given well_row and well_col
%  (where the top-left well is row 1, col 1, per the Axis specification)
    well_string(1) = char(well_row + 64); % Offsets 1 to 'A'
    well_string(2) = char(well_col + 48); % Offsets 1 to '1'
