function nonlinear_energy = neo(input_signal)
%% nonlinear_energy = neo(input_signal)
%
% Applies the nonlinear energy operator to the input signal, defined as
% 
% neo(x[n]) = x[n]^2 - x[n+1]x[n-1]

nonlinear_energy = input_signal(2:end-1).^2 - input_signal(3:end).*input_signal(1:end-2);
