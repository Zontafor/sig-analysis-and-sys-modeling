%% hw7_load_data.m

clear; clc; close all;

S      = load('subj_data.mat');   % contains data, deltaT, etc.
data   = S.data;                  % struct with bp, fv, hr, etco2, ...
Ts     = S.deltaT(1);             % sampling interval [s] (0.5)
fs     = 1/Ts;                    % sampling frequency [Hz]

bp     = data.bp(:);              % systemic arterial BP (demeaned)
fv     = data.fv(:);              % MCA CBFV (demeaned)
hr     = data.hr(:);              % heart rate (demeaned)
etco2  = data.etco2(:);           % end-tidal CO2 (demeaned)

N      = numel(bp);               % samples
t      = (0:N-1).' * Ts;          % time axis
maxlag = 100;                     % lag for (cross-)correlations