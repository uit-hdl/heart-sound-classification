% function [heartRate systolicTimeInterval] = getHeartRateSchmidt(audio_data, Fs, figures)
%
% Derive the heart rate and the sytolic time interval from a PCG recording.
% This is used in the duration-dependant HMM-based segmentation of the PCG
% recording.
%
% This method is based on analysis of the autocorrelation function, and the
% positions of the peaks therein.
%
% This code is derived from the paper:
% S. E. Schmidt et al., "Segmentation of heart sound recordings by a 
% duration-dependent hidden Markov model," Physiol. Meas., vol. 31,
% no. 4, pp. 513-29, Apr. 2010.
%
% Developed by David Springer for comparison purposes in the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
%% INPUTS:
% audio_data: The raw audio data from the PCG recording
% Fs: the sampling frequency of the audio recording
% figures: optional boolean to display figures
%
%% OUTPUTS:
% heartRate: the heart rate of the PCG in beats per minute
% systolicTimeInterval: the duration of systole, as derived from the
% autocorrelation function, in seconds
%
%% Copyright (C) 2016  David Springer
% dave.springer@gmail.com
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [heartRate, systolicTimeInterval,acf,info] = getHeartRateSpringerMod(x, fs, figures)
% MODIFIED VERSION of Schmidts version. Calls the function
% 'findTrueHeartRate', which I wrote to provide a more robust alternative.
% The original version looks witin a fixed search interval for the tallest
% peak, and then derives the heart rate from that.
if nargin < 3
    figures = false;
end

%% Get heatrate:
% From Schmidt:
% "The duration of the heart cycle is estimated as the time from lag zero
% to the highest peaks between 500 and 2000 ms in the resulting
% autocorrelation"
% This is performed after filtering and spike removal:

%% 25-400Hz 4th order Butterworth band pass
x = butterworth_low_pass_filter(x,2,400,fs, false);
x = butterworth_high_pass_filter(x,2,25,fs);

%% Spike removal from the original paper:
x = schmidt_spike_removal(x,fs);

%% Find the homomorphic envelope
homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(x, fs);

%% Find the autocorrelation:
y=homomorphic_envelope-mean(homomorphic_envelope);
[c] = xcorr(y,'coeff');
% we only need the positive axis part of the function:
acf = c(length(homomorphic_envelope)+1:end);

min_index = floor(0.5*fs);
max_index = 2*fs;

% ¤¤¤ MODIFICATION ¤¤¤
%  decision threshold values:
thrHR   = .6;
thrMult = .15;
[heartRate,~,info] = findTrueHeartRate(acf, min_index, max_index, fs, thrHR, thrMult);
info.homomorpic = y;
%% Find the systolic time interval:
% From Schmidt: "The systolic duration is defined as the time from lag zero
% to the highest peak in the interval between 200 ms and half of the heart
% cycle duration"
max_sys_duration = round(((60/heartRate)*fs)/2);
min_sys_duration = round(0.2*fs);

[~, pos] = max(acf(min_sys_duration:max_sys_duration));
systolicTimeInterval = (min_sys_duration+pos)/fs;

% save information about systole fit:
info.y_syst  = acf(min_sys_duration + pos); 
info.sysLeng = systolicTimeInterval;

if(figures)
%     figure('Name', 'Heart rate calculation figure');
    plot((1:numel(acf))/fs,acf);
    hold on;
    true_index = floor(60*fs/heartRate); %#ok<*NASGU>
%     plot(true_index, acf(floor(true_index)),'ro');
    plot((min_sys_duration+pos)/fs, acf(min_sys_duration+pos), 'mo');
    xline([min_index,max_index]/fs)
    xlabel 'time (s)'
%     xlabel('Samples');
%     legend('Autocorrelation', 'Position of max peak used to calculate HR', 'Position of max peak within systolic interval');
end


