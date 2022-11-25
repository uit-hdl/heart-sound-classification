function acf = getACFschmidt(audio_data,Fs)
% convenience function

% This script takes audio in audio_data with corresponding sampling
% frequency Fs and produces the auto correlation function, after first
% filtering and spike removal. I have made no changes to this code. I
% turned it into a function only for my own convenience, so that I could
% more easily study the part of the algorithm that performs heart rate
% estimation based on the ACF.
%% 25-400Hz 4th order Butterworth band pass
audio_data = butterworth_low_pass_filter(audio_data,2,400,Fs, false);
audio_data = butterworth_high_pass_filter(audio_data,2,25,Fs);

%% Spike removal from the original paper:
audio_data = schmidt_spike_removal(audio_data,Fs);

%% Find the homomorphic envelope
homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(audio_data, Fs);

%% Find the autocorrelation:
y=homomorphic_envelope-mean(homomorphic_envelope);
[c] = xcorr(y,'coeff');
acf = c(length(homomorphic_envelope)+1:end);
end