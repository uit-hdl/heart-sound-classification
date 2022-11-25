function A = fmax2AR2(f0,r)
% Allows user to get the AR(2) polynomial based on desired frequency
% content and sharpness. Maps the desired location of the mode of the
% frequency spectrum, f0, to the characteristic polynomial of the
% corresponding AR(2) process. The input r takes values between 0 and 1,
% and determines how spread out the frequency spectrum of the AR(2) process
% will be, with values close to 1 resulting in a signal with more localized
% frequency.

theta = f0*2*pi;
% r detemines the sharpness of the peak, and the degree to which the signal
% is composed only of that frequency.
P = r*exp(theta*[-1i,1i])';

% translate poles into the corresponding characteristic polynomial of the
% AR(2) process
A = poly(P);
end