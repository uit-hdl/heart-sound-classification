function AmpFcn = pulseFcn2ndDeg(c,w,h,A,T,L,Hb,x0)
% Creates a pulse function using B-splines, computed over a given number of
% periods.
% input:
% c  = vector with starting values of the pulses.
% w  = vector with the widths of the pulses.
% A  = amplitudes of the pulses.
% h  = vector of values that determines the sharpness of the peaks.
% L  = number of time steps to compute for
% x0 = initial offset
% output:
% vector with the amplitude function value at [x0,x0+1,...,floor(m*T)]


% DEFAULT VALUES
if nargin==6
    Hb = 1;
    x0 = 1;
elseif nargin==7
    x0 = 1;
end
    
% vector of values for which pulse function is to be evaluated
x = x0:floor(L+(x0-1));

n_pulse = length(c);
d = w/2;   % distance between the knots (before shifting by h)
A = A/0.6; % amplitudes

% CREATE KNOTS FOR EACH PULSE
knotsM = zeros(n_pulse,4);
for i=1:n_pulse
    knots = c(i) + (0:d(i):3*d(i)); 
    knots(2:3) = [knots(2)-h(i),knots(3)+h(i)];
    knotsM(i,:) = knots;
end


% SUM PULSE CONTRIBUTIONS TOGETHER
pulseSum = 0;
for i=1:n_pulse
    pulseSum = pulseSum + A(i)*splineEval(knotsM(i,:),mod(x,T));
end
AmpFcn = Hb + pulseSum;

end