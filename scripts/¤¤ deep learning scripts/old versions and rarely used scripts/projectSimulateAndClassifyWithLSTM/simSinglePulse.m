function y = simSinglePulse(c,w,h,A,T,L,Hb,x0)
% DEFAULT VALUES
if nargin==6
    Hb = 1;
    x0 = 1;
elseif nargin==7
    x0 = 1;
end

ampFcn  = pulseFcn2ndDeg(c,w,h,A(1,:),T,L,Hb(1),x0);
fmaxVec = pulseFcn2ndDeg(c,w,h,A(2,:),T,L,Hb(2),x0);


% SIMULATE THE SIGNAL
sigma = 0.000001;
r = 0.9;

% simulate driving noise:
e = sqrt(sigma)*randn(L, 1); 
% we need a loop since the AR(2) coefficients are now time varying
y = zeros(1,L);

for i=3:L
    % Create the arma:
    Af = fmax2AR2(fmaxVec(i),r);
    y(i) = -y(i-1)*Af(2) - y(i-2)*Af(3) + e(i);
end

% FILTER THROUGH THE AMPLITUDE FUNCTION TO GET THE COMPLETE SIGNAL:
y = y.*ampFcn;

end

