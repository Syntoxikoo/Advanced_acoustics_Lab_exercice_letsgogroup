function Pa = diffFieldSP(pa2,RT,fc)
% calculates the Sound power from measurements done in the reverberation
% chamber at DTU in 1/3 octave bands.
% Inputs are vectors
% 
% Params
%   pa2 = room averaged mean square sound pressure
%   RT = reverberation time in each 1/3 octave band
%   fc = center frequency
% Returns
%   Pa = sound power per 1/3 octave band in dB re 1 pW

if length(pa2)~=24
    error('input arrays should have the legth of 24')
elseif length(RT)~=24
    error('input arrays should have the legth of 24')
elseif length(fc)~=24
    error('input arrays should have the legth of 24')
end


rho = 1.2; 
c = 343; 
V = 245; 
S = 240; % Surface ares of reverberation chamber in m2
w0 = 1e-12; % reference of 1 pW

micCorrection = zeros(1,24);
micCorrection(20:end) = [-0.1,-0.2,-0.6,-1,-1.6]; % B&K 4192 freefield correction

% calculate lambda
lambda = c./fc;


% initialize and calulate sound power

Pa=13.8.*V.*pa2./(rho.*c^2.*RT).*(1+(S.*lambda)./(8.*V));

Pa=10*log10(Pa./w0);

Pa=Pa+micCorrection;

end