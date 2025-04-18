function [f,G] = measGreen(l,dl,S,fileN)
% function to calculate the Green's function from an imported file
% 
% f = frequency vector
% G = Green's function
% 
% l = length of the tube
% dl = distance between microphones at the outlet
% S = cross-sectional area of the tube
% fileN = number of the measurement files to import

if nargin < 4
    error('not enough input variables');
end

% speed of sound
c0=343;

% import frequnecy responses
filenameAB=sprintf('/datas/Lab_2/MEASUREMENT/%d_Frequency Response H1(p_B_source,p_A_source) - Input.txt',fileN);
filenameAC=sprintf('/datas/Lab_2/MEASUREMENT/%d_Frequency Response H1(Pressure C,p_A_source) - Input.txt',fileN);
filenameBC=sprintf('/datas/Lab_2/MEASUREMENT/%d_Frequency Response H1(Pressure C,p_B_source) - Input.txt',fileN);

[f, Hab] = importHxy(filenameAB);
[~, Hac] = importHxy(filenameAC);
[~, Hbc] = importHxy(filenameBC);

% calculate Green's function
w = 2.*pi.*f;
k = w./c0;

%G = sin(k.*dl)./(k.*S).*(Hac.*cos(k.*l)-abs(Hab).^2.*Hbc.*cos(k.*(l+dl)))/(cos(k.*l).^2-2.*real(Hab).*cos(k.*(l+dl))+abs(Hab).^2.*cos(k.*(l+dl)).^2);
G = sin(k.*dl)./(k.*S).*(Hac.*cos(k.*l)-abs(Hab).^2.*Hbc.*cos(k.*(l+dl)))./(cos(k.*l).^2-2.*real(Hab).*cos(k.*(l+dl))+abs(Hab).^2 .*cos(k.*(l+dl)).^2);

end