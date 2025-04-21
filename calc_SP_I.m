% Compute the sound power from the instensity measurement
% Assuming the cube was 1m^2 -> surf = 5 m^2


w0 = 1e-12;
probCorrect = zeros(1,length(Dipole_IS));

probCorrect(end-10:end) = [0 0.1 0.2 0.4 0.6 1.0 1.7 2.5 3.8 5.8 8.3];



SP_datas = zeros(size(datas ));
for ii = 1:size(datas,1)
    % Convert to SP
    SP_datas(ii,:) = datas(ii,:) * surf;
    %to dB
    SP_datas(ii,:) = 10 * log10(SP_datas(ii,:)./w0);
    %Correct
    SP_datas(ii,:) = SP_datas(ii,:) + probCorrect;
end

