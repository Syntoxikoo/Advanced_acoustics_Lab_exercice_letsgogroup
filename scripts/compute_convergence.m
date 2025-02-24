lx = 0.7; ly = 1;
z = 10;
f = linspace(10,10000,1000);

N_modes = [100,100];
xM = 0; yM = 0; % receiver
rS = [0.,0.,0]; %source

G_full = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],1); % use either N_modes array and empty 

% plot(f,20*log10(abs(1j*2*pi*f.*G_full)/2e-5)) # Discuss about which is the best way to compute it


nx = linspace(0,100,100);
ny = linspace(0,100,100);
approx_sol = zeros(length(nx),length(f));
conv = zeros(length(nx),1);
% 

% % Troubleshoot the function 
% figure(1)
% plot(f,20*log10(abs(G_full)/2e-5)) ; grid on 
% title("full response")
% figure; 
% for ii= 1:2
%     disp(round(nx(ii)))
%     approx_sol(ii,:) =  Gf_duct([xM,yM],rS, z, [lx,ly], f, [round(nx(ii)),round(ny(ii))],[],1);
%     plot(f,20*log10(abs(approx_sol(ii,:))/2e-5)); hold on
% end
eps = 1e-9;
for ii= 1:length(nx)
    disp("iter nb : "+round(nx(ii)))
    approx_sol(ii,:) =  Gf_duct([xM,yM],rS, z, [lx,ly], f, [round(nx(ii)),round(ny(ii))],[],1);
    % conv(ii) = mean((abs(G_full - approx_sol(ii,:))).^2);
    conv(ii) = mean((abs(G_full - approx_sol(ii,:)) ./ ((abs(G_full) + eps)).^2));
    if isnan(conv(ii))
        error("the function outputed a non definite value")
    end
end

%% Plot convergence 1D for trial
plot_convergence