lx = 0.7; ly = 1;
z = 10;
f = linspace(10,10000,1000);

N_modes = [100,100];
xM = 0; yM = 0; % receiver
rS = [0.,0.,0]; %source

G_full = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],1); % use either N_modes array and empty 



nx = linspace(0,100,100);
ny = linspace(0,100,100);
approx_sol = zeros(length(nx),length(f));
conv = zeros(length(nx),1);

% for ii= 1:length(nx)
%     disp("iter nb : "+round(nx(ii)))
%     approx_sol(ii,:) =  Gf_duct([xM,yM],rS, z, [lx,ly], f, [round(nx(ii)),round(ny(ii))],[],1);
%     diff = abs(G_full-approx_sol(ii,:));
%     conv(ii) = mean((diff).^2);
%     if isnan(conv(ii))|| isinf(conv(ii))
%         warning("Non-definite value detected at iteration %d", ii);
%         break;
%     end
% end
% save convergence_z10.mat conv

load convergence_z0.mat
%% Plot convergence 1D for trial
plot_convergence

%%
% plot(f,20*log10(abs(1j*2*pi*f.*G_full)/2e-5)) # Discuss about which is the best way to compute it
