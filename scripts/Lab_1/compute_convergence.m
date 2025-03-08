% Compute for freq
lx = 0.7; ly = 1;
z = 10;
f = linspace(10,1000,1000);

N_modes = [8,8];
xM = 0; yM = 0; % receiver
rS = [0.,0.,0]; %source

G_full = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],1,false); % use either N_modes array and empty 


nx = 0:N_modes(1);
ny = 0:N_modes(2);
approx_sol = zeros(length(nx),length(f));
convf = zeros(length(nx),1);

for ii= 1:length(nx)
        approx_sol(ii,:) =  Gf_duct([xM,yM],rS, z, [lx,ly], f, [round(nx(ii)),round(ny(ii))],[],1,false);
        diff = abs(G_full-approx_sol(ii,:));
        convf(ii) = mean((diff).^2);
        if isnan(convf(ii))|| isinf(convf(ii))
            warning("Non-definite value detected at iteration %d", ii);
            break;
        end
end

lx = 0.7; ly = 1;
z = linspace(0,10,1000);
f = 1000;

N_modes = [8,8];
xM = 0; yM = 0; % receiver
rS = [0.,0.,0]; %source

G_full = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],1,false); % use either N_modes array and empty 



nx = 0:N_modes(1);
ny = 0:N_modes(2);
approx_sol = zeros(length(nx),length(z));
convz = zeros(length(nx),1);

for ii= 1:length(nx)
        approx_sol(ii,:) =  Gf_duct([xM,yM],rS, z, [lx,ly], f, [round(nx(ii)),round(ny(ii))],[],1,false);
        diff = abs(G_full.'-approx_sol(ii,:));
        convz(ii) = mean((diff).^2);
        if isnan(convz(ii))|| isinf(convz(ii))
            warning("Non-definite value detected at iteration %d", ii);
            break;
        end
end





% %% Conv 2D
% lx = 0.7; ly = 1;
% z = 0;
% f = linspace(10,10000,1000);
% 
% N_modes = [60,60];
% xM = 0; yM = 0; % receiver
% rS = [0.,0.,0]; %source
% 
% G_full = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],1,false); % use either N_modes array and empty 
% 
% nx = linspace(0,60,60);
% ny = linspace(0,60,60);
% 
% conv = zeros(length(nx),length(ny),1);
% 
% for ii= 1:length(nx)
%     for jj = 1:length(ny)
%         disp("iter nb : "+round(nx(ii))+round(ny(jj)))
%         approx_sol =  Gf_duct([xM,yM],rS, z, [lx,ly], f, [round(nx(ii)),round(ny(jj))],[],1,false);
%         diff = abs(G_full-approx_sol);
%         conv(ii,jj) = mean((diff).^2);
%         if isnan(conv(ii))|| isinf(conv(ii))
%             warning("Non-definite value detected at iteration %d", ii);
%             break;
%         end
%     end
% end


%%
% modal_density(0.7,1.,3,true)