% Import files from Lab_4 and save mat file




for i=1:3
    [ fc , Dipole(i,:)]=importPulse(sprintf('Dipole_%d',i)); 
    [ ~ , highImp(i,:)]=importPulse(sprintf('High_imp_%d',i));
    [ ~ , lowImp(i,:)]=importPulse(sprintf('Low_imp_%d',i));
end

[ ~ , Dipole_wall]=importPulse('Dipole_wall');
[ ~ , highImp_wall]=importPulse('High_imp_wall');
[ ~ , lowImp_wall]=importPulse('Low_imp_wall');
[ ~ , RT]=importPulse('RT');

save('datas/Lab_4/Lab_4_data',"RT","lowImp_wall","highImp_wall","Dipole_wall","lowImp","highImp","Dipole","fc");