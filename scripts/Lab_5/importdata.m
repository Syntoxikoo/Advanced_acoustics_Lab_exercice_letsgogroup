% Import files from Lab_5 and save mat file

[ f_IS , Dipole_IS]=importPI('Intensity Spectrum_DIPOLE');
[ ~ , Dipole_fan_IS]=importPI('Intensity Spectrum_DIPOLE+FAN');
[ ~ , Fan_IS]=importPI('Intensity Spectrum_FAN');
[ ~ , Omni_IS]=importPI('Intensity Spectrum_OMNI');
[ ~ , Omni_1_IS]=importPI('Intensity Spectrum_OMNI_1');

[ f_PI , Dipole_PI]=importPI('PI Index Spectrum(Signal 1,Signal 2)_DIPOLE');
[ ~ , Dipole_fan_PI]=importPI('PI Index Spectrum(Signal 1,Signal 2)_DIPOLE+FAN');
[ ~ , Fan_PI]=importPI('PI Index Spectrum(Signal 1,Signal 2)_FAN');
[ ~ , Omni_PI]=importPI('PI Index Spectrum(Signal 1,Signal 2)_OMNI');
[ ~ , Omni_1_PI]=importPI('PI Index Spectrum(Signal 1,Signal 2)_OMNI_1');

save('datas/Lab_5/Lab_5_data',"Dipole_IS","Dipole_fan_IS","Fan_IS","Omni_IS","Omni_1_IS","f_IS","Dipole_PI","Dipole_fan_PI","Fan_PI","Omni_PI","Omni_1_PI","f_PI");