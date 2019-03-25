function [meta_rate, noerror]=metabolic_data(vo2, vco2)

%This is from Vince Chiu
meta_rate = (16.477*vo2+4.484*vco2)/60; %Brockway equation. 

if vo2>0 && vo2<3000 && vco2>0 && vco2<3000 && meta_rate <3000 %reject out-of-range breath values
    noerror=1;
else
    noerror=0;
end