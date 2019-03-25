function [t,i,fullrate,breathtimes]=run_on_timer(t,i, fullrate, breathtimes)

%This is the function that I need to have run every .1 seconds or so. 
%I need to pass it everything that comes out of it to get carried or used
%in next instance of timer running. 
timepercond=12; %Should normally be 120 but just in case they don't quite 
%line up, you have some extra time.
if t.BytesAvailable>0
    [VO2, VCO2]=dataparse(fread(t,t.BytesAvailable));
    if [VO2, VCO2]~=[0,0]
        [metrate, goodbreath]=metabolic_data(VO2, VCO2);
        time_pass=toc;
        if time_pass<timepercond && goodbreath==1
            i=i+1;
            [fullrate, breathtimes]=store_rate(fullrate, metrate, breathtimes, time_pass, i);
        elseif time_pass>timepercond
            if goodbreath==1
                i=i+1;
                [fullrate, breathtimes]=store_rate(fullrate, metrate, breathtimes, time_pass, i);
            end
         %  conditiondone=1;
        end  
    end
end

