function [SSa, y_bar]=metabolic_average(fullrates, breathtimes, AverageLapSpeed)

side_by_side=[fullrates', breathtimes'];
side_trim=side_by_side(all(side_by_side,2),:); %This line of code is from 
%https://www.mathworks.com/matlabcentral/answers/40528-delete-rows-which-contain-at-least-1-zero

t=side_trim(:,2);
y_meas=side_trim(:,1);

%%  Constant Metabolic Cost Estimation
% After some overground tests, it was determined that a first-order fit
% over a span of 4 minutes may not be ideal for estimating the steady-state
% metabolic rate. Instead, average the metabolic data over the last minute
% of the condition and use that as the estimate for the steady-state rate.

estimated_dot_E = mean(y_meas(end-5:end))/AverageLapSpeed;      %Average last minute of metabolics
y_bar = estimated_dot_E*ones(length(y_meas),1); %Make a line of constant slope

SSa=estimated_dot_E;


end