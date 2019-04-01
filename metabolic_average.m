function [SSa, y_bar]=metabolic_average(fullrates, breathtimes)

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

estimated_dot_E = mean(y_meas(end-5:end));      %Average last minute of metabolics
y_bar = estimated_dot_E*ones(length(y_meas),1); %Make a line of constant slope

%find the error between the best-fit predicted response and the
%measurement vector
mean_squared_error = ((y_bar-y_meas)'*(y_bar-y_meas))/n_samp;

SSa=estimated_dot_E;


end