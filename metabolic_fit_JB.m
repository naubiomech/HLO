function [SSa, y_bar]=metabolic_fit_JB(fullrates, breathtimes)

side_by_side=[fullrates', breathtimes'];
side_trim=side_by_side(all(side_by_side,2),:); %This line of code is from 
%https://www.mathworks.com/matlabcentral/answers/40528-delete-rows-which-contain-at-least-1-zero

tau=42; %seconds. This can be edited if you wish

t=side_trim(:,2);
y_meas=side_trim(:,1);

%Somehow, I need to insert into here the calculations to get the steady
%state metabolic rate estimated using the two minutes of data. 

%ALL OF THE FOLLOWING CODE IS ALMOST EXACTLY COPY PASTE OF THE CODE FOUND ONLINE IN THE
%SUPPLEMENTARY MATERIAL SECTION OF THE ZHANG 2017 SCIENCE ARTICLE

%%  Constant Metabolic Cost Estimation
% This function is a least-squares system identification of a constant
% input to a first-order dynamic system .
% This script finds a constant actual metabolic rate 'estimated_dot_E' that result in
% the least squared error between the predicted system response y_bar and
% the series of measurements y_meas.

%INPUTS
%t: a vector of times (in seconds) associated with each measurement of the system output
%y_meas: the measurements of the system output at each time in t
%tau: the time constant of the system (in seconds)

%OUTPUTS
% estimated_dot_E: the estimated metabolic rate of the period.
% y_bar: the best fit system outputs predicted by the polynomial
% relationship

% The system is modeled as a first order discrete dynamical system such that
% y(i) = (1-dt/tau)*y(i-1) + (dt/tau)*u(c,p)
% dt = t(i) - t(i-1)

%The system is equivalent to the forward Euler integration of a first-order
%continuos system of the form
% y_dot = (1/tau)*(u-y)

%The function in this script identifies the vector x_star that results in
%the least squared-error between y_bar and y_meas.

%x_star is the optimal solution to a specially formulated matrix equation
% A*x = y_meas, where x = [y_bar(1) u]'

% x_star can be found using the pseudo-inverse of A
% x_star = pinv(A)*y_meas .

% y_bar can be found by using x_star
% A*x_star = y_bar

% Adapted from the function of Wyatt Felt and C. David Remy at
% https://cn.mathworks.com/matlabcentral/fileexchange/51328-instantaneuos-cost-mapping
% Copyright@Juanjuan Zhang, Steven H Collins 11/30/2016

%Generate the matrix A
n_samp = length(t);
A = zeros(n_samp,2);
A(1,:) = [1,0];
for i = 2:length(t)
    for j = 1:2
        dt = t(i)-t(i-1);
        if j == 1
            A(i ,j) = A(i-1,j)*(1-dt/tau);
        else
            A(i ,j) = A(i-1,j)*(1-dt/tau) + (dt/tau);
        end
    end
end

%solve for the optimal parameters
x_star = pinv(A)*y_meas;
%solve for the best-fit predicted response
y_bar = A*x_star;
%find the error between the best-fit predicted response and the
%measurement vector
mean_squared_error = ((y_bar-y_meas)'*(y_bar-y_meas))/n_samp;
%solve for the optimal parameters
estimated_dot_E = x_star(2);
%disp(num2str(x_star(1)))
%disp(num2str(x_star(2)))

SSa=estimated_dot_E;


end