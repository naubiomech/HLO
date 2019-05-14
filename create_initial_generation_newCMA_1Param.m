function [N, xmean, sigma, lambda, pc, ps, B, D, C, mu, weights, mueff, ...
    invsqrtC, eigeneval, chiN, cc, cs, c1, cmu, damps, x]=create_initial_generation_newCMA_1Param(...
    NumberofParams, InitGuess, sigma_start, ConditionsPerGen, pc_init, ps_init,...
    B_init, D_init, C_init, Peak_torque, Min_torque )

N=NumberofParams
xmean=InitGuess %Initial guess for the parameters
sigma=sigma_start
lambda=ConditionsPerGen
pc=pc_init
ps=ps_init
B=B_init
D=D_init
C=C_init

%Beginning of code part
mu=lambda/2
weights = log(mu+1/2)-log(1:mu)'; % mu*1 array for weighted recombination
mu = floor(mu);
weights = weights/sum(weights);     % Normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2); % Variance-effectiveness of sum w_i x_i

invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
eigeneval = 0;                      % Track update of B and D
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % Expectation of
%   ||N(0,I)|| == norm(randn(N,1))


% Strategy parameter setting: Adaptation
cc = (4+mueff/N) / (N+4 + 2*mueff/N);  % Time constant for cumulation for C
cs = (mueff+2) / (N+mueff+5);  % t-const for cumulation for sigma control
c1 = 2 / ((N+1.3)^2+mueff);    % Learning rate for rank-one update of C
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));  % and for rank-mu update
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % Damping for sigma
% usually close to 1

 % Generate and evaluate lambda offspring
 %With the new plan of 8 controllers per generation, do 7 evenly spaced,
 %and then the mean. 
 theta = 0; %initialize theta. 
 spacing = linspace(-0.25,0.75,lambda-1);
    for k=1:lambda-1
        %x(:, k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C)
        %x(:, k) = apply_constraints(x(:,k), Peak_torque ); % Apply hard constraints
        % to the random generated parameters. See apply_constraints.m for
        % seperate definition of the function.
        
        %Replacing the above 2 lines that normally get used with our new
        %function to get a pretty, circular distribution around that
        %initial point. 
        x(:,k) = xmean + sigma*spacing(k);
%         x(:,k) = xmean + sigma*[cosd(theta)]%; sind(theta)];
%         theta = theta + 360/(lambda-1) + 5;
        x(:, k) = apply_constraints(x(:,k), Peak_torque, Min_torque ); % Apply hard constraints
        
    end
    randnums = randperm(lambda-1);
    xrand = x(:,randnums); %Randomize the results
    Max = max(xrand);
    Min = min(xrand);
    x(:,1) = Max;
    x(:,2) = Min;
    x(:,3) = xmean;
    x(:,1:3) = x(:,randperm(3)); %Randomize the first 3 parameters again
    x(:,4:lambda)= xrand(~ismember(xrand,[Max,Min]));
     
end

%So we now have a circle of values to test. 