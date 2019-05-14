function [xmean, mu, weights, ps, pc, C, c1, cmu, cc, sigma, cs, damps,...
    chiN, eigeneval, invsqrtC, counteval, x, lambda, Params_full, N, B, D, mueff]=...
    create_next_generation_newCMA_1Param(orderedconds, xmean, mu, weights,ps, pc, C, c1,...
    cmu, cc, sigma, cs, damps, chiN, eigeneval, invsqrtC, counteval, x, lambda, N, B, D, mueff, Peak_torque, Min_torque)

 % Sort by fitness and compute weighted mean into xmean
    orderedconds = orderedconds(any(orderedconds,2),:);   %Remove rows of all zeros
    metrates = orderedconds(:,1);
    idx = orderedconds(:,2);
    best_params_from_previous_gen = orderedconds(1, 3:length(orderedconds(1,:)));
    best_params_from_previous_as_col = best_params_from_previous_gen'; %Column vector
    xold = xmean;
    xmean = x(:,idx(1:mu))*weights;   % Recombination, new mean value
    
  % Cumulation: Update evolution paths
    ps = (1-cs)*ps ...
        + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(N+1);
    pc = (1-cc)*pc ...
        + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
    
    % Adapt covariance matrix C
    artmp = (1/sigma) * (x(:,idx(1:mu))-repmat(xold,1,mu));
    C = (1-c1-cmu) * C ...                  % Regard old matrix
        + c1 * (pc*pc' ...                 % Plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % Minor correction if hsig==0
        + cmu * artmp * diag(weights) * artmp'; % Plus rank mu update
    
   % Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
    
    % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
    if counteval - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
        eigeneval = counteval;
        C = triu(C) + triu(C,1)'; % Enforce symmetry
        [B,D] = eig(C);           % Eigen decomposition, B==normalized eigenvectors
        D = sqrt(diag(D));        % D is a vector of standard deviations now
        invsqrtC = B * diag(D.^-1) * B';
    end
    
    %Create the next generation the old way (creating tons to see!)
      %for k=1:lambda
      for k=1:lambda*10    
        x(:, k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C)
        x(:, k) = apply_constraints(x(:,k), Peak_torque, Min_torque); % Apply hard constraints
        % to the random generated parameters. See apply_constraints.m for
        % seperate definition of the function.
      end
      
      %New method to create the next generation. 
      cvar=C;
      sigma=sigma;
      avg=xmean; 
      [eigenvec, eigenval ] = eig(cvar);

        % Get the index of the largest eigenvector
        [largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
        largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

        % Get the largest eigenvalue
        largest_eigenval = max(max(eigenval));

        % Get the smallest eigenvector and eigenvalue
        if(largest_eigenvec_ind_c == 1)
            %smallest_eigenval = max(eigenval(:,2))
            %smallest_eigenvec = eigenvec(:,2);
            smallest_eigenval = max(eigenval(:,1))
            smallest_eigenvec = eigenvec(:,1);
        else
            smallest_eigenval = max(eigenval(:,1))
            smallest_eigenvec = eigenvec(1,:);
        end

        % Calculate the angle between the x-axis and the largest eigenvector
        %angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
        angle = atan2(largest_eigenvec(1), largest_eigenvec(1));
        
        % This angle is between -pi and pi.
        % Let's shift it such that the angle is between 0 and 2pi
        if(angle < 0)
            angle = angle + 2*pi;
        end
        
        theta_grid = linspace(0,2*pi);
        phi = angle;
        X0=avg(1);
        %Y0=avg(2);
        Y0=length(orderedconds)/2;
        a = largest_eigenval*sigma;
        b = smallest_eigenval*sigma;

        % the ellipse in x and y coordinates 
        ellipse_x_r  = a*cos( theta_grid );
        ellipse_y_r  = b*sin( theta_grid );

        thetas_for_testpoints = linspace(0, 360, lambda-1);
        thetas_for_testpoints = thetas_for_testpoints(1:lambda-2);
%         thetas_for_testpoints = linspace(0,1,lambda-1);
%         thetas_for_testpoints = thetas_for_testpoints(2:lambda-1);
        %lambda-2 because also want to test the mean and last times best
        %one. 
        
%         xtestpoints = a*thetas_for_testpoints;
%         ytestpoints = b*thetas_for_testpoints;
        xtestpoints = a*cosd(thetas_for_testpoints); 
        ytestpoints = b*sind(thetas_for_testpoints); 
        
        
        %Define a rotation matrix
        R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

        %let's rotate the ellipse to some angle phi
        r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

        testpoints_rot = [xtestpoints; ytestpoints]'*R + X0;  
        Max = max(testpoints_rot);  % Maximum setpoint
        Min = min(testpoints_rot);  % Minimum setpoint
        
        % Draw the error ellipse
        figure;
        plot((r_ellipse(:,1) + X0)*Peak_torque,(r_ellipse(:,2) + Y0),'b-')
        hold on;

        data1=x(1,:);
        %data2=x(2,:);
       
        % Plot the original data
        plot(data1*Peak_torque, 1:length(data1), 'k.');
        hold on;

        %Plot the test points
        plot((testpoints_rot(:,1))*Peak_torque,(testpoints_rot(:,2)),'rx')
        
      clearvars x
      WhereToPlace = ones(1,lambda);    %Need to make a vector of where to place new conditions
      WhereToPlace(3) = 0;      % 3rd parameter to test should be best from previous gen
      WhereToPlace(lambda) = 0; % Last parameter to test should be mean from previous gen
      WhereToPlace = logical(WhereToPlace); %Make this array a logical array
      PossibleLoc = 1:lambda;   % Vector from 1 to lambda indicating position of conditions
      NewParamLocations = PossibleLoc(WhereToPlace);    % Where to put the new parameters to test
      OldParamLocations = PossibleLoc(~WhereToPlace);   % Where to put the previous gen parameters
      x = zeros(1,lambda);      % Presize the conditions array
      x(1,1) = Max(1);
      x(1,2) = Min(1);
      x(1, NewParamLocations(3:end))=testpoints_rot(~ismember(testpoints_rot(:,1),[Max,Min]),1)';
      
      %To test again the best one from last time
      x(1,OldParamLocations(1:end-1)) = best_params_from_previous_as_col;
      
      %x(:,lambda-1)=best_params_from_previous_as_col;
      %To test what is predicted as the best for this one, xmean
      %x(:, lambda)=xmean;
      x(1,OldParamLocations(end))= xmean;
      
      randnums = randperm(3);
      for i = 1:length(randnums)
          space_hold(i) = x(randnums(i)); %Randomize the first 3 or 5 data points
      end
      x(1:3) = space_hold;
      
      disp('Before constraints')
      x
     % pause
      
      for i=1:length(x(1,:))
      x(:,i) = apply_constraints(x(:,i), Peak_torque, Min_torque);
      end
      
      
      disp('after constraints')
      x
     % pause
      x=x';
      Params_full=x;
     
      %Plot the final test points (some same, some constrained) 
       plot(x(:,1)*Peak_torque,1:length(x(:,1)),'gx')
       
     x=x' %Need in that format to work 
end