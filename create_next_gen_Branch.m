function [x,xmean] = create_next_gen_Branch(orderedconds)
% orderedconds is an array containing the ordered metabolic data, condition number, 
% setpoint, and ranking of the increment and decrement data. The first three setpoints 
% are from the increment, the second set of three from the decrement.
% Average the ratings from each set of set points, determine the new median
% set point, and define the next set of set points.

Setpoints = orderedconds(:,3);                %Ordered set points
Rankings = orderedconds(:,4);                 %The rankings of ordered setpoints
LowIdx = Setpoints == min(Setpoints);         %Logical of rows of low set points
HighIdx = Setpoints == max(Setpoints);        %Logical of rows of high set points
MidIdx = ~(LowIdx+HighIdx);                   %Logical of rows of middle set points
LowSetpoint = min(Setpoints);                 %Low set point value
HighSetpoint = max(Setpoints);                %High set point value
MidSetpoint = Setpoints(MidIdx);              %Middle set point value   
MidSetpoint = MidSetpoint(1);

LowRanks = Rankings(LowIdx);                  %Rankings of low set point
HighRanks = Rankings(HighIdx);                %Rankings of high set point
MidRanks = Rankings(MidIdx);                  %Rankings of middle set point

AvgLowRank = mean(LowRanks);                  %Average ranking of low set point
AvgHighRank = mean(HighRanks);                %Average ranking of high set point
AvgMidRank = mean(MidRanks);                  %Average ranking of middle set point

RankedSetpoints = sortrows([LowSetpoint AvgLowRank;...
    MidSetpoint AvgMidRank; HighSetpoint AvgHighRank],2);

if length(unique(RankedSetpoints(1,2))) == 3 %If 3 unique values pick the first (ranking = [1, 1.5])
    xmean = RankedSetpoints(1,1);
elseif RankedSetpoints(1,2) == RankedSetpoints(2,2) %If ranking of first two equal
    xmean = (RankedSetpoints(1,1)+RankedSetpoints(2,1))/2; %Median is the average of two best setpoints
elseif RankedSetpoints(1,2) == RankedSetpoints(2,2) && ...
        RankedSetpoints(1,2) == RankedSetpoints(3,2) %If ranking for all is equal
    % Check for consistent set points
    % Check for lowest overall metabolic rate
    
else % As a final catch hand it off to the user to decide which setpoint to use
    disp(['   Setpoint ',' Avg Rank ','  Rank 1 ','  Rank 2 ','  Met Rate 1 ','  Met Rate 2 ']);
    disp(RankedSetpoints);
    Row = input('Input the ROW of the setpoint you want to use as the median in the next walk: ');
    xmean = RankedSetpoints(Row,1)
end
