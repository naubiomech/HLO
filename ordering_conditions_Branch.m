function orderedconds = ordering_conditions_Branch(SSdata)
% First three setpoints are incremented from low to high, second set of
% three setpoints are decremented from high to low. Rank each set of
% setpoints to find the one with minimum metabolic cost. 

firstSet = SSdata(1:3,:);
secondSet = SSdata(4:6,:);
orderedFirstSet = sortrows(firstSet, 1);
orderedSecondSet = sortrows(secondSet, 1);
orderedconds = [orderedFirstSet; orderedSecondSet];
orderedconds = [orderedconds, [1:3,1:3]'];