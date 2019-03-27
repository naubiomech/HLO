function [SS, y_bar]=fake_metabolic_fit_JB(ParamsForCondition)

if length(ParamsForCondition)==2    
    x=ParamsForCondition(1);
    y=ParamsForCondition(2);
    cost=(x-12)^2+(y-67)^2;
elseif length(ParamsForCondition)==1
    x=ParamsForCondition(1);
    cost=(x-12)^2;
else
    error("Number of parameters defined not compatible with current script.");
end
SS=cost;
y_bar=linspace(0,10);

end