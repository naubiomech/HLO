function [SS, y_bar]=fake_metabolic_fit_JB(ParamsForCondition);

x=ParamsForCondition(1);
%y=ParamsForCondition(2);
cost=(x-12)^2;%+(y-67)^2;
SS=cost;
y_bar=linspace(0,10);

end