
function [dCurve,pEmp] = durationCurve_vs3(inputSeries)
[r c] = size(inputSeries);
dCurve = sort(inputSeries,1,'descend');    
pEmp = cumsum(ones(r,1))/(r);
end