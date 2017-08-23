[data,TXT,RAW]=xlsread('POPSIZE2500_train.xls','Sheet1');

[r, c] = size(data);
obj_feval = zeros(r, 4);

for i = 1:r
    optpars = data(i,:);
    [obj_feval(i,:)] = Cascade_HUB_vs17_valid(optpars);
    i
end
    