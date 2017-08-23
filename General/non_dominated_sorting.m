[obj_fcn_val]=xlsread('POPSIZE2500_valid2.xls','4D_Par_of_validated');
[dec_var]=xlsread('POPSIZE2500_valid2.xls','4D_Par_dv');

firmE_4D = 1883 - obj_fcn_val(:,1);
avgE_4D = 1883 - obj_fcn_val(:,2);
floodhazard_4D = obj_fcn_val(:,3);
alteration_4D = obj_fcn_val(:,4);


D = size(obj_fcn_val);
n = D(1);
m = D(2);
E = size(dec_var);
l = E(2);


% computing the rank of domination for each member
rank = zeros(n,1);

for i = 1:n
    for k = 1:n
        DOM1 = zeros(1,m);
        DOM2 = zeros(1,m);
        DOM = 0;
        for j = 1:m
            if obj_fcn_val(i,j) >= obj_fcn_val(k,j)
                DOM1(j) = 1;
            end
            if obj_fcn_val(i,j) > obj_fcn_val(k,j)
                DOM2(j) = 1;
            end
        end
        if sum(DOM1) == m
            DOM = 1;
        end
        if sum(DOM2) == 0
            DOM = 0;
        end
        if DOM == 1
            rank(i) = rank(i) + 1;
        end
    end
end


% selecting members with rank of domination = 0
B = rank == 0;
N = sum(B);

Index = find(rank==0);

newPareto = zeros(N,4);
new_dec_var = zeros(N,l);


% Filling new matrizes with non-dominated members
for h = 1:N
    p = Index(h);
    newPareto(h,:) = obj_fcn_val(p,:);
    new_dec_var(h,:) = dec_var(p,:);
end

firmE = 1883 - newPareto(:,1);
avgE = 1883 - newPareto(:,2);
floodhazard = newPareto(:,3);
alteration = newPareto(:,4);

