[data,TXT,RAW]=xlsread('InOutModel_rev5.xls','inputs');

    inputDate = datetime(2011,12,31) + caldays(1:1461);
    inputYear = data(:,1);
    inputMonth = data(:,2);
    year_index = data(:,3);
    month_index = data(:,4);
    dayofmonth = data(:,31);
    
% POWER GENERATION    
    % Observed energy generation [MW]
    obs_gen_G = data(:,17)/1000/24;
    obs_gen_T = data(:,18)/1000/24;
    obs_gen_P2 = data(:,21)/1000/24;
    obs_gen_P3 = data(:,22)/1000/24;
    
    obs_gen_total = obs_gen_G + obs_gen_T + obs_gen_P2 + obs_gen_P3;
    
% NATURAL HYDROGRAPH
    % Observed inflows [m3/s]
    in_desviaciones = data(:,5);
    in_tenche = data(:,6);
    in_guadalupe = data(:,7);
    in_concepcion = data(:,8);
    in_grande = data(:,9);
    in_porce2 = data(:,10);
    in_porce3 = data(:,11);
    HydrographTS_nat = in_guadalupe + in_grande + in_porce2 + in_porce3; 
    
% DISCHARGE PORCE 3
    % Observed discharge
    [data_out,TXT3,RAW3]=xlsread('InOutModel_rev5.xls','Output_template');
    out_porce3_total = data_out(:,2);     
    
    
    months_count = max(month_index);
    year_count = max(year_index);

    
% OBJECTIVE FUNCTIONS    
    % 1. Objective function: Firm energy (MAX)
    % Average power at monthly partition:
        avg_power_month_MW = zeros(months_count-120,1);    
        for i = 121:months_count 
           nnn = month_index == i;
           avg_power_month_MW(i-120) = mean(obs_gen_total(nnn));
        end

        % Firm energy
        [dCurve_P_monthly,pEmp] = durationCurve_vs2(avg_power_month_MW);  
        p95_e = find(pEmp > 0.95);
        firm_pwr_2012_2015 = dCurve_P_monthly(p95_e(1));
    
    
    % 2. Objective function: Average energy (MAX)
        avg_pwr_2012_2015 = mean(obs_gen_total(3653:5113));
    
    
    % 3. Objective function: Flood hazard (MIN)
    % concerning yearly maxima: create RT-Qmax relationship and compute area under graph in between defined thresholds of Q in accordance with GOTTA

        %Find yearly maxima
        max_HydrographTS = zeros(year_count-10,2);
        for n = 11 : year_count
            nnnn = year_index == n;
            z = max(out_porce3_total(nnnn));
            max_HydrographTS(n-10,1) = z;
        end

        % Calculate return period 
        % RT = n / m with n number of years of measurement and m relative ranking
        max_HydrographTS(:,1) = sort(max_HydrographTS(:,1),'descend');
        n = size(max_HydrographTS(:,1));
        n = n(1,1);
        for m=1:n
            max_HydrographTS(m,2)= n / m;
        end 

        % RT-Qmax relationship
        Qbankfull = 250; % HAS TO BE REVISED IN ACCORDANCE WITH GOTTA

        q = 1;
        while max_HydrographTS(q,1) > Qbankfull
            q = q+1;
        end

        max_HydrographTS_cut = [];    
        max_HydrographTS_cut(:,1) = max_HydrographTS(1:q,1);
        max_HydrographTS_cut(:,2) = max_HydrographTS(1:q,2);

        RTbankfull = interp1(max_HydrographTS_cut(:,1),max_HydrographTS_cut(:,2),Qbankfull);

        % Area under graph RT vs. Qmax
        FHareas = zeros(q-1,1);
        FHareas(q-1) = (Qbankfull+max_HydrographTS_cut(q-1,1))/2 * (max_HydrographTS_cut(q-1,2)-RTbankfull);

        for i = 1:q-2
            FHareas(i) = (max_HydrographTS_cut(i,1)+max_HydrographTS_cut(i+1,1))/2 * (max_HydrographTS_cut(i,2)-max_HydrographTS_cut(i+1,2));
        end
        area_RTvsQmax_2012_2015 = sum(FHareas) - Qbankfull*(max_HydrographTS_cut(1,2)-RTbankfull);
        

    % 4. Objective function: Flow alteration (MIN)
        % Monthly FDC (over entire modelling period)
                fdc_shift = zeros(8,12);
                out_porce3_total_2012_2015 = out_porce3_total(3653:5113);
                HydrographTS_nat_2012_2015 = HydrographTS_nat(3653:5113);
                [FDC_all,pEmp] = durationCurve_vs2([out_porce3_total_2012_2015,HydrographTS_nat_2012_2015]); % HydrographTS_nat required after checking with EPM

                for m = 1:12
                    nn = inputMonth(3653:5113) == m;
                    [dCurve,pEmp] = durationCurve_vs2([out_porce3_total_2012_2015(nn),HydrographTS_nat_2012_2015(nn)]);

                    [fdc_deficit_m,fdc_surplus_m] = FDC_Shift(dCurve(:,2),dCurve(:,1),pEmp,[0.05,0.10,0.75,0.95,0.99]);

                    fdc_shift(1:4,m) = fdc_deficit_m;
                    fdc_shift(5:8,m) = fdc_surplus_m;
                end

                % Alteration as simple sum of all alteration values
                alteration_monthly = sum(fdc_shift);
                alteration_2012_2015 = sum(alteration_monthly);
            