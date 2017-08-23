%% Hub code connecting the elements of the hydropower dam cascade in the Nechí catchment (Antioquia, Colombia)
% Author: Josias Ritter (The Nature Conservancy), ritterjosias@gmail.com
% July 2016


% RESERVOIR/PP CODES:
%====================
% _M:                 : Miraflores
% _G:                 : Guatron (Troneras)
    % _GT:            : Troneras
    % _G3:            : Guadalupe 3
    % _G4:            : Guadalupe 4
% _T:                 : La Tasajera (Riogrande2)
% _P2:                : Porce 2
% _P3:                : Porce 3


% INPUTS (ALL LOADED BY 'InputReader_vs5.m' FROM 'InOutModel_rev5.xls' AND SAVED IN 'inputdata.mat'):
%====================================================================================================
% nat_inflow          : time series of natural (i.e. unimpaired) reservoir inflows [m3/s]
% monthly_mean        : vector of total natural inflow in monthly averaged values (e.g. monthy_mean(3) is the Multiannual average system inflow in March) [m3/s]
% out_niquia          : time series of discharge through Niquía (needed for supplyRequirementTS of Riogrande2) [m3/s]
% out_riogrande1      : time series of discharge through Niquía (needed for supplyRequirementTS of Riogrande2) [m3/s]
% eFlowReq            : time series of instream flow req. between PP inlet and outlet [m3/s]
% supplyReq           : time series of water demands NOT hydropower, eg. water leaving reservoir in other directions [m3/s]
% initial_storage     : storage at time step 1 [m3]
% maxVol              : technical maximum USEFUL Volume of reservoir (corresponding maxReservoirLevel) [m3] 
% maxReservoirLevel   : reservoir top of conservation (technical max.)  |  THESE PARAMETERS MUST BE RELATIVE
% turbDatum           : tailwater elevation                             |  TO THE SAME DATUM
% installedPower      : Generation capacity (total) [MW]
% headLoss            : Head losses in the load pipe. e.g. 0.1 (is 10%)
% numengines          : number of engines 
% maxflow             : max. discharge per turbine [m3/s]
% ET                  : Open water evaporation rate [mm/day]
% inputMonth          : time series of the corresponding month
% inputYear           : time series of the corresponding year
% month_index         : running count of months over modelling period
% year_index          : running count of years over modelling period
% obs_volume          : time series of observed useful volume of the reservoir. For calibration purposes [m3]
% firmE_req           : annual requirement of firm energy [kWh/year]
% obs_E_price         : time series of observed mean bolsa energy price after international transactions [COP/kWh]

% OTHER INPUTS:
%==============
% VSE                 : Elevation-Volume-Surface curve [m, m3, m2]
% targetGen           : Function computing the daily generation target and allocating it between the four PPs [in % of installed power]
% turbineData_standarized : For Francis and Kaplan turbines. Table returning normalised power output depending on discharge and head (and their maximums)
% Pelton_standardised  : For Pelton turbines. Table returning normalised power output per meter of head depending on discharge and maxflow

% MAIN OUTPUTS (evaluations of objective functions as outputs to GODLIKE):
%=========================================================================
% Evaluation values of objective functions (all to be minimised):
    % obj_feval(1) = maxGen - firm_power_month       : Firm energy
    % obj_feval(2) = maxGen - avg_power_MW           : Average energy
    % obj_feval(3) = area_RTvsQmax                   : Flood hazard
    % obj_feval(4) = alteration                      : Flow alteration

% OTHER OUTPUTS:
%===============
% volumeTS            : time series of useful reservoir storage
% elevationTS         : time series of water level in the reservoir
% surfaceTS           : time series of open water surface of reservoir
% ET_ActualTS         : time series of evaporation from reservoir [m3/s]
% supplyTS            : time series of water NOT USED for electricity generation, except environmental flow between dam and turbines [m3/s]
% turbinedFlowTS      : time series of water used for electricity generation [m3/s]
% turbinedFlowTS_pred : time series of estimated turbined water for calculation of WB of downstream reservoir before computing generation target [m3/s]
% volumeTS_pred       : time series of estimated reservoir volume resulting from WB using turbinedFlowTS_pred of upstream plants before computing generation target [m3]
% powerTS             : time series of average daily power [MW]
% spillTS             : time series of spills [m3/s]
% eFlowRTS            : time series of % of fulfilment of instream flow requirement between PP inlet and outlet
% supplyRTS           : time series of % of fullfilment of water demands NOT USED for electricity generation, except environmental flow between dam and turbines
% energyDispRTS       : time series of % of fullfilment of electricity target generation
% HydrographTS        : time series of the impaired hydrograph downstream of Porce 3
% fdc_shift           : change in flow duration curve downstream of Porce 3
% avg_power_month_MW
% num_alerts
% num_criticals
% RT-Qmax curve

% PARAMETERS FOR OPTIMISATION (values as inputs from GODLIKE):
%==========================================================
% drywet_thrshld [0;250]        : threshold of mean discharge over the last 7 days that determines if we are in dry or wet season
% parQ [0;23.6]                   : release curve of Miraflores (in Cascade_HUB)
% G_filling_crit [0;1]          : in dry season, flow from M to GT is maxflow_M when filling of Troneras is lower than this value (in Cascade_HUB)
% p1 [0;15]                     : determining weight decline in inflow forecast (ranking method)
% op1 [0;1]                     : weight of w_cascadefill_ind in computation of targetGen_TOTAL
% op2 [0;1]                     : relative weight of price_ind in computation of targetGen_TOTAL.
% op3 [0;1]                     : relative weight of forecast_ind in computation of targetGen_TOTAL. Weight of lanina_ind is computed depending on op1 & op2 & op3
% flood_targetfill_P3 [0;1]     : target filling of Porce 3 in rainy season
% flood_targetfill_P2 [0;1]     : target filling of Porce 2 in rainy season
% FirmEtarget [0;1]             : target firm energy; i.e. minimum targetGen_TOTAL to create higher targets if average energy in current month is in danger of beeing too low (for increasing firm energy) 
% alloc_weight [0;1]            : weights between indicators in while loop of tagetGen_TOTAL allocation to individual PPs
% minFlowP3_w [0;2]             : factor for the minimum outflow from Porce 3, expressed as the historical observed minimum discharge of the natural Hydrograph in the current month 

% ADDITIONAL COMMENTS:
%=====================
% All volumes as useful storages


function[obj_feval] = Cascade_HUB_vs17_valid(optpars)

load inputdata;

drywet_thrshld = optpars(1);
parQ = optpars(2);
G_filling_crit = optpars(3);

volumeTS_M(3653) = obs_volume_M(3653);
volumeTS_G(3653) = obs_volume_G(3653);
volumeTS_T(3653) = obs_volume_T(3653);
volumeTS_P2(3653) = obs_volume_P2(3653);
volumeTS_P3(3653) = obs_volume_P3(3653);

for t = 3653:r
%% Water Balance predictions BEFORE computation of generation targets and release decisions

% 1. Miraflores (independent of target generation, so release already computed)

    % Evaporation and Water balance before release decision
    surfaceTS_M(t) = interp1(volumen_M,surface_M,volumeTS_M(t)); % surface del embalse al final del paso anterior
    ET_ActualTS_M(t) = ET_M(t)/1000 * surfaceTS_M(t) / 86400;   % ET_M in mm/day!
    volumeTS_M(t+1) = max(0,volumeTS_M(t) + (nat_inflow_M(t) - ET_ActualTS_M(t)) * 86400);  
   
    % Túnel de Tenche
    Q_available_M(t) = max(0,(maxflow_M - in_desviaciones(t) - in_concepcion(t) + eFlowReq_M(t))); % First use inflows from deviations and Concepción to fill Túnel de Tenche
   
    % Release decision simplified by function of release through Túnel de Tenche depending on filling of Miraflores vs filling of Troneras 
    vol_ratioTS(t) = (volumeTS_M(t)/maxvol_M) / (volumeTS_G(t)/maxvol_G);
    
    % Two-area plane functions
    if vol_ratioTS(t) > 1
        turbinedFlowTS_M(t) = min((maxflow_M*(volumeTS_M(t)/maxvol_M) + (parQ-maxflow_M)*(volumeTS_G(t)/maxvol_G)),Q_available_M(t)); % parQ to be varied between [0:23.6] (m3/s) (Capacity of tunel de tenche)
    else turbinedFlowTS_M(t) = min((parQ*(volumeTS_M(t)/maxvol_M)),Q_available_M(t));
    end
   
    if (volumeTS_G(t)/maxvol_G) < 0.2
        turbinedFlowTS_M(t) = Q_available_M(t);
    end
    
    % Dry season boost for entire system
    if inputMonth(t) < 5 || inputMonth(t) == 12
        if drywet_ind(t) < drywet_thrshld 
           if (volumeTS_G(t)/maxvol_G) < G_filling_crit
           turbinedFlowTS_M(t) = Q_available_M(t);
           end
        end
    end
    
    if turbinedFlowTS_M(t) > volumeTS_M(t+1)/86400
        turbinedFlowTS_M(t) = volumeTS_M(t+1)/86400;
    end
    
    volumeTS_M(t+1) = volumeTS_M(t+1) - turbinedFlowTS_M(t) * 86400;
    
    % seeking for spills
    if volumeTS_M(t+1) > maxvol_M
        spillTS_M(t) = (volumeTS_M(t+1) - maxvol_M) / 86400;
        volumeTS_M(t+1) = maxvol_M;
    end
    
    elevationTS_M(t+1) = interp1(volumen_M,cota_M,volumeTS_M(t+1));
    if isnan(elevationTS_M(t+1)) == 1
        elevationTS_M(t+1) = max(cota_M);
    end

    % Actual discharge through Túnel de Tenche
    Q_TunelTenche(t) = min(maxflow_M,maxflow_M - Q_available_M(t) + turbinedFlowTS_M(t) + spillTS_M(t));
    
    if Q_TunelTenche(t) < 0
        Q_TunelTenche(t) = 0;
    end
    
    elevationTS_M(t+1) = interp1(volumen_M,cota_M,volumeTS_M(t+1));
    
    
% 2. Guatron (Consisting of dam and turbines Troneras, plus run-of-river Guadalupe 3 & 4)

    % Evaporation and Water balance in Troneras before release decision
    surfaceTS_G(t) = interp1(volumen_G,surface_G,volumeTS_G(t)); % surface del embalse al final del paso anterior
    ET_ActualTS_G(t) = ET_G(t)/1000 * surfaceTS_G(t) / 86400;   % ET_G in mm/day!
    volumeTS_G(t+1) = max(0,volumeTS_G(t) + (Q_TunelTenche(t) + nat_inflow_G(t) - ET_ActualTS_G(t)) * 86400); 
    
    % There is no actual eflow requirement. Either way, for eflow of El Salto cascade, the water of Quebrada Cañagordas and El Cañal are assumed to be
    % sufficient to fulfill the eFlow requirement
    
   
% 3. Riogrande2/La Tasajera (assuming flows towards Niquía and Riogrande1 as fixed supplyTS)
    
    % Evaporation and Water balance in Riogrande2 before release decision
    surfaceTS_T(t) = interp1(volumen_T,surface_T,volumeTS_T(t)); % surface del embalse al final del paso anterior
    ET_ActualTS_T(t) = ET_T(t)/1000 * surfaceTS_T(t) / 86400;   % ET_T in mm/day!
    volumeTS_T(t+1) = max(0,volumeTS_T(t) + (nat_inflow_T(t) - ET_ActualTS_T(t)) * 86400);
    
    % Priority 1. Water supply towards Niquía, as it is domestic water supply 
    if volumeTS_T(t+1) > (supplyReq_T_niquia(t) * 86400)
        supplyTS_T_niquia(t) = supplyReq_T_niquia(t);
    else
        supplyTS_T_niquia(t) = volumeTS_T(t+1)/86400;           
    end
    
    supplyRTS_T_niquia(t) = supplyTS_T_niquia(t)/supplyReq_T_niquia(t); 
    volumeTS_T(t+1) = volumeTS_T(t+1) - supplyTS_T_niquia(t) * 86400;
    
    % Priority 2. Ecoflow towards Riogrande 1 (in River Grande)
    if volumeTS_T(t+1) > (eFlowReq_T(t) * 86400) 
        eFlowTS_T(t) = eFlowReq_T(t);
    else
        eFlowTS_T(t) = volumeTS_T(t+1)/86400;         
    end
    
    eFlowRTS_T(t) = eFlowTS_T(t)/eFlowReq_T(t);
    volumeTS_T(t+1) = volumeTS_T(t+1) - eFlowTS_T(t) * 86400;
    
    
% 4. Porce 2

    % Evaporation and Water balance before release decision
    surfaceTS_P2(t) = interp1(volumen_P2,surface_P2,volumeTS_P2(t)); % surface del embalse al final del paso anterior
    ET_ActualTS_P2(t) = ET_P2(t)/1000 * surfaceTS_P2(t) / 86400;   % ET_P2 in mm/day!
    turbinedFlowTS_pred_T(t) = maxflow_T * volumeTS_T(t+1) / maxvol_T;   % estimated inflow from La Tasajera required for predicted available water volume for energy generation
    volumeTS_pred_P2(t+1) = max(0,volumeTS_P2(t) + (nat_inflow_P2(t) + supplyTS_T_niquia(t) + eFlowTS_T(t) + turbinedFlowTS_pred_T(t) - ET_ActualTS_P2(t)) * 86400);

    
% 5. Porce 3

    % Evaporation and Water balance before release decision
    surfaceTS_P3(t) = interp1(volumen_P3,surface_P3,volumeTS_P3(t)); % surface del embalse al final del paso anterior
    ET_ActualTS_P3(t) = ET_P3(t)/1000 * surfaceTS_P3(t) / 86400;   % ET_P3 in mm/day!
    turbinedFlowTS_pred_G(t) = maxflow_GT * volumeTS_G(t+1) / maxvol_G;   % estimated inflow from Guatron required for predicted available water volume for energy generation
    turbinedFlowTS_pred_P2(t) = maxflow_P2 * volumeTS_pred_P2(t+1) / maxvol_P2;   % estimated inflow from Porce 2 required for predicted available water volume for energy generation
    volumeTS_pred_P3(t+1) = max(0,volumeTS_P3(t) + (nat_inflow_P3(t) + turbinedFlowTS_pred_P2(t) + turbinedFlowTS_pred_G(t) - ET_ActualTS_P3(t) - eFlowReq_P3(t)) * 86400);
    nethead_pred_P3 = (interp1(volumen_P3,cota_P3,(volumeTS_P3(t)+volumeTS_pred_P3(t+1))/2) - turbDatum_P3) * (1-headLoss_P3);
    
    
%% function targetGen
    % global target generation (and its allocation to power plants) based on storages, demands, forecast etc. at every time step t
    % providing as outputs: targetGen_TOTAL, targetGen_TOT_adj, targetGen_G, targetGen_T, targetGen_P2, targetGen_P3
    
    [targetGen_TOTAL, targetGen_TOT_adj, targetGen_G, targetGen_T, targetGen_P2, targetGen_P3] = ...
        targetGen_vs7_valid(t,volumeTS_M,maxvol_M,volumeTS_G,maxvol_G,volumeTS_T,maxvol_T,volumeTS_pred_P2,maxvol_P2,volumeTS_pred_P3,maxvol_P3, ...
        nethead_pred_P3,installedPower_GT,installedPower_G3,installedPower_G4,installedPower_T,installedPower_P2,installedPower_P3,maxGen, ... 
        inputMonth,Qin_monthly_mean,monthly_Qmin_nat,enso,price_ratio_monthly_mean,drywet_ind, ...
        optpars);    
    
    % saves to time series
    targetGenTS_TOTAL(t) = targetGen_TOTAL;
    targetGenTS_TOT_adj(t) = targetGen_TOT_adj;
    targetGenTS_G(t) = targetGen_G;
    targetGenTS_T(t) = targetGen_T;
    targetGenTS_P2(t) = targetGen_P2;
    targetGenTS_P3(t) = targetGen_P3;

    
%% Release decisions

% Model 1. Miraflores independent of release decision and therefore already
% computed in WB (see above)


% Model 2. Guatron (Consisting of dam and turbines Troneras, plus run-of-river Guadalupe 3 & 4)

    % Power generation
    
    if targetGen_G > 0
        
        availableHead = interp1(volumen_G,cota_G,(volumeTS_G(t)+volumeTS_G(t+1))/2);
        if isnan(availableHead)
            availableHead = max(cota_G);  % the reservoir is spilling above the max tecnical level
        end

        head_loss = (availableHead - turbDatum_GT) * headLoss_GT + (-1) * turbDatum_G3 * headLoss_G3 + (-1) * turbDatum_G4 * headLoss_G4; % for simplicity assumed  as a percent of the total head. Hazen Williams could be a nice upgrade!
        nHead = (availableHead - turbDatum_GT) + (-1) * turbDatum_G3 + (-1) * turbDatum_G4 - head_loss;

        optimumRequiredWater = targetGen_G * (installedPower_GT + installedPower_G3 + installedPower_G4) * 1e6 / (9800 * nHead * 0.89) * 86400;  % 0.89 is target efficiency on a typical day. It's possible to achive because of the available turndown ratio ([1.5 * numEngines]:1) and demand distribution.
        optimumRequiredWater =  min(optimumRequiredWater,maxflow_G3 * numengines_G3 * 86400);   % G3 has lowest maxflow

        avWater = volumeTS_G(t+1);
        if avWater < 0
           avWater = 0;
        end

        if avWater >= optimumRequiredWater

            % Troneras
            for i = 1 : numengines_GT
                power_GT(i) = i * F_GT((optimumRequiredWater / 86400 / i) ,((availableHead-turbDatum_GT) * (1-headLoss_GT)));
                if power_GT(i) > (i/numengines_GT * installedPower_GT)
                    power_GT(i) = NaN;
                end
            end

            a =~ isnan(power_GT);         
            if max(a) > 0
                powerTS_GT(t) = max(power_GT);
            else
                powerTS_GT(t) = 0;   % is technically imposible to operate at given target
            end

            % Guadalupe 3
            for i = 1 : numengines_G3
                power_G3(i) = i * interp1(discharge_G3,powerOutputP_per_head_G3,optimumRequiredWater/86400/i) * ((-1)*turbDatum_G3) * (1-headLoss_G3);
            end

            a =~ isnan(power_G3);         
            if max(a) > 0
                powerTS_G3(t) = max(power_G3);
            else
                powerTS_G3(t) = 0;   % is technically imposible to operate at given target
            end        

            % Guadalupe 4
            for i = 1 : numengines_G4
                power_G4(i) = i * F_G4(optimumRequiredWater / 86400 / i ,((-1)*turbDatum_G4) * (1-headLoss_G4));
                if power_G4(i) > (i/numengines_G4 * installedPower_G4)
                    power_G4(i) = NaN;
                end
            end

            a =~ isnan(power_G4);         
            if max(a) > 0
                powerTS_G4(t) = max(power_G4);
            else
                powerTS_G4(t) = 0;   % is technically imposible to operate at given target
            end        

            % Guatron total
            powerTS_G(t) = powerTS_GT(t) + powerTS_G3(t) + powerTS_G4(t);

            if powerTS_G(t) > 0      % if it is technically not possible to operate at given target, hold on to the water!
                volumeTS_G(t+1) = volumeTS_G(t+1) - optimumRequiredWater; 
                turbinedFlowTS_G(t) = optimumRequiredWater/86400;
            else
                powerTS_G(t) = 0;
                powerTS_GT(t) = 0;
                powerTS_G3(t) = 0;
                powerTS_G4(t) = 0;
                turbinedFlowTS_G(t) = 0;
            end

        elseif avWater > 0  % Not enough available water to reach targetGen

            % Troneras
            for i = 1 : numengines_GT
                power_GT(i) = i * F_GT(avWater / 86400 / i ,((availableHead-turbDatum_GT) * (1-headLoss_GT)));
                if power_GT(i) > (i/numengines_GT * installedPower_GT)
                    power_GT(i) = NaN;
                end
            end

            a =~ isnan(power_GT);         
            if max(a) > 0
                powerTS_GT(t) = max(power_GT);
            else
                powerTS_GT(t) = 0;   % is technically imposible to operate at given target
            end

            % Guadalupe 3
            for i = 1 : numengines_G3
                power_G3(i) = i * interp1(discharge_G3,powerOutputP_per_head_G3,avWater/86400/i) * ((-1)*turbDatum_G3 * (1-headLoss_G3));
            end

            a =~ isnan(power_G3);         
            if max(a) > 0
                powerTS_G3(t) = max(power_G3);
            else
                powerTS_G3(t) = 0;   % is technically imposible to operate at given target
            end        

            % Guadalupe 4
            for i = 1 : numengines_G4
                power_G4(i) = i * F_G4(avWater / 86400 / i ,((-1)*turbDatum_G4) * (1-headLoss_G4));
                if power_G4(i) > (i/numengines_G4 * installedPower_G4)
                    power_G4(i) = NaN;
                end
            end

            a =~ isnan(power_G4);         
            if max(a) > 0
                powerTS_G4(t) = max(power_G4);
            else
                powerTS_G4(t) = 0;   % is technically imposible to operate at given target
            end        

            % Guatron total
            powerTS_G(t) = powerTS_GT(t) + powerTS_G3(t) + powerTS_G4(t);

            if powerTS_G(t) > 0      % if it is technically not possible to operate at given target, hold on to the water!
                volumeTS_G(t+1) = volumeTS_G(t+1) - avWater; 
                turbinedFlowTS_G(t) = avWater/86400;
            else
                powerTS_G(t) = 0;
                powerTS_GT(t) = 0;
                powerTS_G3(t) = 0;
                powerTS_G4(t) = 0;
                turbinedFlowTS_G(t) = 0;
            end

        else  
            powerTS_G(t) = 0;
        end

        energyDispRTS_G(t) = powerTS_G(t)/(targetGenTS_G(t) * (installedPower_GT + installedPower_G3 + installedPower_G4));
    
    else energyDispRTS_G(t) = 1;
    end
    
    
    % seeks for overflows:
    
    if volumeTS_G(t+1) > maxvol_G
        spillTS_G(t) = (volumeTS_G(t+1) - maxvol_G) / 86400;
        volumeTS_G(t+1) = maxvol_G;
    end
    
    elevationTS_G(t+1) = interp1(volumen_G,cota_G,volumeTS_G(t+1));
    if isnan(elevationTS_G(t+1)) == 1
        elevationTS_G(t+1) = max(cota_G);
    end

    
% Model 3. Riogrande2/La Tasajera

    % Power generation
    
    if targetGen_T > 0
    
        availableHead = interp1(volumen_T,cota_T,(volumeTS_T(t)+volumeTS_T(t+1))/2);
        if isnan(availableHead)
            availableHead = max(cota_T);  % the reservoir is spilling above the max tecnical level
        end

        head_loss = (availableHead - turbDatum_T) * headLoss_T; % for simplicity assumed as a percent of the total head. Hazen Williams could be a nice upgrade!
        nHead = availableHead - head_loss - turbDatum_T;

        optimumRequiredWater = targetGen_T * installedPower_T * 1e6 / (9800 * nHead * 0.85) * 86400;  % 0.85 is target efficiency in a typical day for Pelton. It's possible to achive because of the available turndown ratio ([1.5 * numEngines]:1) and demand distribution.
        optimumRequiredWater =  min(optimumRequiredWater,maxflow_T * numengines_T * 86400);

        avWater = volumeTS_T(t+1);
        if avWater < 0
           avWater = 0;
        end

        if avWater >= optimumRequiredWater

            for i = 1 : numengines_T
                power_T(i) = i * interp1(discharge_T,powerOutputP_per_head_T,optimumRequiredWater/86400/i) * nHead;
            end

            a =~ isnan(power_T);         
            if max(a) > 0
                powerTS_T(t) = max(power_T);
            else
                powerTS_T(t) = 0;   % is technically imposible to operate at given target
            end

            if powerTS_T(t) > 0             % if it is technically not possible to operate at given target, hold on to the water!
                volumeTS_T(t+1) = volumeTS_T(t+1) - optimumRequiredWater; 
                turbinedFlowTS_T(t) = optimumRequiredWater/86400;
            else
                powerTS_T(t) = 0;
                turbinedFlowTS_T(t) = 0;
            end

        elseif avWater > 0  % Not enough available water to reach targetGen

            for i = 1 : numengines_T
                power_T(i) = i * interp1(discharge_T,powerOutputP_per_head_T,avWater/86400/i) * nHead;
            end

            a =~ isnan(power_T);         
            if max(a) > 0
                powerTS_T(t) = max(power_T);
            else
                powerTS_T(t) = 0;   % is technically imposible to operate at given target
            end

            if powerTS_T(t) > 0             % if it is technically not possible to operate at given target, hold on to the water!
                volumeTS_T(t+1) = volumeTS_T(t+1) - avWater; 
                turbinedFlowTS_T(t) = avWater/86400;
            else
                powerTS_T(t) = 0;
                turbinedFlowTS_T(t) = 0;
            end

        else  
            powerTS_T(t) = 0;
        end

        energyDispRTS_T(t) = powerTS_T(t)/(targetGenTS_T(t) * installedPower_T);
    else energyDispRTS_T(t) = 1;
    end
    
    
    % seeks for overflows:
    if volumeTS_T(t+1) > maxvol_T
        spillTS_T(t) = (volumeTS_T(t+1) - maxvol_T) / 86400;
        volumeTS_T(t+1) = maxvol_T;
    end
    
    elevationTS_T(t+1) = interp1(volumen_T,cota_T,volumeTS_T(t+1));
    if isnan(elevationTS_T(t+1)) == 1
        elevationTS_T(t+1) = max(cota_T);
    end
    

% Model 4. Porce 2
    
    % true WB before release decision:
    
    volumeTS_P2(t+1) = max(0,volumeTS_P2(t) + (nat_inflow_P2(t) + supplyTS_T_niquia(t) + eFlowTS_T(t) + turbinedFlowTS_T(t) + spillTS_T(t) - ET_ActualTS_P2(t)) * 86400);
    
    % Power generation:
    
    if targetGen_P2 > 0
        
        availableHead = interp1(volumen_P2,cota_P2,(volumeTS_P2(t)+volumeTS_P2(t+1))/2);
        if isnan(availableHead)
            availableHead = max(cota_P2);  % the reservoir is spilling above the max tecnical level
        end

        head_loss = (availableHead - turbDatum_P2) * headLoss_P2; % for simplicity assumed  as a percent of the total head. Hazen Williams could be a nice upgrade!
        nHead = availableHead - head_loss - turbDatum_P2;
        if nHead > max(netHead_P2)
           nHead = max(netHead_P2);
        end

        optimumRequiredWater = targetGen_P2 * installedPower_P2 * 1e6 / (9800 * nHead * 0.89) * 86400;  % 0.89 is target efficiency in a typical day. It's possible to achive because of the available turndown ratio ([1.5 * numEngines]:1) and demand distribution.
        optimumRequiredWater =  min(optimumRequiredWater,maxflow_P2 * numengines_P2 * 86400);

        avWater = volumeTS_P2(t+1);
        if avWater < 0
           avWater = 0;
        end

        if avWater >= optimumRequiredWater

            for i = 1 : numengines_P2
                power_P2(i) = i * F_P2(optimumRequiredWater / 86400 / i ,nHead);
                if power_P2(i) > (i/numengines_P2 * installedPower_P2)
                    power_P2(i) = NaN;
                end
            end

            a =~ isnan(power_P2);         
            if max(a) > 0
                powerTS_P2(t) = max(power_P2);
            else
                powerTS_P2(t) = 0;   % is technically imposible to operate at given target
            end

            if powerTS_P2(t) > 0     % if it is technically not possible to operate at given target, hold on to the water!
                volumeTS_P2(t+1) = volumeTS_P2(t+1) - optimumRequiredWater; 
                turbinedFlowTS_P2(t) = optimumRequiredWater/86400;
            else
                powerTS_P2(t) = 0;
                turbinedFlowTS_P2(t) = 0;
            end

        elseif avWater > 0  % Not enough available water to reach targetGen

            for i = 1 : numengines_P2
                power_P2(i) = i * F_P2(avWater / 86400 / i ,nHead);
                if power_P2(i) > (i/numengines_P2 * installedPower_P2)
                    power_P2(i) = NaN;
                end
            end

            a =~ isnan(power_P2);         
            if max(a) > 0
                powerTS_P2(t) = max(power_P2);
            else
                powerTS_P2(t) = 0;   % is technically imposible to operate at given target
            end

            if powerTS_P2(t) > 0     % if it is technically not possible to operate at given target, hold on to the water!
                volumeTS_P2(t+1) = volumeTS_P2(t+1) - avWater; 
                turbinedFlowTS_P2(t) = avWater/86400;
            else
                powerTS_P2(t) = 0;
                turbinedFlowTS_P2(t) = 0;
            end

        else  
            powerTS_P2(t) = 0;
        end

        energyDispRTS_P2(t) = powerTS_P2(t)/(targetGenTS_P2(t) * installedPower_P2);
    else energyDispRTS_P2(t) = 1;
    end
    
    
    % seeks for overflows:

    if volumeTS_P2(t+1) > maxvol_P2
        spillTS_P2(t) = (volumeTS_P2(t+1) - maxvol_P2) / 86400;
        volumeTS_P2(t+1) = maxvol_P2;
    end
    
    elevationTS_P2(t+1) = interp1(volumen_P2,cota_P2,volumeTS_P2(t+1));
    if isnan(elevationTS_P2(t+1)) == 1
        elevationTS_P2(t+1) = max(cota_P2);
    end
    
    
% Model 5. Porce 3
    
    % true WB before release decision:
    
    volumeTS_P3(t+1) = max(0,volumeTS_P3(t) + (nat_inflow_P3(t) + turbinedFlowTS_P2(t) + turbinedFlowTS_G(t) + spillTS_P2(t) + spillTS_G(t) - ET_ActualTS_P3(t)) * 86400);
    
    if volumeTS_P3(t+1) > (eFlowReq_P3(t) * 86400)
        eFlowTS_P3(t) = eFlowReq_P3(t);
    else
        eFlowTS_P3(t) = volumeTS_P3(t+1)/86400;         
    end
    
    eFlowRTS_P3(t) = eFlowTS_P3(t)/eFlowReq_P3(t);
    volumeTS_P3(t+1) = volumeTS_P3(t+1) - eFlowTS_P3(t) * 86400;
    
    % Power generation:
    
    if targetGen_P3 > 0
        
        availableHead = interp1(volumen_P3,cota_P3,(volumeTS_P3(t)+volumeTS_P3(t+1))/2);
        if isnan(availableHead)
            availableHead = max(cota_P3);  % the reservoir is spilling above the max tecnical level 
        end

        head_loss = (availableHead - turbDatum_P3) * headLoss_P3; % for simplicity assumed  as a percent of the total head. Hazen Williams could be a nice upgrade!
        nHead = availableHead - head_loss - turbDatum_P3;
        if nHead > max(netHead_P3)
           nHead = max(netHead_P3);
        end

        optimumRequiredWater = targetGen_P3 * installedPower_P3 * 1e6 / (9800 * nHead * 0.89) * 86400;  % 0.89 is target efficiency in a typical day. It's possible to achive because of the available turndown ratio ([1.5 * numEngines]:1) and demand distribution.
        optimumRequiredWater =  min(optimumRequiredWater,maxflow_P3 * numengines_P3 * 86400);

        avWater = volumeTS_P3(t+1);
        if avWater < 0
           avWater = 0;
        end


        if avWater >= optimumRequiredWater

            for i = 1 : numengines_P3
                power_P3(i) = i * F_P3(optimumRequiredWater / 86400 / i ,nHead);
                if power_P3(i) > (i/numengines_P3 * installedPower_P3)
                    power_P3(i) = NaN;
                end
            end

            a =~ isnan(power_P3);         
            if max(a) > 0
                powerTS_P3(t) = max(power_P3);
            else
                powerTS_P3(t) = 0;   % is technically imposible to operate at given target
            end

            if powerTS_P3(t) > 0     % if it is technically not possible to operate at given target, hold on to the water!
                volumeTS_P3(t+1) = volumeTS_P3(t+1) - optimumRequiredWater; 
                turbinedFlowTS_P3(t) = optimumRequiredWater/86400;
            else
                powerTS_P3(t) = 0;
                turbinedFlowTS_P3(t) = 0;
            end

        elseif avWater > 0  % Not enough available water to reach targetGen

            for i = 1 : numengines_P3
                power_P3(i) = i * F_P3(avWater / 86400 / i ,nHead);
                if power_P3(i) > (i/numengines_P3 * installedPower_P3)
                    power_P3(i) = NaN;
                end
            end

            a =~ isnan(power_P3);         
            if max(a) > 0
                powerTS_P3(t) = max(power_P3);
            else
                powerTS_P3(t) = 0;   % is technically imposible to operate at given target
            end

            if powerTS_P3(t) > 0     % if it is technically not possible to operate at given target, hold on to the water!
                volumeTS_P3(t+1) = volumeTS_P3(t+1) - avWater; 
                turbinedFlowTS_P3(t) = avWater/86400;
            else
                powerTS_P3(t) = 0;
                turbinedFlowTS_P3(t) = 0;
            end

        else  
            powerTS_P3(t) = 0;
        end

        energyDispRTS_P3(t) = powerTS_P3(t)/(targetGenTS_P3(t) * installedPower_P3);
    else energyDispRTS_P3(t) = 1;
    end
    
    
    % seeks for overflows:
 
    if volumeTS_P3(t+1) > maxvol_P3
        spillTS_P3(t) = (volumeTS_P3(t+1) - maxvol_P3) / 86400;
        volumeTS_P3(t+1) = maxvol_P3;
    end
    
    elevationTS_P3(t+1) = interp1(volumen_P3,cota_P3,volumeTS_P3(t+1));
    if isnan(elevationTS_P3(t+1)) == 1
        elevationTS_P3(t+1) = max(cota_P3);
    end
    
    
% Overall cascade performance    
    
    powerTS_TOTAL(t) = powerTS_G(t) + powerTS_T(t) + powerTS_P2(t) + powerTS_P3(t);
    energyDispRTS_TOTAL(t) = powerTS_TOTAL(t) / (targetGenTS_TOTAL(t)*(installedPower_G3+installedPower_G4+installedPower_GT+installedPower_P2+installedPower_P3+installedPower_T));
    energyDispRTS_TOT_adj(t) = powerTS_TOTAL(t) / (targetGenTS_TOT_adj(t)*(installedPower_G3+installedPower_G4+installedPower_GT+installedPower_P2+installedPower_P3+installedPower_T));
    
    HydrographTS(t) = turbinedFlowTS_P3(t) + spillTS_P3(t) + eFlowTS_P3(t);
    
    
end 

%% Evaluation

obj_feval = zeros(1,4);
months_count = max(month_index);

% 1. Objective function: Firm energy (MAX)
    % Average power at monthly partition:
    avg_power_month_MW = zeros(months_count,1);    
    for i = 121:months_count 
       nnn = month_index == i;
       avg_power_month_MW(i) = mean(powerTS_TOTAL(nnn));
    end

    % Firm energy
    [dCurve_P_monthly,pEmp] = durationCurve_vs3(avg_power_month_MW(121:months_count));  
    firm_power_month = prctile(dCurve_P_monthly,5);
    
obj_feval(1,1) = maxGen - firm_power_month; % Since objective function has to be minimised


% 2. Objective function: Average energy (MAX)
    % Average power at yearly partition:
    avg_power_year_MW = zeros(year_count,1);    
    for i = 1:year_count
       nnn = year_index == i;
       avg_power_year_MW(i) = mean(powerTS_TOTAL(nnn));
    end
    
    [dCurve_P_yearly,pEmp] = durationCurve_vs3(avg_power_year_MW);
    total_energy_gwh = mean(powerTS_TOTAL)*24*(365*max(year_index)+3)/1000;
    avg_power_MW = mean(powerTS_TOTAL(3653:5113));
    
    obj_feval(1,2) = maxGen - avg_power_MW; % Since objective function has to be minimised
    
    
% 3. Objective function: Flood hazard (MIN)
    % concerning yearly maxima: create RT-Qmax relationship and compute area under graph in between defined thresholds of Q in accordance with GOTTA
    
    %Find yearly maxima
    max_HydrographTS = zeros(year_count-10,2);
    for n = 11 : year_count
        nnnn = year_index == n;
        z = max(HydrographTS(nnnn));
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
    % Area under graph RT vs. Qmax
    Qbankfull = 250; % HAS TO BE REVISED IN ACCORDANCE WITH GOTTA
    
    FHareas = zeros(year_count-10,1);
    for i = 1:year_count-11
        FHareas(i) = (max_HydrographTS(i,1)+max_HydrographTS(i+1,1))/2 * (max_HydrographTS(i,2)-max_HydrographTS(i+1,2));
    end
    FHareas(year_count-10) = max_HydrographTS(year_count-10,1)/2 * max_HydrographTS(year_count-10,2);
    area_RTvsQmax = sum(FHareas) - Qbankfull*max_HydrographTS(1,2);

    if Qbankfull > min(max_HydrographTS)
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
    area_RTvsQmax = sum(FHareas) - Qbankfull*(max_HydrographTS_cut(1,2)-RTbankfull);
    end
    
obj_feval(1,3) = area_RTvsQmax; % area under graph RT vs Qmax(yearly) above bankfull discharge
    

% 4. Objective function: Flow alteration (MIN)
    % FDC over entire modelling period for altered (simulated) and for natural hydrograph
    [FDC_all,pEmp] = durationCurve_vs3([HydrographTS,HydrographTS_nat]); % HydrographTS_nat required after checking with EPM
        
    % Monthly FDC (over 2012-2015)
            fdc_shift = zeros(8,12);
            
                for m = 1:12
                    nn = inputMonth == m & month_index > 120;
                    [dCurve,pEmp] = durationCurve_vs2([HydrographTS(nn),HydrographTS_nat(nn)]);

                    [fdc_deficit_m,fdc_surplus_m] = FDC_Shift(dCurve(:,2),dCurve(:,1),pEmp,[0.05,0.10,0.75,0.95,0.99]);

                    fdc_shift(1:4,m) = fdc_deficit_m;
                    fdc_shift(5:8,m) = fdc_surplus_m;
                end
            
            % Alerts or critical events of streamflow regime alteration:            
            alerts_deficit_surplus = fdc_shift > 0.2;
            critical_deficit_surplus = fdc_shift > 0.5;

            eco_deficits = fdc_shift(1:4,:);
            eco_surplus = fdc_shift(5:8,:);

            num_alerts = sum(sum(alerts_deficit_surplus));
            num_criticals = sum(sum(critical_deficit_surplus));
            
            % Alteration as simple sum of all alteration values
            alteration_monthly = sum(fdc_shift);
            alteration = sum(alteration_monthly);
            
obj_feval(1,4) = alteration; % sum of alteration values in ELOHA dashboard

end