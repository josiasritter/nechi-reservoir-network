clear all

dischargeP = [0:0.1:1];

efficiencyP(1) = 0;
efficiencyP(2) = 0.69;
efficiencyP(3) = 0.83;
efficiencyP(4) = 0.87;
efficiencyP(5) = 0.88;
efficiencyP(6) = 0.89;
efficiencyP(7) = 0.9;
efficiencyP(8) = 0.9;
efficiencyP(9) = 0.9;
efficiencyP(10) = 0.89;
efficiencyP(11) = 0.88;

powerOutputP_per_head = dischargeP .* efficiencyP * 9800 / 1000000; % [MW/ m / m3/s]
% by multiplying this with the maxflow per turbine and the actual Net Head,it returns the power output per turbine [MW]

save Pelton_standardised