Training_2002-2011:
-------------------

Cascade_HUB_vs17_train:
Operational model of the cascade in training configuration used for optimisation of the 12 decision variables. Also contains the four objective functions for evaluation of operational strategies created by the optimisation algorithm (we used GODLIKE: https://github.com/rodyo/FEX-GODLIKE).

targetGen_vs7_train:
Behavioural component of the operational model. This represents the operational strategy depending on 12 decision variables to be optimised. This component is called every time step by Cascade_HUB_vs17_train and makes the release decisions (depending on the decision variable values).

-------------------------------------------------------------------------

Validation_2012-2015:
---------------------

Cascade_HUB_vs17_valid:
Operational model of the cascade in validation configuration used for testing the optimised operational strategies (each consisting of 12 decision variables) on an "unseen" data period.

targetGen_vs7_valid:
Behavioural component of the operational model in validation configuration. It converts the optimised operational strategy to release decisions. This component is called every time step by Cascade_HUB_vs17_valid and makes the release decisions depending on the given operational strategy.

validation_loop:
Loop over the total of operational strategies r to be tested by Cascade_HUB_vs17_valid on the "unseen" data period. Evaluation of each strategy with regard to all 4 objective functions.

-------------------------------------------------------------------------

General:
--------

FDC_Shift:
Called by Cascade_HUB_vs17_valid to calculate the shift of the flow-duration-curve (FDC) downstream of the cascade in comparison to the FDC of a natural reference hydrograph

Peltoncurve:
Assumed Pelton-efficiency-curve depending on discharge (normalized)

turbineData_standardized:
Assumed Francis-efficiency-curve depending on discharge and hydraulic head (normalized)

durationCurve_vs3:
Computes the duration curve from given hydrographs and energy time series

non-dominated-sorting:
Sorting algorithm of the multi-objective optimisation algorithm NSGA-GA II

real_objfnceval_2012_2015:
Evaluation of the present (real) operational strategy of the operating company over the testing period by means of 4 objective functions.


---------------------------------------------------------------------------------



Data Availability:
------------------

Time series data (2000-2015) freely available from Compañia de Expertos en Mercados (http://informacioninteligente10.xm.com.co/hidrologia/Paginas/HistoricoHidrologia.aspx)
Daily values of:
- Natural reservoir inflows from drainage areas (does not include inflows from upstream dams)
- Reservoir storages
- Reservoir spills
- Power generation


For technical data on reservoirs and turbines please contact the operating company Empresas Públicas de Medellín (EPM) or the corresponding author of the article Josias Ritter (ritterjosias@gmail.com)