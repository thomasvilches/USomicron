
dc = [1;map(y->73+y,0:11);map(y->125+y,0:21);map(y->180+y,0:27);map(y-> 211+y,0:38);map(y->288+y,0:14);map(y->337+y,0:66);map(y->446+y,0:30);map(y->477+y,0:5);map(y->491+y,0:21);map(y->560+y,0:11);map(y->591+y,0:11);map(y->615+y,0:14);map(y->688+y,0:24);map(y->776+y,0:14);map(y->840+y,0:14)]#491
rc = [0.935;map(y->0.935-(0.07/12)*y,1:12);map(y->0.865-(0.12/22)*y,1:22);map(y-> 0.745+(0.075/28)*y,1:28);map(y-> 0.82-(0.10/39)*y,1:39);map(y-> 0.72+(0.235/15)*y,1:15);map(y-> 0.955-(0.135/67)*y,1:67);map(y-> 0.82-(0.41/31)*y,1:31);map(y-> 0.41+(0.05/6)*y,1:6);map(y-> 0.46-(0.16/22)*y,1:22);map(y-> 0.30+(0.05/12)*y,1:12);map(y-> 0.35+(0.025/12)*y,1:12);map(y-> 0.375-(0.04/15)*y,1:15);map(y-> 0.335-(0.037/25)*y,1:25);map(y-> 0.298+(0.0052/15)*y,1:15);map(y-> 0.3032-(0.028/15)*y,1:15)]#12

# #Calibration
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,0,1,rc,dc,[911;-1;-3;-4;-5;-7;-9],[0.0;0.0;0.0;0.0;0.0],true,425,70,0)

## Worst-case scenario
# vaccinate up to the end in the current pace March 31st
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,0,1,rc,dc,[941;-1;-3;-5;-7;-9;-11],[0.0;0.0;0.0;0.0;0.0],true,425,70,0,[[0.76; 0.82; 0.90], [0.677;0.82;0.90]])

#vaccinate Octuber, Nov, Dec
#50% flu coverage
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,4,2,rc,dc,[760;760;852;941;-5;-7;-9],[0.6;0.5;0.38;0.54;0.75]./2,true,425,70,0,[[0.76; 0.82; 0.90], [0.677;0.82;0.90]])

#vaccinate Octuber, Nov, Dec
# Flu coverage
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,4,3,rc,dc,[760;760;852;941;-5;-7;-9],[0.6;0.5;0.38;0.54;0.75],true,425,70,0,[[0.76; 0.82; 0.90], [0.677;0.82;0.90]])

#vaccinate Octuber, Nov, Dec
#Stop vaccinating in octuber
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,0,4,rc,dc,[760;941;-3;-5;-7;-9;-11],[0.0;0.0;0.0;0.0;0.0],true,425,70,0,[[0.76; 0.82; 0.90], [0.677;0.82;0.90]])




## Worst-case scenario
# vaccinate up to the end in the current pace March 31st
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,0,1,rc,dc,[941;-1;-3;-5;-7;-9;-11],[0.0;0.0;0.0;0.0;0.0],true,425,70,0,[[0.76; 0.82; 0.90], [0.677;0.82;0.90]], 0.5)

#vaccinate Octuber, Nov, Dec
#50% flu coverage
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,4,2,rc,dc,[760;760;852;941;-5;-7;-9],[0.6;0.5;0.38;0.54;0.75]./2,true,425,70,0,[[0.76; 0.82; 0.90], [0.677;0.82;0.90]], 0.5)

#vaccinate Octuber, Nov, Dec
# Flu coverage
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,4,3,rc,dc,[760;760;852;941;-5;-7;-9],[0.6;0.5;0.38;0.54;0.75],true,425,70,0,[[0.76; 0.82; 0.90], [0.677;0.82;0.90]], 0.5)

#vaccinate Octuber, Nov, Dec
#Stop vaccinating in octuber
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,0,4,rc,dc,[760;941;-3;-5;-7;-9;-11],[0.0;0.0;0.0;0.0;0.0],true,425,70,0,[[0.76; 0.82; 0.90], [0.677;0.82;0.90]], 0.5)
