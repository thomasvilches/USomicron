
#######
### New features will be implemented ###
#####

#### new calibration
dc = [1;map(y->79+y,0:11);map(y->130+y,0:14);map(y->175+y,0:27);map(y-> 220+y,0:36);map(y->294+y,0:14);map(y->336+y,0:53);map(y->417+y,0:15)]#]]
rc = [1.0;map(y->1.0-(0.1/12)*y,1:12);map(y->0.9-(0.10/15)*y,1:15);map(y-> 0.8+(0.06/28)*y,1:28);map(y-> 0.86-(0.2114/37)*y,1:37);map(y-> 0.6486+(0.1584/15)*y,1:15);map(y-> 0.807-(0.173/53)*y,1:54);map(y-> 0.634+(0.039/16)*y,1:16)]#]]
run_param_scen_cal(true,0.106,"usa",10,50,1,1,1,1,1,92,999,194,7,999,1,2,14,rc,dc,478,true,999,1,1)
run_param_scen_cal(true,0.106,"usa",10,50,1,1,1,1,1,92,999,194,7,999,2,2,14,rc,dc,478,false,999,1,1)



#### new calibration
dc = [1;map(y->79+y,0:11);map(y->130+y,0:14);map(y->175+y,0:27);map(y-> 220+y,0:36);map(y->294+y,0:14);map(y->336+y,0:53);map(y->417+y,0:15)]#]]
rc = [1.0;map(y->1.0-(0.1/12)*y,1:12);map(y->0.9-(0.10/15)*y,1:15);map(y-> 0.8+(0.06/28)*y,1:28);map(y-> 0.86-(0.2114/37)*y,1:37);map(y-> 0.6486+(0.1584/15)*y,1:15);map(y-> 0.807-(0.173/53)*y,1:54);map(y-> 0.634+(0.039/16)*y,1:16)]#]]

for r = 0.5:0.1:1.0, tt = 1.0:0.02:1.1
    run_param_scen_cal(false,0.106,"usa",10,50,1,1,1,1,1,92,999,194,7,446,1,2,14,rc,dc,862,true,999,1,1,r,tt,0.5,1)
end


for r = 0.5:0.1:1.0, tt = 1.0:0.02:1.1
    run_param_scen_cal(false,0.106,"usa",10,50,1,1,1,1,1,92,999,194,7,446,1,2,14,rc,dc,862,true,999,1,1,r,tt,0.75,1)
end


for r = 0.5:0.1:1.0, tt = 1.0:0.02:1.1
    run_param_scen_cal(false,0.106,"usa",10,50,1,1,1,1,1,92,999,194,7,446,1,2,14,rc,dc,862,true,999,1,1,r,tt,1.0,1)
end


