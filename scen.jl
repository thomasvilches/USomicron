
#for the national estimations - NYC parameter - PRESYMP - New matrix of waning immunity - fairly okay
dc = [1;map(y->73+y,0:11);map(y->125+y,0:21);map(y->180+y,0:27);map(y-> 211+y,0:38);map(y->288+y,0:14);map(y->337+y,0:66);map(y->446+y,0:30);map(y->477+y,0:5);map(y->491+y,0:21);map(y->558+y,0:9)]#491
rc = [0.935;map(y->0.935-(0.07/12)*y,1:12);map(y->0.865-(0.12/22)*y,1:22);map(y-> 0.745+(0.075/28)*y,1:28);map(y-> 0.82-(0.10/39)*y,1:39);map(y-> 0.72+(0.235/15)*y,1:15);map(y-> 0.955-(0.135/67)*y,1:67);map(y-> 0.82-(0.41/31)*y,1:31);map(y-> 0.41+(0.05/6)*y,1:6);map(y-> 0.46-(0.16/22)*y,1:22);map(y-> 0.30+(0.04/10)*y,1:10)]#12
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,0,rc,dc,957,true,425,70,0)
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,1,rc,dc,957,true,425,70,0)
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,2,rc,dc,957,true,425,70,0)
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,3,rc,dc,957,true,425,70,0)

