
#for the national estimations - NYC parameter - PRESYMP
dc = [1;map(y->73+y,0:11);map(y->127+y,0:21);map(y->180+y,0:27);map(y-> 211+y,0:38);map(y->288+y,0:14);map(y->335+y,0:66);map(y->442+y,0:26);map(y->486+y,0:17)]
rc = [0.935;map(y->0.935-(0.07/12)*y,1:12);map(y->0.865-(0.12/22)*y,1:22);map(y-> 0.745+(0.075/28)*y,1:28);map(y-> 0.82-(0.10/39)*y,1:39);map(y-> 0.72+(0.22/15)*y,1:15);map(y-> 0.94-(0.149/67)*y,1:67);map(y-> 0.791-(0.31/27)*y,1:27);map(y-> 0.481-(0.14/18)*y,1:18)]
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,1,rc,dc,528,true,425,70,0)
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,2,rc,dc,528,false,425,70,0)



#for the national estimations - NYC parameter - PRESYMP - hosp parameter = 1
dc = [1;map(y->73+y,0:11);map(y->127+y,0:21);map(y->180+y,0:27);map(y-> 211+y,0:38);map(y->288+y,0:14);map(y->335+y,0:66);map(y->442+y,0:26);map(y->486+y,0:17)]
rc = [0.935;map(y->0.935-(0.07/12)*y,1:12);map(y->0.865-(0.12/22)*y,1:22);map(y-> 0.745+(0.075/28)*y,1:28);map(y-> 0.82-(0.10/39)*y,1:39);map(y-> 0.72+(0.22/15)*y,1:15);map(y-> 0.94-(0.149/67)*y,1:67);map(y-> 0.791-(0.31/27)*y,1:27);map(y-> 0.481-(0.14/18)*y,1:18)]
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,3,rc,dc,528,true,425,70,0,1.0)
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,4,rc,dc,528,false,425,70,0,1.0)



#for the national estimations - NYC parameter - PRESYMP - New matrix of waning immunity
dc = [1;map(y->73+y,0:11);map(y->127+y,0:21);map(y->180+y,0:27);map(y-> 211+y,0:38);map(y->288+y,0:14);map(y->335+y,0:66);map(y->446+y,0:30);map(y->477+y,0:5);map(y->490+y,0:20)]#491
rc = [0.935;map(y->0.935-(0.07/12)*y,1:12);map(y->0.865-(0.12/22)*y,1:22);map(y-> 0.745+(0.075/28)*y,1:28);map(y-> 0.82-(0.10/39)*y,1:39);map(y-> 0.72+(0.235/15)*y,1:15);map(y-> 0.955-(0.15/67)*y,1:67);map(y-> 0.805-(0.4/31)*y,1:31);map(y-> 0.405+(0.06/6)*y,1:6);map(y-> 0.46-(0.1166667/21)*y,1:21)]#12
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,5,rc,dc,528,true,425,70,0)
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,6,rc,dc,528,false,425,70,0)
