
#for the national estimations - NYC parameter
dc = [1;map(y->48+y,0:6);map(y->77+y,0:11);map(y->127+y,0:21);map(y->177+y,0:27);map(y-> 222+y,0:38);map(y->288+y,0:14);map(y->330+y,0:62);map(y->442+y,0:26);map(y->489+y,0:17)]
rc = [1.0;map(y->1.0+(0.06/7)*y,1:7);map(y->1.06-(0.1/12)*y,1:12);map(y->0.96-(0.16/22)*y,1:22);map(y-> 0.8+(0.075/28)*y,1:28);map(y-> 0.875-(0.12/39)*y,1:39);map(y-> 0.755+(0.30/15)*y,1:15);map(y-> 1.055-(0.21/63)*y,1:63);map(y-> 0.845-(0.33/27)*y,1:27);map(y-> 0.515-(0.18/18)*y,1:18)]
run_param_scen_cal(true,0.096,"usa",50,1,1,1,1,1,1,rc,dc,518,true,425,70,1)
run_param_scen_cal(true,0.096,"usa",50,1,1,1,1,1,2,rc,dc,518,false,425,70,1) 

#getting good

#for the national estimations - NYC parameter
dc = [1;map(y->77+y,0:11);map(y->127+y,0:21);map(y->173+y,0:27);map(y-> 222+y,0:38);map(y->288+y,0:14);map(y->331+y,0:62);map(y->442+y,0:26)]
rc = [0.935;map(y->0.935-(0.08/12)*y,1:12);map(y->0.855-(0.16/22)*y,1:22);map(y-> 0.695+(0.075/28)*y,1:28);map(y-> 0.77-(0.10/39)*y,1:39);map(y-> 0.67+(0.25/15)*y,1:15);map(y-> 0.92-(0.17/63)*y,1:63);map(y-> 0.75-(0.31/27)*y,1:27)]
run_param_scen_cal(true,0.121,"usa",15,1,1,1,1,1,5,rc,dc,440,true,425,70,1)
run_param_scen_cal(true,0.121,"usa",15,1,1,1,1,1,6,rc,dc,440,false,425,70,1)


#for the national estimations - NYC parameter - PRESYMP
dc = [1;map(y->73+y,0:11);map(y->127+y,0:21);map(y->180+y,0:27);map(y-> 217+y,0:38);map(y->288+y,0:14);map(y->331+y,0:62);map(y->442+y,0:26)]
rc = [0.935;map(y->0.935-(0.07/12)*y,1:12);map(y->0.865-(0.12/22)*y,1:22);map(y-> 0.745+(0.075/28)*y,1:28);map(y-> 0.82-(0.10/39)*y,1:39);map(y-> 0.72+(0.28/15)*y,1:15);map(y-> 1.0-(0.13/63)*y,1:63);map(y-> 0.75-(0.31/27)*y,1:27)]
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,7,rc,dc,440,true,425,70,1)
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,8,rc,dc,440,false,425,70,1)


#for the national estimations - NYC parameter - PRESYMP
dc = [1;map(y->73+y,0:11);map(y->127+y,0:21);map(y->180+y,0:27);map(y-> 211+y,0:38);map(y->288+y,0:14);map(y->335+y,0:66);map(y->442+y,0:26);map(y->486+y,0:17)]
rc = [0.935;map(y->0.935-(0.07/12)*y,1:12);map(y->0.865-(0.12/22)*y,1:22);map(y-> 0.745+(0.075/28)*y,1:28);map(y-> 0.82-(0.10/39)*y,1:39);map(y-> 0.72+(0.22/15)*y,1:15);map(y-> 0.94-(0.149/67)*y,1:67);map(y-> 0.791-(0.31/27)*y,1:27);map(y-> 0.481-(0.14/18)*y,1:18)]
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,9,rc,dc,518,true,425,70,0)
run_param_scen_cal(true,0.121,"usa",30,1,1,1,1,1,10,rc,dc,518,false,425,70,0)



#for the national estimations
dc = [1;map(y->48+y,0:6);map(y->84+y,0:11);map(y->127+y,0:21);map(y->177+y,0:27);map(y-> 222+y,0:38);map(y->280+y,0:13);map(y->330+y,0:66);map(y->442+y,0:26);map(y->488+y,0:19)]
rc = [1.0;map(y->1.0+(0.06/7)*y,1:7);map(y->1.06-(0.1/12)*y,1:12);map(y->0.96-(0.16/22)*y,1:22);map(y-> 0.8+(0.075/28)*y,1:28);map(y-> 0.875-(0.12/39)*y,1:39);map(y-> 0.755+(0.25/14)*y,1:14);map(y-> 1.005-(0.15/67)*y,1:67);map(y-> 0.855-(0.34/27)*y,1:27);map(y-> 0.505-(0.155556/20)*y,1:20)]
run_param_scen_cal(true,0.096,"usa",50,1,1,1,1,1,1,rc,dc,528,true,425,70,1) 
run_param_scen_cal(true,0.096,"usa",50,1,1,1,1,1,2,rc,dc,528,false,425,70,1) 

