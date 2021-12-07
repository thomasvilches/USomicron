
#############################

dc = [1;map(y-> 95+y,0:3);map(y->131+y,0:9);map(y->166+y,0:13);map(y->198+y,0:30);map(y->288+y,0:19);map(y->335+y,0:3)]#;map(y->98+y,0:22);map(y->195+y,0:4);map(y->224+y,0:59);map(y->298+y,0:19)]
rc = [1.0;map(y-> 1.0-(0.054/4)*y,1:4);map(y-> 0.946-(0.046/10)*y,1:10);map(y-> 0.90+(0.05/14)*y,1:14);map(y-> 0.95-(0.2945/31)*y,1:31);map(y-> 0.6555+(0.3945/20)*y,1:20);map(y-> 1.05-(0.165/4)*y,1:4)]#;map(y->1.0-(0.143/23)*y,1:23);map(y->0.857+(0.04/5)*y,1:5);map(y->0.897-(0.29/60)*y,1:60);map(y->0.607+(0.30/20)*y,1:20)]
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,1,3,14,rc,dc,435,true)
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,2,3,14,rc,dc,435,false)

#### vac change the behaviour
dc = [1;map(y-> 95+y,0:3);map(y->131+y,0:9);map(y->166+y,0:13);map(y->197+y,0:43);map(y->288+y,0:19);map(y->338+y,0:23)]#;map(y->288+y,0:19);map(y->335+y,0:3)]#;map(y->98+y,0:22);map(y->195+y,0:4);map(y->224+y,0:59);map(y->298+y,0:19)]
rc = [1.0;map(y-> 1.0-(0.054/4)*y,1:4);map(y-> 0.946-(0.046/10)*y,1:10);map(y-> 0.90+(0.05/14)*y,1:14);map(y-> 0.95-(0.418/44)*y,1:44);map(y-> 0.532+(0.32/20)*y,1:20);map(y-> 0.852-(0.32/24)*y,1:24)]#;map(y-> 1.05-(0.165/4)*y,1:4)]#;map(y->1.0-(0.143/23)*y,1:23);map(y->0.857+(0.04/5)*y,1:5);map(y->0.897-(0.29/60)*y,1:60);map(y->0.607+(0.30/20)*y,1:20)]
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,5,2,14,rc,dc,435,true)
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,6,2,14,rc,dc,435,false)


#### nonvac change the behaviour slowly
dc = [1;map(y-> 95+y,0:3);map(y->131+y,0:9);map(y->166+y,0:13);map(y->200+y,0:44);map(y->300+y,0:19);map(y->339+y,0:9)]#;map(y->288+y,0:19);map(y->338+y,0:23)]#;map(y->288+y,0:19);map(y->335+y,0:3)]#;map(y->98+y,0:22);map(y->195+y,0:4);map(y->224+y,0:59);map(y->298+y,0:19)]
rc = [1.0;map(y-> 1.0-(0.054/4)*y,1:4);map(y-> 0.946-(0.046/10)*y,1:10);map(y-> 0.90+(0.045/14)*y,1:14);map(y-> 0.945-(0.538/45)*y,1:45);map(y-> 0.407+(0.073/20)*y,1:20);map(y-> 0.48-(0.059/10)*y,1:10)]#;map(y-> 0.532+(0.32/20)*y,1:20);map(y-> 0.852-(0.32/24)*y,1:24)]#;map(y-> 1.05-(0.165/4)*y,1:4)]#;map(y->1.0-(0.143/23)*y,1:23);map(y->0.857+(0.04/5)*y,1:5);map(y->0.897-(0.29/60)*y,1:60);map(y->0.607+(0.30/20)*y,1:20)]
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,3,2,14,rc,dc,435,true,215,1)
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,4,2,14,rc,dc,435,false,215,1)


#### nonvac doesn't change the behaviour, vaccinated individuals go to a higher level, but oscilates with nonvac
dc = [1;map(y-> 95+y,0:3);map(y->131+y,0:9);map(y->166+y,0:13);map(y->199+y,0:31);map(y->295+y,0:14);map(y->339+y,0:9)]#;map(y->300+y,0:19);map(y->339+y,0:9)]#;map(y->288+y,0:19);map(y->338+y,0:23)]#;map(y->288+y,0:19);map(y->335+y,0:3)]#;map(y->98+y,0:22);map(y->195+y,0:4);map(y->224+y,0:59);map(y->298+y,0:19)]
rc = [1.0;map(y-> 1.0-(0.054/4)*y,1:4);map(y-> 0.946-(0.046/10)*y,1:10);map(y-> 0.90+(0.045/14)*y,1:14);map(y-> 0.945-(0.32/32)*y,1:32);map(y-> 0.625+(0.30/15)*y,1:15);map(y-> 0.925-(0.19/10)*y,1:10)]#;map(y-> 0.407+(0.073/20)*y,1:20);map(y-> 0.48-(0.059/10)*y,1:10)]#;map(y-> 0.532+(0.32/20)*y,1:20);map(y-> 0.852-(0.32/24)*y,1:24)]#;map(y-> 1.05-(0.165/4)*y,1:4)]#;map(y->1.0-(0.143/23)*y,1:23);map(y->0.857+(0.04/5)*y,1:5);map(y->0.897-(0.29/60)*y,1:60);map(y->0.607+(0.30/20)*y,1:20)]
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,7,2,14,rc,dc,435,true,999,1)
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,8,2,14,rc,dc,435,false,999,1)



#### nonvac doesn't change the behaviour, vaccinated individuals go to a higher level, but oscilates with nonvac, waning immunity! commit: 0d63f03a7
dc = [1;map(y-> 95+y,0:3);map(y->131+y,0:9);map(y->166+y,0:13);map(y->199+y,0:31);map(y->295+y,0:14);map(y->327+y,0:46)]#;map(y->300+y,0:19);map(y->339+y,0:9)]#;map(y->288+y,0:19);map(y->338+y,0:23)]#;map(y->288+y,0:19);map(y->335+y,0:3)]#;map(y->98+y,0:22);map(y->195+y,0:4);map(y->224+y,0:59);map(y->298+y,0:19)]
rc = [1.0;map(y-> 1.0-(0.054/4)*y,1:4);map(y-> 0.946-(0.046/10)*y,1:10);map(y-> 0.90+(0.045/14)*y,1:14);map(y-> 0.945-(0.32/32)*y,1:32);map(y-> 0.625+(0.25/15)*y,1:15);map(y-> 0.875-(0.2611/47)*y,1:47)]#;map(y-> 0.407+(0.073/20)*y,1:20);map(y-> 0.48-(0.059/10)*y,1:10)]#;map(y-> 0.532+(0.32/20)*y,1:20);map(y-> 0.852-(0.32/24)*y,1:24)]#;map(y-> 1.05-(0.165/4)*y,1:4)]#;map(y->1.0-(0.143/23)*y,1:23);map(y->0.857+(0.04/5)*y,1:5);map(y->0.897-(0.29/60)*y,1:60);map(y->0.607+(0.30/20)*y,1:20)]
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,9,2,14,rc,dc,435,true,999,1,1)
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,10,2,14,rc,dc,435,false,999,1,1)


#### nonvac doesn't change the behaviour, vaccinated individuals go to a higher level, but oscilates with nonvac, waning immunity both vaccine and rec
dc = [1;map(y-> 95+y,0:3);map(y->131+y,0:9);map(y->166+y,0:13);map(y->201+y,0:30);map(y->291+y,0:14);map(y->318+y,0:49)]#;map(y->300+y,0:19);map(y->339+y,0:9)]#;map(y->288+y,0:19);map(y->338+y,0:23)]#;map(y->288+y,0:19);map(y->335+y,0:3)]#;map(y->98+y,0:22);map(y->195+y,0:4);map(y->224+y,0:59);map(y->298+y,0:19)]
rc = [1.0;map(y-> 1.0-(0.054/4)*y,1:4);map(y-> 0.946-(0.046/10)*y,1:10);map(y-> 0.90+(0.045/14)*y,1:14);map(y-> 0.945-(0.31/31)*y,1:31);map(y-> 0.635+(0.28/15)*y,1:15);map(y-> 0.955-(0.31/50)*y,1:50)]#;map(y-> 0.407+(0.073/20)*y,1:20);map(y-> 0.48-(0.059/10)*y,1:10)]#;map(y-> 0.532+(0.32/20)*y,1:20);map(y-> 0.852-(0.32/24)*y,1:24)]#;map(y-> 1.05-(0.165/4)*y,1:4)]#;map(y->1.0-(0.143/23)*y,1:23);map(y->0.857+(0.04/5)*y,1:5);map(y->0.897-(0.29/60)*y,1:60);map(y->0.607+(0.30/20)*y,1:20)]
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,11,2,14,rc,dc,435,true,999,1,1)
run_param_scen_cal(true,0.119,"newyorkcity",20,20,1,5,1,1,1,125,999,173,7,999,12,2,14,rc,dc,435,false,999,1,1)





#Alpha 2020-03-15
#Gamma 999
#Delta 2020-08-15
#Iota 2020-08-21
#Beta 999
#### nonvac doesn't change the behaviour, vaccinated individuals go to a higher level, but oscilates with nonvac, waning immunity both vaccine and rec
dc = [1;map(y->75+y,0:(6));map(y->128+y,0:14);map(y->177+y,0:27);map(y-> 220+y,0:34);map(y->292+y,0:14);map(y->336+y,0:52);map(y->417+y,0:15)]#;map(y->122+y,0:15)]
rc = [1.0;map(y->1.0-(0.08/7)*y,1:7);map(y->0.92-(0.10/15)*y,1:15);map(y-> 0.82+(0.07/28)*y,1:28);map(y-> 0.89-(0.205/35)*y,1:35);map(y-> 0.685+(0.15/15)*y,1:15);map(y-> 0.835-(0.17/53)*y,1:53);map(y-> 0.665+(0.036/16)*y,1:16)]#;map(y->0.92-(0.085/16)*y,1:16)]
run_param_scen_cal(true,0.109,"usa",10,50,1,1,1,1,1,92,999,194,7,999,1,2,14,rc,dc,466,true,999,1,1)
run_param_scen_cal(true,0.109,"usa",10,50,1,1,1,1,1,92,999,194,7,999,2,2,14,rc,dc,466,false,999,1,1)




dc = [1;map(y->75+y,0:(6));map(y->128+y,0:14);map(y->177+y,0:27);map(y-> 220+y,0:34);map(y->292+y,0:14);map(y->336+y,0:52);map(y->417+y,0:15)]#;map(y->122+y,0:15)]
rc = [1.0;map(y->1.0-(0.08/7)*y,1:7);map(y->0.92-(0.10/15)*y,1:15);map(y-> 0.82+(0.07/28)*y,1:28);map(y-> 0.89-(0.205/35)*y,1:35);map(y-> 0.685+(0.15/15)*y,1:15);map(y-> 0.835-(0.17/53)*y,1:53);map(y-> 0.665+(0.036/16)*y,1:16)]#;map(y->0.92-(0.085/16)*y,1:16)]
run_param_scen_cal(false,0.109,"usa",10,50,1,1,1,1,1,92,999,194,7,446,3,2,14,rc,dc,852,true,999,1,1,0.0,1.0,false)
run_param_scen_cal(false,0.109,"usa",10,50,1,1,1,1,1,92,999,194,7,446,3,2,14,rc,dc,852,true,999,1,1,0.2,1.2,false)


dc = [1;map(y->75+y,0:(6));map(y->128+y,0:14);map(y->177+y,0:27);map(y-> 220+y,0:34);map(y->292+y,0:14);map(y->336+y,0:52);map(y->417+y,0:15)]#;map(y->122+y,0:15)]
rc = [1.0;map(y->1.0-(0.08/7)*y,1:7);map(y->0.92-(0.10/15)*y,1:15);map(y-> 0.82+(0.07/28)*y,1:28);map(y-> 0.89-(0.205/35)*y,1:35);map(y-> 0.685+(0.15/15)*y,1:15);map(y-> 0.835-(0.17/53)*y,1:53);map(y-> 0.665+(0.036/16)*y,1:16)]#;map(y->0.92-(0.085/16)*y,1:16)]
run_param_scen_cal(false,0.109,"usa",10,50,1,1,1,1,1,92,999,194,7,446,4,2,14,rc,dc,852,true,999,1,1,0.0,1.0,true)
run_param_scen_cal(false,0.109,"usa",10,50,1,1,1,1,1,92,999,194,7,446,4,2,14,rc,dc,852,true,999,1,1,0.2,1.2,true)



for r = 0.0:0.04:0.2, tt = 1.0:0.04:1.2
    run_param_scen_cal(false,0.109,"usa",10,50,1,1,1,1,1,92,999,194,7,446,3,2,14,rc,dc,852,true,999,1,1,r,tt,true)
end


for r = 0.0:0.04:0.2, tt = 1.0:0.04:1.2
    run_param_scen_cal(false,0.109,"usa",10,50,1,1,1,1,1,92,999,194,7,446,4,2,14,rc,dc,852,true,999,1,1,r,tt,false)
end
