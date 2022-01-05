# COVID-19 Agent Based Model
Impact of accelerating booster vaccination amidst Omicron surge in the United States

## Model details:
A stochastic, age-stratified agent-based computational model for the transmission dynamics of COVID-19. The computational model simulates autonomous agents (representing individuals in a human population) and their interactions within a constrained virtual environment. Agents follow the natural history of disease, including epidemiological stages of susceptible, infected and incubating, asymptomatic, presymptomatic, and symptomatic with either mild, severe, or critical illness, recovered, and dead.

Model features include:

- Age structured with realistic contact dynamics
- Asymptomatic, Presymptomatic transmission
- Isolation of mild/severe cases
- The average number of daily contacts can be changed to fit to data
- Four strains (Original, Alpha, Iota and Delta), corresponding to when they were identified in the United States
- Waning of immunity induced by the vaccine and recovery.

## How to download and run

Prerequisites: Julia 1.0.4, access to a cluster or a high-compute workstation. 

1) Download or **clone** the entire repository and navigate to the folder.
2) Launch Julia and cctivate the project by: `julia --project=.`. Double check if the project environment is set correct by entering `Pkg` mode by typing `]`. 
3) Instantiate the project by typing `] instantiate`.
4) Include the file *simulations_cluster.jl* using the command
```
include("simulations_cluster.jl")
```
Note, that in our version of this file we connect to our compute cluster using the `Slurm` cluster software. The user may want to simply use `addprocs` to run locally on their computer, run everything in a serial manner (takes long), or use a compute cluster with the help of `ClusterManagers`. The simulations/scenarios can be launched by executing 

```
run_param_scen_cal(calibrating,b,h,init1,init2,init3,init4,init5,init6,time2,time3,time4,time5,time6,index,status,days_after,v1,v2,modeltime,vaccination,when_relax,turnon,waning,red,trans,redred,nb,ddkids,rate_kids,boost_inc_day,boost_inc,change_elig,ba)
```

to run the scenarios. The arguments are

- calibrating: Boolean \- if *true*, the results's file name will not contain the reduction and trasmissibility of Omicron. Keep true.
- b: Float64 \- probability of transmission for presymptomatic cases.
- h: Int64 \- previous herd immunity in the population (either 5, 10, 20, 30 or 50%).
- init1-6: Int64 \- number of initial infected for strains.
- time1-6:Int64 \- day of introduction of strains. (In this model we are not using 3rd and 6th strains)
- vaccination: Bool \- apply vaccination?.
- index: Int64 \- index to differ different files (see *Model output*).
- status: Integer \- On April 3, vaccinated individuals are allowed to go back to the normal number of contacts. This argument controls if they need to have first or second dosis (set it to '1' or '2'). If you do not want to allow them to go back to normal, set it to '3'.
- days_after: Int64 \- how many days after the desire dosis the individuals are allowed to go back to normal behavior.
- v1: vector with changes in contact pattern.
- v2: vector when the change (v1) will happen.
- modeltime: Int64 \- simulated time (number of days).
- when_relax: Int64 \- time starting applying double dose in all age groups (keep it 999).
- turnon and waning: Int64 \- if vaccinated people have oscillations in contact pattern and if waning is on. Keep it at 1 and 1.
- red: Float64 \- Immunity escape of Omicron (0.8).
- trans: Float64 \- Transmissibility of Omicron compared to Delta (1.35).
- redred: Float64 \- Booster mitigating of reduction (0.75).
- nb: Int64 \- Number of booster per person.
- ddkids and rate_kids: Int64 and Float64 \- day to start increasing the doses for kids, the proportion for increase. Not used here.
- boost_inc_day: Int64 \- Day for booster rate increase.
- boost_inc: Float64 \- Increase proportion (2x, 3x).
- change_elig: Int64 \- Day of change in the eligibility for booster.
- ba: Vector{Int64} \- New eligibility (per vaccine).

The scenarios are in

```
include("scen.jl")
```

## Model output

***First, make sure that the address 'main_folder' inside function 'create_folder' in file [simulations_cluster](simulations_cluster.jl) points to a valid directory in your system.*** 

The model is parameterized to fit to data from US from Octuber to August. The [incidence data](data_us.csv) was taken from [NY times](https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv). The function 'run_param_scen_cal' will generate a folder, inside the pointed directory, named **results\__b_\_herd\_immu\__h_\__vaccine\__index_\__state_** in which all the ***variables*** are the ones cited in the previous section. ***b*** is the probability, but the '.' is replaced by '\_'. If you ran the ***scen.jl*** file.

Inside this folder, one will find different data files. The most important ones are the ones named
**simlevel_\*\_inc\_\*\*.dat** in which

- \* stands for **lat**, **lat2-6**, **hos**, **hos2-6**, **icu**, **icu2-6**, **ded** , **ded2-6**. Which are the number of infections, Non-ICU hospitalizations, ICU hospitalizations and deaths generated by each one of the strains. It may be necessary to scale the number of hospitalizations and deaths by the reported ones when analyzing the data.


The files contain a *modeltime* x *number of simulations* (by default 318x500) matrix. However, the first row is the heading of the file and the first column of it is the timeline. The other columns are the incidence of a given outcome in the given day of simulation.
