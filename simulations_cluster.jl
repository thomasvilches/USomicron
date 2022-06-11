using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using ClusterManagers
using Dates
using DelimitedFiles

## load the packages by covid19abm

#using covid19abm

#addprocs(2, exeflags="--project=.")


#@everywhere using covid19abm

addprocs(SlurmManager(250), N=16, topology=:master_worker, exeflags = "--project=.")
@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere include("covid19abm.jl")
@everywhere const cv=covid19abm


function run(myp::cv.ModelParameters, nsims=1000, folderprefix="./")
    println("starting $nsims simulations...\nsave folder set to $(folderprefix)")
    dump(myp)
   
    # will return 6 dataframes. 1 total, 4 age-specific 
    cdr = pmap(1:nsims) do x                 
            cv.runsim(x, myp)
    end      

    println("simulations finished")
    println("total size of simulation dataframes: $(Base.summarysize(cdr))")
    ## write the infectors     

    ## write contact numbers
    #writedlm("$(folderprefix)/ctnumbers.dat", [cdr[i].ct_numbers for i = 1:nsims])    
    ## stack the sims together
    allag = vcat([cdr[i].a  for i = 1:nsims]...)
    working = vcat([cdr[i].work for i = 1:nsims]...)
    ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cdr[i].g3 for i = 1:nsims]...)
    ag4 = vcat([cdr[i].g4 for i = 1:nsims]...)
    ag5 = vcat([cdr[i].g5 for i = 1:nsims]...)
    ag6 = vcat([cdr[i].g6 for i = 1:nsims]...)
    ag7 = vcat([cdr[i].g7 for i = 1:nsims]...)
    ag8 = vcat([cdr[i].g8 for i = 1:nsims]...) 
    ag9 = vcat([cdr[i].g9 for i = 1:nsims]...)  

    mydfs = Dict("all" => allag, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4, "ag5" => ag5, "ag6" => ag6,"ag7" => ag7,"ag8" => ag8,"ag9" => ag9, "working"=>working)
    #mydfs = Dict("all" => allag, "working"=>working, "kids"=>kids)
    #mydfs = Dict("all" => allag)
    
    ## save at the simulation and time level
    ## to ignore for now: miso, iiso, mild 
    #c1 = Symbol.((:LAT, :ASYMP, :INF, :PRE, :MILD,:IISO, :HOS, :ICU, :DED), :_INC)
    #c2 = Symbol.((:LAT, :ASYMP, :INF, :PRE, :MILD,:IISO, :HOS, :ICU, :DED), :_PREV)
    
    c1 = Symbol.((:LAT, :ASYMP, :PRE, :MILD, :INF, :HOS, :ICU, :DED,:LAT2, :ASYMP2,:PRE2, :MILD2, :INF2, :HOS2, :ICU2, :DED2,:LAT3, :ASYMP3, :PRE3, :MILD3, :INF3, :HOS3, :ICU3, :DED3,:LAT4, :PRE4, :ASYMP4, :MILD4, :INF4, :HOS4, :ICU4, :DED4,:LAT5, :ASYMP5, :PRE5, :MILD5, :INF5, :HOS5, :ICU5, :DED5,:LAT6, :PRE6, :ASYMP6, :MILD6, :INF6, :HOS6, :ICU6, :DED6), :_INC)
    #c2 = Symbol.((:LAT, :MILD, :INF, :HOS, :ICU, :DED,:LAT2, :MILD2, :INF2, :HOS2, :ICU2, :DED2,:LAT3, :MILD3, :INF3, :HOS3, :ICU3, :DED3,:LAT4, :MILD4, :INF4, :HOS4, :ICU4, :DED4,:LAT5, :MILD5, :INF5, :HOS5, :ICU5, :DED5,:LAT6, :MILD6, :INF6, :HOS6, :ICU6, :DED6), :_PREV)
    #c2 = Symbol.((:HOS, :ICU,:HOS2, :ICU2, :HOS3, :ICU3, :HOS4, :ICU4, :HOS5, :ICU5, :HOS6, :ICU6), :_PREV)
    
    #c2 = Symbol.((:LAT, :HOS, :ICU, :DED,:LAT2, :HOS2, :ICU2, :DED2), :_PREV)
    for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        # simulation level, save file per health status, per age group
        #for c in vcat(c1..., c2...)
        for c in vcat(c1...)
        #for c in vcat(c2...)
            udf = unstack(df, :time, :sim, c) 
            fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
            CSV.write(fn, udf)
        end
        println("saving dataframe time level: $k")
        # time level, save file per age group
        #yaf = compute_yearly_average(df)       
        #fn = string("$(folderprefix)/timelevel_", k, ".dat")   
        #CSV.write(fn, yaf)       
    end
    
    cov1 = [cdr[i].cov1 for i=1:nsims]
    cov2 = [cdr[i].cov2 for i=1:nsims]
    cov12 = [cdr[i].cov12 for i=1:nsims]
    cov22 = [cdr[i].cov22 for i=1:nsims]
    
    writedlm(string(folderprefix,"/cov.dat"),[cov1 cov2 cov12 cov22])

    vac_m = [cdr[i].n_moderna for i=1:nsims]
    vac_m_w = [cdr[i].n_moderna_w for i=1:nsims]

    vac_p = [cdr[i].n_pfizer for i=1:nsims]
    vac_p_w = [cdr[i].n_pfizer_w for i=1:nsims]

    vac_j = [cdr[i].n_jensen for i=1:nsims]
    vac_j_w = [cdr[i].n_jensen_w for i=1:nsims]

    vac_m_2 = [cdr[i].n_moderna_2 for i=1:nsims]
    vac_m_w_2 = [cdr[i].n_moderna_w_2 for i=1:nsims]

    vac_p_2 = [cdr[i].n_pfizer_2 for i=1:nsims]
    vac_p_w_2 = [cdr[i].n_pfizer_w_2 for i=1:nsims]

    vac_j_2 = [cdr[i].n_jensen_2 for i=1:nsims]
    vac_j_w_2 = [cdr[i].n_jensen_w_2 for i=1:nsims]

    vac_m_3 = [cdr[i].n_moderna_3 for i=1:nsims]
    vac_m_w_3 = [cdr[i].n_moderna_w_3 for i=1:nsims]

    vac_p_3 = [cdr[i].n_pfizer_3 for i=1:nsims]
    vac_p_w_3 = [cdr[i].n_pfizer_w_3 for i=1:nsims]

    vac_j_3 = [cdr[i].n_jensen_3 for i=1:nsims]
    vac_j_w_3 = [cdr[i].n_jensen_w_3 for i=1:nsims]

    vac_m_4 = [cdr[i].n_moderna_4 for i=1:nsims]
    vac_m_w_4 = [cdr[i].n_moderna_w_4 for i=1:nsims]

    vac_p_4 = [cdr[i].n_pfizer_4 for i=1:nsims]
    vac_p_w_4 = [cdr[i].n_pfizer_w_4 for i=1:nsims]

    vac_j_4 = [cdr[i].n_jensen_4 for i=1:nsims]
    vac_j_w_4 = [cdr[i].n_jensen_w_4 for i=1:nsims]

    nvacgiven = [cdr[i].nvacgiven for i=1:nsims]
    n5plus = [cdr[i].n5plus for i=1:nsims]



    writedlm(string(folderprefix,"/vaccine_all.dat"),[vac_p vac_m vac_j vac_p_2 vac_m_2 vac_j_2 vac_p_3 vac_m_3 vac_j_3 vac_p_4 vac_m_4 vac_j_4])
    writedlm(string(folderprefix,"/vaccine_working.dat"),[vac_p_w vac_m_w vac_j_w vac_p_w_2 vac_m_w_2 vac_j_w_2 vac_p_w_3 vac_m_w_3 vac_j_w_3 vac_p_w_4 vac_m_w_4 vac_j_w_4])

    writedlm(string(folderprefix,"/year_of_death.dat"),hcat([cdr[i].vector_dead for i=1:nsims]...))
    writedlm(string(folderprefix,"/nvac_given.dat"),[nvacgiven n5plus])


    writedlm(string(folderprefix,"/lat_vac_1.dat"),hcat([cdr[i].lat for i=1:nsims]...))
    writedlm(string(folderprefix,"/lat_vac_2.dat"),hcat([cdr[i].lat2 for i=1:nsims]...))
    writedlm(string(folderprefix,"/lat_vac_3.dat"),hcat([cdr[i].lat3 for i=1:nsims]...))
    writedlm(string(folderprefix,"/lat_unvac_r.dat"),hcat([cdr[i].lat4 for i=1:nsims]...))
    writedlm(string(folderprefix,"/lat_unvac_nr.dat"),hcat([cdr[i].lat5 for i=1:nsims]...))
    writedlm(string(folderprefix,"/lat_vac_1_r.dat"),hcat([cdr[i].lat6 for i=1:nsims]...))
    writedlm(string(folderprefix,"/lat_vac_2_r.dat"),hcat([cdr[i].lat7 for i=1:nsims]...))
    writedlm(string(folderprefix,"/lat_vac_3_r.dat"),hcat([cdr[i].lat8 for i=1:nsims]...))

    writedlm(string(folderprefix,"/hos_vac_1.dat"),hcat([cdr[i].hos for i=1:nsims]...))
    writedlm(string(folderprefix,"/hos_vac_2.dat"),hcat([cdr[i].hos2 for i=1:nsims]...))
    writedlm(string(folderprefix,"/hos_vac_3.dat"),hcat([cdr[i].hos3 for i=1:nsims]...))
    writedlm(string(folderprefix,"/hos_unvac_r.dat"),hcat([cdr[i].hos4 for i=1:nsims]...))
    writedlm(string(folderprefix,"/hos_unvac_nr.dat"),hcat([cdr[i].hos5 for i=1:nsims]...))
    writedlm(string(folderprefix,"/hos_vac_1_r.dat"),hcat([cdr[i].hos6 for i=1:nsims]...))
    writedlm(string(folderprefix,"/hos_vac_2_r.dat"),hcat([cdr[i].hos7 for i=1:nsims]...))
    writedlm(string(folderprefix,"/hos_vac_3_r.dat"),hcat([cdr[i].hos8 for i=1:nsims]...))

    writedlm(string(folderprefix,"/icu_vac_1.dat"),hcat([cdr[i].icu for i=1:nsims]...))
    writedlm(string(folderprefix,"/icu_vac_2.dat"),hcat([cdr[i].icu2 for i=1:nsims]...))
    writedlm(string(folderprefix,"/icu_vac_3.dat"),hcat([cdr[i].icu3 for i=1:nsims]...))
    writedlm(string(folderprefix,"/icu_unvac_r.dat"),hcat([cdr[i].icu4 for i=1:nsims]...))
    writedlm(string(folderprefix,"/icu_unvac_nr.dat"),hcat([cdr[i].icu5 for i=1:nsims]...))
    writedlm(string(folderprefix,"/icu_vac_1_r.dat"),hcat([cdr[i].icu6 for i=1:nsims]...))
    writedlm(string(folderprefix,"/icu_vac_2_r.dat"),hcat([cdr[i].icu7 for i=1:nsims]...))
    writedlm(string(folderprefix,"/icu_vac_3_r.dat"),hcat([cdr[i].icu8 for i=1:nsims]...))

    writedlm(string(folderprefix,"/ded_vac_1.dat"),hcat([cdr[i].ded for i=1:nsims]...))
    writedlm(string(folderprefix,"/ded_vac_2.dat"),hcat([cdr[i].ded2 for i=1:nsims]...))
    writedlm(string(folderprefix,"/ded_vac_3.dat"),hcat([cdr[i].ded3 for i=1:nsims]...))
    writedlm(string(folderprefix,"/ded_unvac_r.dat"),hcat([cdr[i].ded4 for i=1:nsims]...))
    writedlm(string(folderprefix,"/ded_unvac_nr.dat"),hcat([cdr[i].ded5 for i=1:nsims]...))
    writedlm(string(folderprefix,"/ded_vac_1_r.dat"),hcat([cdr[i].ded6 for i=1:nsims]...))
    writedlm(string(folderprefix,"/ded_vac_2_r.dat"),hcat([cdr[i].ded7 for i=1:nsims]...))
    writedlm(string(folderprefix,"/ded_vac_3_r.dat"),hcat([cdr[i].ded8 for i=1:nsims]...))

    return mydfs
end


function create_folder(ip::cv.ModelParameters,province="usa",calibrating = true,letter = "A")
    
    #RF = string("heatmap/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vaccine_ef), "." => "_"))","_herd_immu_","$(ip.herd)","_$strategy","cov_$(replace(string(ip.cov_val)))") ## 
    main_folder = "/data/thomas-covid/USbooster_scenarios"
    #main_folder = "."
    if calibrating
        RF = string(main_folder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_$(ip.file_index)_$(letter)_$(province)") ##  
    else
        RF = string(main_folder,"/results_prob_$(ip.file_index)_$(province)") ##  
    end
    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end



function run_param_scen_cal(calibrating::Bool,b::Float64,province::String="usa",ic1::Int64=1,ic2::Int64=1,ic3::Int64=1,ic4::Int64=1,ic5::Int64=1,ic6::Int64=1,index::Int64 = 0,idxtime::Int64 = 1,rc=[0.0],dc=[0],mt::Vector{Int64}=[973;-1;-3;-5;-7;-8],vac_cov::Vector{Float64}=[0.8;0.8;0.8;0.8;0.8],vac::Bool=true,tbn::Int64 = 9999,ro::Int64 = 1,dr::Int64=0,hospar::Float64 = 3.1,nsims::Int64=500)
    
    letters = ["A";"B";"C";"D"]
    #b = bd[h_i]
    #ic = init_con[h_i]
    @everywhere ip = cv.ModelParameters(β=$b,fsevere = 1.0,fmild = 1.0,vaccinating = $vac,
    start_several_inf=true,initialinf3=$ic3,initialinf6=$ic6,initialinf=$ic1,initialinf2=$ic2,initialinf5=$ic5,initialinf4=$ic4,
    status_relax = 2, relax_after = 14,file_index = $index,
    modeltime=$mt, prov = Symbol($province),
    time_change_contact = $dc,
    change_rate_values = $rc,time_back_to_normal = $tbn,relax_over = $ro, reduce_days = $dr,
    hosp_red = $hospar,
    scenario = $index,
    intervention_prob = $vac_cov)

    folder = create_folder(ip,province,calibrating,letters[idxtime])

    #println("$v_e $(ip.vaccine_ef)")
    run(ip,nsims,folder)
   
end