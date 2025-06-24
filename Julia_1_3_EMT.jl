# ***************************************************
# *************** Run EMT simulation ****************
# ***************************************************
if (phase=="Max")
    # ******** Save files to read in matlab *********
    save_to_matlab( "Generated_files//tmp//tmp_1.mat",
                    Vars_res[1:1,:] ,
                    y_theta        ,   n1_ratio        ,   n2_ratio       ,
                    Pload_ini_res  ,   Pload_fin_res   ,   Qload_ini_res  ,   Qload_fin_res  ,   p_ref_ed_res   ,   q_ref_ed_res   ,   p_ref_agc_res  )

    # --------------- Run Matlab and read results --------------
    println("Running Simulink Simulation")
    sleep(3)
    eng.eval("run('Matlab_2_Ref_Loop.m')", nargout=0)
    sleep(3)
    global results = matopen("Generated_files//tmp//"*phase*"_"*string(Optim_iter)*"_tmp.mat");
    (Time_ML, Vars_ML , ML_obj) = [read(results,i) for i in ["Time_ML",  "Vars_ML", "ML_obj"]];
    close(results)
    println("Simulation readed. Accuracy between Sim & Opt: ", 100*(abs.(ML_obj-JL_obj_res)/(ML_obj)), " %.")
    include("1_4_Plot_series.jl")
    println("Optim. score: ",round.( JL_obj_res , digits=4 ) )
    println("Sim.   score: ",round.( ML_obj     , digits=4 ) )

    
    if (ML_obj >= global_best_score) 
        global best_solution        = true          ; 
        global global_best_score    = copy(ML_obj)  ;
        savefig(plot_emt,BASE_SAVE_FOLDER*"//1_Mesh_Refinement//2_best_plots//"*string(Optim_iter)*"_"*phase*".svg")
    else
        global best_solution = false ; 
    end

    # ------------------ Save ws solution -----------------------
    save_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*string(Optim_iter)*"_"*phase*"_"*string(ms_count), 
                                            break_times       , t_p_sec           , Vars_res            , Derivs_res        ,
                                            Pload_ini_res     , Pload_fin_res     , Qload_ini_res       , Qload_fin_res     ,  p_ref_ed_res      , q_ref_ed_res      , p_ref_agc_res     , 
                                            ibr_conf_num      , y_theta           , n1_ratio            , n2_ratio          , copy(ML_obj))
    # -------------- Save best global solution ------------------
    if (best_solution==true)
        save_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*string(Optim_iter)*"_"*phase, 
                                                break_times       , t_p_sec           , Vars_res            , Derivs_res        ,
                                                Pload_ini_res     , Pload_fin_res     , Qload_ini_res       , Qload_fin_res     ,  p_ref_ed_res      , q_ref_ed_res      , p_ref_agc_res     , 
                                                ibr_conf_num      , y_theta           , n1_ratio            , n2_ratio          , global_best_score)
    end

    global (Time_ML,Vars_ML) = (0.0 , 0.0)  ; # avoid memory issues
    global plot_emt = 0                     ; # avoid memory issues
elseif (phase=="Min")
    # ******** Save files to read in matlab *********
    for optim_iter in 1:Optim_iter
        global Vars_mnv        = copy(Vars_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])  ;
        global Pload_ini_mnv   = copy(M_Pload_ini[:,optim_iter] )                           ;
        global Qload_ini_mnv   = copy(M_Qload_ini[:,optim_iter] )                           ;
        global Pload_fin_mnv   = copy(M_Pload_fin[:,optim_iter] )                           ;
        global Qload_fin_mnv   = copy(M_Qload_fin[:,optim_iter] )                           ;
        global p_ref_ed_mnv    = copy(M_p_ref_ed[:,optim_iter]  )                           ;
        global q_ref_ed_mnv    = copy(M_q_ref_ed[:,optim_iter]  )                           ;
        global p_ref_agc_mnv   = copy(M_p_ref_agc[:,optim_iter] )                           ;
        global tmp_file_name = "Generated_files//tmp//tmp_"*string(optim_iter)*".mat"
        save_to_matlab( tmp_file_name,
                        Vars_mnv[1:1,:] ,
                        y_theta_res    ,   n1_ratio_res    ,   n2_ratio_res   ,
                        Pload_ini_mnv  ,   Pload_fin_mnv   ,   Qload_ini_mnv  ,   Qload_fin_mnv  ,   p_ref_ed_mnv   ,   q_ref_ed_mnv   ,   p_ref_agc_mnv  )
    end

    # --------------- Run Matlab and read results --------------
    println("Running Simulink Simulation")
    sleep(3)
    eng.eval("run('Matlab_2_Ref_Loop.m')", nargout=0)
    sleep(3)
    global results = matopen("Generated_files//tmp//"*phase*"_"*string(Optim_iter)*"_tmp.mat");
    (Time_ML, Vars_ML , ML_obj) = [read(results,i) for i in ["Time_ML",  "Vars_ML", "ML_obj"]];
    close(results)
    println("Simulation readed. Accuracy between Sim & Opt: ", 100*(abs.(ML_obj-JL_obj_res)/(ML_obj)), " %.")
    include("1_4_Plot_Series.jl")
    println("Optim. score: ",round.( JL_obj_res , digits=4 ) )
    println("Sim.   score: ",round.( ML_obj , digits=4 ) )
    
    if (ML_obj <= global_best_score) 
        global best_solution        = true          ; 
        global global_best_score    = copy(ML_obj)  ;
        savefig(plot_emt,BASE_SAVE_FOLDER*"//1_Mesh_Refinement//2_best_plots//"*string(Optim_iter)*"_"*phase*"_"*string(ID_configs[ibr_conf_num])*".svg")
    else
        global best_solution = false ; 
    end

    # -------------------------- Save config solution ------------------------
    for optim_iter in 1:Optim_iter
        global Vars_mnv        = copy(Vars_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])  ;
        global Derivs_mnv      = copy(Derivs_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:]);
        global Pload_ini_mnv   = copy(M_Pload_ini[:,optim_iter] )                           ;
        global Qload_ini_mnv   = copy(M_Qload_ini[:,optim_iter] )                           ;
        global Pload_fin_mnv   = copy(M_Pload_fin[:,optim_iter] )                           ;
        global Qload_fin_mnv   = copy(M_Qload_fin[:,optim_iter] )                           ;
        global p_ref_ed_mnv    = copy(M_p_ref_ed[:,optim_iter]  )                           ;
        global q_ref_ed_mnv    = copy(M_q_ref_ed[:,optim_iter]  )                           ;
        global p_ref_agc_mnv   = copy(M_p_ref_agc[:,optim_iter] )                           ;
        save_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*string(Optim_iter)*"_"*phase*"_mnv_"*string(optim_iter)*"_"*string(ID_configs[ibr_conf_num]), 
                                break_times_mnv    , t_p_sec_mnv   , Vars_mnv      , Derivs_mnv    ,
                                Pload_ini_mnv      , Pload_fin_mnv , Qload_ini_mnv , Qload_fin_mnv ,  p_ref_ed_mnv , q_ref_ed_mnv  , p_ref_agc_mnv ,
                                ibr_conf_num_res   , y_theta_res   , n1_ratio_res  , n2_ratio_res  , copy(ML_obj))
    # ------------ Save best GLOBAL (ACROSS CONFIGS) solution ------------------
        if (best_solution == true)
            save_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*string(Optim_iter)*"_"*phase*"_mnv_"*string(optim_iter), 
                                    break_times_mnv    , t_p_sec_mnv   , Vars_mnv      , Derivs_mnv    ,
                                    Pload_ini_mnv      , Pload_fin_mnv , Qload_ini_mnv , Qload_fin_mnv ,  p_ref_ed_mnv , q_ref_ed_mnv  , p_ref_agc_mnv ,
                                    ibr_conf_num_res   , y_theta_res   , n1_ratio_res  , n2_ratio_res  , global_best_score)
        end
    end

    global (Time_ML,Vars_ML) = (0.0 , 0.0)  ; # avoid memory issues
    global plot_emt = 0                     ; # avoid memory issues
end
