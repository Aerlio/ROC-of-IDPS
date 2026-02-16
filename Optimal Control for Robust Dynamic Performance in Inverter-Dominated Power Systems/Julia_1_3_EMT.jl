global time_NLP = time() - time_NLP_start ;
if (phase=="Max")
    global call_name = string(Optim_iter)*"_"*phase*"_"*string(ms_count)

    # ***********************************************************************
    # ******************* Save tmp files and run Simulink *******************
    # ***********************************************************************

    # --------------------- Run Matlab ----------------------
    save_to_matlab( "Generated_files//tmp//tmp_1.mat",
                    Vars_res[1:1,:] ,
                    y_theta        ,   n1_ratio        ,   n2_ratio       ,
                    Pload_ini_res  ,   Pload_fin_res   ,   Qload_ini_res  ,   Qload_fin_res  ,   p_ref_ed_res   ,   q_ref_ed_res   ,   p_ref_agc_res  )
    println("Running Simulink Simulation")
    sleep(3)
    eng.eval("run('Matlab_2_Ref_Loop.m')", nargout=0)
    sleep(3)
    # --------------- Read simulation results ---------------
    global results = matopen("Generated_files//tmp//"*phase*"_"*string(Optim_iter)*"_tmp.mat");
    (ts_EMT, Vars_EMT , EMT_obj, EMT_time) = [read(results,i) for i in ["ts_EMT",  "Vars_EMT", "EMT_obj", "EMT_time"]];
    close(results)

    # ------------------ Save simulation ---------------------
    global EMT_results = matopen(BASE_SAVE_FOLDER*"//EMT//"*call_name*".mat", "w")
    write(EMT_results, "ts_EMT"    , ts_EMT    );
    write(EMT_results, "Vars_EMT"  , Vars_EMT  );
    write(EMT_results, "EMT_obj"   , EMT_obj   );
    write(EMT_results, "EMT_time"  , EMT_time  );
    close(EMT_results)

    # ----------------- Save for plotting --------------------
    global Opt_results = matopen(BASE_SAVE_FOLDER*"//Optimization//"*call_name*"_for_plotting.mat", "w")
    write(Opt_results, "p_ref_res"  , p_ref_res );
    write(Opt_results, "q_ref_res"  , q_ref_res );
    write(Opt_results, "p_dist_res" , p_dist_res);
    write(Opt_results, "q_dist_res" , q_dist_res);
    write(Opt_results, "Pibr_res"   , Pibr_res  );  
    write(Opt_results, "Qibr_res"   , Qibr_res  );  
    write(Opt_results, "pload_res"  , pload_res ); 
    write(Opt_results, "qload_res"  , qload_res ); 
    write(Opt_results, "w_dq_res"   , w_dq_res  ); 
    write(Opt_results, "ixmag_res"  , ixmag_res );
    write(Opt_results, "vmag_res"   , vmag_res  ); 
    write(Opt_results, "score_res"  , score_res );
    close(Opt_results)

    # ----------------- Save optimization --------------------
    save_results(BASE_SAVE_FOLDER*"//Optimization//"*call_name, 
                                            break_times       , t_p_sec           , Vars_res            , Derivs_res        ,
                                            Pload_ini_res     , Pload_fin_res     , Qload_ini_res       , Qload_fin_res     ,  p_ref_ed_res      , q_ref_ed_res      , p_ref_agc_res     , 
                                            ibr_conf_num      , y_theta           , n1_ratio            , n2_ratio          , 
                                            OPT_obj_res , time_NLP , EMT_obj , EMT_time)


    # ***********************************************************************
    # ***********************************************************************
    # ***********************************************************************
    
    println("Simulation readed. Accuracy between Sim & Opt: ", 100*(abs.(EMT_obj-OPT_obj_res)/(EMT_obj)), " %.")
    println("Optim. score: ",round.( OPT_obj_res , digits=4 ) )
    println("Sim.   score: ",round.( EMT_obj     , digits=4 ) )

    # ***********************************************************************
    # ***************** Save as best solution if applicable******************
    # ***********************************************************************

    if (EMT_obj >= global_best_score) 
        global best_solution        = true          ; 
        global global_best_score    = copy(EMT_obj)  ;

        save_results(BASE_SAVE_FOLDER*"//Optimization//Selected_Results//"*string(Optim_iter)*"_"*phase, 
                                            break_times       , t_p_sec           , Vars_res            , Derivs_res        ,
                                            Pload_ini_res     , Pload_fin_res     , Qload_ini_res       , Qload_fin_res     ,  p_ref_ed_res      , q_ref_ed_res      , p_ref_agc_res     , 
                                            ibr_conf_num      , y_theta           , n1_ratio            , n2_ratio          , 
                                            OPT_obj_res , time_NLP , EMT_obj , EMT_time)
    else
        global best_solution = false ; 
    end

    include("Julia_1_4_Plot_Series.jl")  # Plot results
    if best_solution savefig(plot_emt,BASE_SAVE_FOLDER*"//Plots//Selected_Results_Plots//"*string(Optim_iter)*"_"*phase*".svg") end
    global (ts_EMT,Vars_EMT) = (0.0 , 0.0)  ; # avoid memory issues
    global plot_emt = 0                     ; # avoid memory issues

elseif (phase=="Min")

    global call_name = string(Optim_iter)*"_"*phase*"_"*string(ID_configs[ibr_conf_num])
    # ***********************************************************************
    # ******************* Save tmp files and run Simulink *******************
    # ***********************************************************************

    # --------------------- Prepare files for Matlab ----------------------
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

    # --------------------- Run Matlab ----------------------
    println("Running Simulink Simulation")
    sleep(3)
    eng.eval("run('Matlab_2_Ref_Loop.m')", nargout=0)
    sleep(3)
    # --------------- Read simulation results ---------------
    global results = matopen("Generated_files//tmp//"*phase*"_"*string(Optim_iter)*"_tmp.mat");
    (ts_EMT, Vars_EMT , EMT_obj, EMT_time) = [read(results,i) for i in ["ts_EMT",  "Vars_EMT", "EMT_obj", "EMT_time"]];
    close(results)
    # ------------------ Save simulation ---------------------
    global EMT_results = matopen(BASE_SAVE_FOLDER*"//EMT//"*call_name*".mat", "w")
    write(EMT_results, "ts_EMT"    , ts_EMT    );
    write(EMT_results, "Vars_EMT"  , Vars_EMT  );
    write(EMT_results, "EMT_obj"   , EMT_obj   );
    write(EMT_results, "EMT_time"  , EMT_time  );
    close(EMT_results)
    
    # ----------------- Save for plotting --------------------
    global Opt_results = matopen(BASE_SAVE_FOLDER*"//Optimization//"*call_name*"_for_plotting.mat", "w")
    write(Opt_results, "p_ref_res"  , p_ref_res );
    write(Opt_results, "q_ref_res"  , q_ref_res );
    write(Opt_results, "p_dist_res" , p_dist_res);
    write(Opt_results, "q_dist_res" , q_dist_res);
    write(Opt_results, "Pibr_res"   , Pibr_res  );  
    write(Opt_results, "Qibr_res"   , Qibr_res  );  
    write(Opt_results, "pload_res"  , pload_res ); 
    write(Opt_results, "qload_res"  , qload_res ); 
    write(Opt_results, "w_dq_res"   , w_dq_res  ); 
    write(Opt_results, "ixmag_res"  , ixmag_res );
    write(Opt_results, "vmag_res"   , vmag_res  ); 
    write(Opt_results, "score_res"  , score_res );
    close(Opt_results)

    # ----------------- Save optimization --------------------
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
        save_results(BASE_SAVE_FOLDER*"//Optimization//"*string(Optim_iter)*"_"*phase*"_mnv_"*string(optim_iter)*"_"*string(ID_configs[ibr_conf_num]), 
                                break_times_mnv    , t_p_sec_mnv   , Vars_mnv      , Derivs_mnv    ,
                                Pload_ini_mnv      , Pload_fin_mnv , Qload_ini_mnv , Qload_fin_mnv ,  p_ref_ed_mnv , q_ref_ed_mnv  , p_ref_agc_mnv ,
                                ibr_conf_num_res   , y_theta_res   , n1_ratio_res  , n2_ratio_res  , 
                                OPT_obj_res , time_NLP , EMT_obj , EMT_time)
    end
    
    # ***********************************************************************
    # ***********************************************************************
    # ***********************************************************************

    println("Simulation readed. Accuracy between Sim & Opt: ", 100*(abs.(EMT_obj-OPT_obj_res)/(EMT_obj)), " %.")
    println("Optim. score: ",round.( OPT_obj_res , digits=4 ) )
    println("Sim.   score: ",round.( EMT_obj , digits=4 ) )
    
    # ***********************************************************************
    # ***************** Save as best solution if applicable******************
    # ***********************************************************************

    if (EMT_obj <= global_best_score) 
        global best_solution        = true           ; 
        global global_best_score    = copy(EMT_obj)  ;

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
            save_results(BASE_SAVE_FOLDER*"//Optimization//Selected_Results//"*string(Optim_iter)*"_"*phase*"_mnv_"*string(optim_iter), 
                                    break_times_mnv    , t_p_sec_mnv   , Vars_mnv      , Derivs_mnv    ,
                                    Pload_ini_mnv      , Pload_fin_mnv , Qload_ini_mnv , Qload_fin_mnv ,  p_ref_ed_mnv , q_ref_ed_mnv  , p_ref_agc_mnv ,
                                    ibr_conf_num_res   , y_theta_res   , n1_ratio_res  , n2_ratio_res  , 
                                    OPT_obj_res , time_NLP , EMT_obj , EMT_time)
        end

    else
        global best_solution = false ; 
    end

    
    include("Julia_1_4_Plot_Series.jl")  # Plot results
    if best_solution savefig(plot_emt,BASE_SAVE_FOLDER*"//Plots//Selected_Results_Plots//"*string(Optim_iter)*"_"*phase*".svg") end
    global (ts_EMT,Vars_EMT) = (0.0 , 0.0)  ; # avoid memory issues
    global plot_emt = 0                     ; # avoid memory issues
end
