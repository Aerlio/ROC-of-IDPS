Random.seed!(4321*Optim_iter)

global phase            = "Min"                     ;
global H                = copy(H_min)               ; # Min phase time horizon
global ref_upd_time     = copy(ref_upd_time_min)    ; # Min phase AGC.upd.time

# --------- Loading previous solutions ---------
for optim_iter in 1:Optim_iter
    global phase_read = string(optim_iter)*"_Max"
    global (_               , _                 , mnv_Vars_res      , _                 ,
            mnv_Pload_ini   , mnv_Pload_fin     , mnv_Qload_ini     , mnv_Qload_fin     , mnv_p_ref_ed  , mnv_q_ref_ed  , mnv_p_ref_agc , 
            _               , mnv_y_theta       , mnv_n1_ratio      , mnv_n2_ratio      , mnv_max_emt_score ) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*phase_read);   
    if (optim_iter==1)
        global M_Vars_ini     = copy(reshape(copy(mnv_Vars_res[1,:]),1,Nx))             ;
        global M_Pload_ini    = copy(mnv_Pload_ini)                                     ;
        global M_Pload_fin    = copy(mnv_Pload_fin)                                     ;
        global M_Qload_ini    = copy(mnv_Qload_ini)                                     ;
        global M_Qload_fin    = copy(mnv_Qload_fin)                                     ;
        global M_p_ref_ed     = copy(mnv_p_ref_ed)                                      ;
        global M_q_ref_ed     = copy(mnv_q_ref_ed)                                      ;
        global M_p_ref_agc    = copy(mnv_p_ref_agc)                                     ;
    else
        global M_Vars_ini     = cat(M_Vars_ini  , copy(reshape(mnv_Vars_res[1,:],1,Nx)) , dims=1) ;
        global M_Pload_ini    = cat(M_Pload_ini , copy(mnv_Pload_ini)                   , dims=2) ;
        global M_Pload_fin    = cat(M_Pload_fin , copy(mnv_Pload_fin)                   , dims=2) ;
        global M_Qload_ini    = cat(M_Qload_ini , copy(mnv_Qload_ini)                   , dims=2) ;
        global M_Qload_fin    = cat(M_Qload_fin , copy(mnv_Qload_fin)                   , dims=2) ;
        global M_p_ref_ed     = cat(M_p_ref_ed  , copy(mnv_p_ref_ed)                    , dims=2) ;
        global M_q_ref_ed     = cat(M_q_ref_ed  , copy(mnv_q_ref_ed)                    , dims=2) ;
        global M_p_ref_agc    = cat(M_p_ref_agc , copy(mnv_p_ref_agc)                   , dims=2) ;
    end
    if (optim_iter == Optim_iter)
        global ws_n1_ratio = copy(mnv_n1_ratio) ;
        global ws_n2_ratio = copy(mnv_n2_ratio) ;
    end
    if optim_iter == 1  global Max_emt_score = mnv_max_emt_score;
    else                global Max_emt_score = min(mnv_max_emt_score , Max_emt_score); end
end


#********* Loop Initialization 1/2*********
global global_best_score        = 1e2   ; # Best score accross all configurations
global stop_criteria            = false ; # Criteria to stop looking for initial points (based on patience)

global any_feasible_found       = false ; # Indicates if any feasible solution has been found in entire min. phase
global local_feasible_found     = false ; # Indicates if any feasible solution has been found for given config.
global mesh_feasible_found      = false ; # Indicates if any feasible solution has been found for given config. and mesh

#********** Look into configurations that achieved feasibility with last maximisation phase scenarios *************
global keep_hope_for_configs = false;
if (keep_hope_for_configs == true)
    global Configs_Search = [7,6,5,4,3,2,1]
else
    if (Optim_iter == 1)
        global Configs_Search = [7,6,5,4,3,2,1]
    else
        global best_files   = readdir(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//1_plots//") 
        global last_feas_configs = filter(f -> startswith(f, string(Optim_iter-1)*"_Min_"), best_files)
        global last_feas_configs = [last_feas_configs[i][7:end-4] for i in 1:length(last_feas_configs)]
        global Configs_Search = [findfirst(configs_csv[2:end,1].==last_feas_configs[i]) for i in 1:length(last_feas_configs)]
    end
end

println("Looking into configuration: ", ID_configs[Configs_Search])

for config_i in Configs_Search

    # ****************** Load ibr config num ******************
    global ibr_conf_num = copy(config_i) ;
    load_ibrs(ibr_conf_num)              ;
    

    #********* Loop Initialization 2/2*********
    global mesh_loop                = 0     ; # Mesh counter
    global local_feasible_found     = false ; # Indicates if any feasible solution has been found for given configuration
    global mesh_feasible_found      = false ; # Indicates if any feasible solution has been found for given config. and mesh
    global continue_mesh_refinement = true  ;
    while (continue_mesh_refinement==true)
        global mesh_loop = mesh_loop + 1 ;

        # Solve without saturation for coarse mesh and if last mesh didnt found a feasible solution
        global include_saturation = true ; # <- Indicates if saturation is being included or not 

        println("Mesh_loop: ", mesh_loop ," | Phase: ", phase," | Configuration: ", id_config)
        #----------- Refine mesh ------------                            
        global (break_times_mnv,break_times) = custom_mesh(mesh_loop, phase , Optim_iter , ref_upd_time , H )
        
        #*************  Hermitte-Simpson collocation & precalculations ************ #
        global (t_ini,t_fin) = ( 0.0 , break_times_mnv[end] )                                                                                   ;
        global mid_points_mnv    = ((break_times_mnv[1:end-1]+break_times_mnv[2:end])/2)[(break_times_mnv[2:end]-break_times_mnv[1:end-1]).!=0] ;
        global t_p_sec_mnv       = sort(cat(break_times_mnv, mid_points_mnv, dims=1))                                                           ;
        global K_mnv             = size(t_p_sec_mnv,1)                  ;
        global K_upd_mnv         = zeros(K_mnv)                         ;
        for k in 1:(K_mnv-1)
            if   (t_p_sec_mnv[k] <= ref_upd_time) && (t_p_sec_mnv[k+1] <= ref_upd_time) global K_upd_mnv[k]   = 1                ; 
            else                                                                        global K_upd_mnv[k]   = 2                ; end
            if k==(K_mnv-1)                                                             global K_upd_mnv[end] = K_upd_mnv[end-1] ; end 
        end
        global K_upd_mnv = Int.(K_upd_mnv)
        global K_upd     = repeat( K_upd_mnv , Optim_iter)

        global t_p_sec       = repeat(t_p_sec_mnv,Optim_iter) ;
        global theta_DQ      = w_b*w_0*t_p_sec                ;
        global K = size(t_p_sec,1)
        
        #***************************************************************************************
        #************************************  Create model ************************************
        #***************************************************************************************

        println("Writing Optim. Program...")
        global model_traj      = Model(KNITRO.Optimizer) ;
        for optim_iter in 1:Optim_iter 
            # ----------------- Load exog. inputs (from max.) -----------------
            if optim_iter == 1
                global Pload_traj  =  repeat(M_Pload_fin[:,optim_iter]',K_mnv,1);
                global Qload_traj  =  repeat(M_Qload_fin[:,optim_iter]',K_mnv,1);
            else
                global Pload_traj  =  cat(Pload_traj,repeat(M_Pload_fin[:,optim_iter]',K_mnv,1),dims=1);
                global Qload_traj  =  cat(Qload_traj,repeat(M_Qload_fin[:,optim_iter]',K_mnv,1),dims=1);
            end
            #----------------- Reference variables (from max.) ----------------
            if optim_iter == 1
                global p_ref        =               cat(repeat(M_p_ref_ed[:,optim_iter]'     ,                  sum(1 .-(K_upd_mnv.==2))),  
                                                        repeat(M_p_ref_agc[:,optim_iter]'    , K_mnv           -sum(1 .-(K_upd_mnv.==2))), dims=1)           ;
                global q_ref        =                   repeat(M_q_ref_ed[:,optim_iter]'     , K_mnv          )                                              ;
            else
                global p_ref        =  cat(p_ref,   cat(repeat(M_p_ref_ed[:,optim_iter]'     ,                  sum(1 .-(K_upd_mnv.==2))),  
                                                        repeat(M_p_ref_agc[:,optim_iter]'    , K_mnv           -sum(1 .-(K_upd_mnv.==2))), dims=1) , dims=1) ;
                global q_ref        =  cat(q_ref,       repeat(M_q_ref_ed[:,optim_iter]'     , K_mnv          )                                    , dims=1) ;
            end
        end
        
        # -------------- IBR Parameters (Tunable) -------------
        global n1_ratio     =  @variable(model_traj, n1_ratio[n_ibr=1:Nibr] );
        global n2_ratio     =  @variable(model_traj, n2_ratio[n_ibr=1:Nibr] );
        # ----------------- Build model ------------------------
        global model_traj = model_generator(model_traj, true) 

        # ****************** Constraints ***********************

        # -------------- IBR Parameters (Tunable) -------------
        @constraint(model_traj , n1_ratio   .<=  n1_ratio_max ) ;
        @constraint(model_traj , n1_ratio   .>=  n1_ratio_min ) ;
        @constraint(model_traj , n2_ratio   .<=  n2_ratio_max ) ;
        @constraint(model_traj , n2_ratio   .>=  n2_ratio_min ) ;

        # ******** INITIAL STATE interlink with ED ********
        for optim_iter in 1:Optim_iter 
            for n_x in 1:Nx
                if typeof(Vars[1,n_x])==VariableRef
                    if (n_x>=hpq_ini) && (n_x<=hpq_end)
                        fix.(Vars[1+K_mnv*(optim_iter-1),n_x], M_Vars_ini[optim_iter,n_x-hpq_ini+ix_ini:n_x-hpq_end+ix_end] ; force=true) # hpq_0 = ix_0 if GFL IBR
                    else
                        fix.(Vars[1+K_mnv*(optim_iter-1),n_x], M_Vars_ini[optim_iter,n_x] ; force=true)
                    end
                end
            end
        end

        # **************** Collocation ********************
        for n_x in 1:Nx
            if typeof(Vars[1,n_x])==VariableRef
                for k in 1:(K-1)
                    if (in.(t_p_sec[k],Ref(break_times))==false)
                        global hk = t_p_sec[k+1] - t_p_sec[k-1];
                        @constraint(model_traj, Vars[k+1,n_x]   - Vars[k-1,n_x] - (hk/6)*Derivs[k-1,n_x] - (4*hk/6)*Derivs[k,n_x]  - (hk/6)*Derivs[k+1,n_x] ==  0)
                        @constraint(model_traj, Vars[k,n_x] - 0.5*Vars[k-1,n_x] - 0.5*Vars[k+1,n_x]      - (hk/8)*Derivs[k-1,n_x]  + (hk/8)*Derivs[k+1,n_x] ==  0)
                    elseif (t_p_sec[k] == t_p_sec[k+1])
                        @constraint(model_traj, Vars[k+1,n_x] - Vars[k,n_x] == 0  );
                    end
                end
            end
        end

        # ----------------- Objective ------------------------
        global object_value = @variable(model_traj, object_value       ) # <- dummy variable for obj >= obj_s for all s in S
        @constraint(model_traj,[optim_iter=1:Optim_iter],  object_value >= score[K_mnv*optim_iter])
        @constraint(model_traj, object_value >= 0.0                    )
        @constraint(model_traj, object_value <= 1.0                    )
        println("Calling optimiser... N° of variables: ",sum(typeof.(Vars).==VariableRef))

        # ---------------- Optimizer call --------------------
        global mesh_feasible_found       = false ; # Indicates if any feasible solution has been found for given initial point and mesh
        global round_count               = 0     ; # Round count: [1,2]. 1: Solve w/o adding noise to warmstart point, 2: Solve adding noise to warmstart point
        global msh_loop_obj              = []    ; # Collects objective f° values found
        global msh_loop_stop_criteria    = false ; # Mesh-Stopping criteria (based in round count [max round count = 2])
        while (msh_loop_stop_criteria==false)
            global round_count = round_count + 1 ;

            #-----------------------------------------------------
            #------------------ Warmstarting ---------------------
            #-----------------------------------------------------
            if (local_feasible_found==false)
                for optim_iter in 1:Optim_iter
                    if optim_iter==1
                        global ws_Vars      = repeat(M_Vars_ini[optim_iter:optim_iter,:],K_mnv,1)   ;
                    else
                        global ws_Vars      = cat(ws_Vars,repeat(M_Vars_ini[optim_iter:optim_iter,:],K_mnv,1),dims=1)   ;
                    end
                end
                global ws_Vars[:,end-Nibr:end-1] = ws_Vars[:,end-Nibr:end-1] .+ theta_DQ ;
            else
                global ws_Vars = zeros(K,Nx);
                for optim_iter in 1:Optim_iter 
                    for n_x in 1:Nx 
                        global ws_Vars[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,n_x] = linear_interpolation(
                                                    Interpolations.deduplicate_knots!(t_p_itp[optim_iter]), Vars_itp[optim_iter][:,n_x],
                                                                                                                extrapolation_bc=Linear())(t_p_sec[1+K_mnv*(optim_iter-1):K_mnv*optim_iter]);  
                    end
                end
            end

            # --------------------- Add noise (if round count = 2) ----------------------
            if round_count>1
                println("Adding noise to starting point")
                global noise_level = 0.1 ;
                global ws_Vars[:,end-Nibr:end-1] = ws_Vars[:,end-Nibr:end-1] .- theta_DQ                    ;
                global ws_Vars[:,:] = ws_Vars[:,:] .+ noise_level .* randn(size(ws_Vars[:,:])) .* maximum(abs.(ws_Vars[:,:]),dims=1)    ;
                global ws_Vars[:,end-Nibr:end-1] = ws_Vars[:,end-Nibr:end-1] .+ theta_DQ                    ;
                global ws_n1_ratio  = ws_n1_ratio  .+ noise_level .* randn(size(ws_n1_ratio)) .* ws_n1_ratio ;
                global ws_n2_ratio  = ws_n2_ratio  .+ noise_level .* randn(size(ws_n2_ratio)) .* ws_n2_ratio ;
            end

            # ********** Set warmstart **********
            set_start_value.( n1_ratio  , ws_n1_ratio  )        ;
            set_start_value.( n2_ratio  , ws_n2_ratio  )        ;
            for n_x in 1:Nx if typeof(Vars[1,n_x])==VariableRef set_start_value.(Vars[:,n_x], copy(ws_Vars[:,n_x])); end end
            set_start_value.( object_value , maximum([ws_Vars[K_mnv*optim_iter,end] for optim_iter in 1:Optim_iter]) )

            
            ######################################## Optimization config. ##########################################
            set_attribute(model_traj, "out_csvinfo"             , 1                         )
            set_attribute(model_traj, "bar_murule"              , 5                         )
            set_attribute(model_traj, "hessopt"                 , 6                         )
            set_attribute(model_traj, "bar_linsys_storage"      , 1                         )
            set_attribute(model_traj, "outlev"                  , 2                         )
            set_attribute(model_traj, "convex"                  , 0                         )
            set_attribute(model_traj, "strat_warm_start"        , 1                         )
            set_attribute(model_traj, "hessian_no_f"            , 1                         )
            set_attribute(model_traj, "linsolver_ordering"      , 0                         )
            set_attribute(model_traj, "algorithm"               , 1                         )
            set_attribute(model_traj, "scale"                   , 3                         )
            set_attribute(model_traj, "linsolver"               , 7                         )
            set_attribute(model_traj, "bar_directinterval"      , 50                        )
            set_attribute(model_traj, "presolve"                , 1                         )
            set_attribute(model_traj, "presolve_level"          , 2                         )
            set_attribute(model_traj, "presolve_initpt"         , 0                         )
            set_attribute(model_traj, "presolveop_redundant"    , 2                         )
            set_attribute(model_traj, "presolveop_substitution" , 2                         )
            set_attribute(model_traj, "presolveop_tighten"      , 2                         )
            set_attribute(model_traj, "ftol"                    , 1e-2                      )
            set_attribute(model_traj, "ftol_iters"              , 10                        )
            if mesh_loop==1
                set_attribute(model_traj, "maxit"               , 500                       )
                set_attribute(model_traj, "bar_initmu"          , 1e2                       )
                set_attribute(model_traj, "xtol"                , 1e-4                      )
                set_attribute(model_traj, "xtol_iters"          , 10                        )
            else
                set_attribute(model_traj, "maxit"               , 250                       )
                set_attribute(model_traj, "bar_initmu"          , [1e-8,1e-2][round_count]  )
                set_attribute(model_traj, "xtol"                , 1e-4                      )
                set_attribute(model_traj, "xtol_iters"          , [10,3][mesh_loop]         )
            end
            println("Looking for optimal solution, ",round_count,"° round")

            # ---------------------- Run optimizer -----------------------
            @objective(model_traj, Min, object_value) ;

            global num_csv_files  = length(readdir(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//0_Optim_csvs"))
            global save_filename  = BASE_SAVE_FOLDER*"//1_Mesh_Refinement//0_Optim_csvs//"*string(Optim_iter)*"_"*phase*"_"*string(config_i)*"_"*string(mesh_loop)*"_"*string(round_count)*"_min_"*string(num_csv_files)*".csv"
            set_attribute(model_traj, "out_csvname", save_filename)
            sleep(1)
            optimize!(model_traj)
            sleep(1)
            
            # ---------- Feasibility stage (look for closest feasibility point) ------------
            if (primal_status(model_traj)!=MOI.FEASIBLE_POINT)
                println("Non feasible solution... starting feasibility stage")
                global n1_ratio_feas    =   copy(value.(n1_ratio))  ;
                global n2_ratio_feas    =   copy(value.(n2_ratio))  ;
                global Vars_feas        =   copy(value.(Vars))      ;
                set_attribute(model_traj, "maxit"       , 50        )
                set_attribute(model_traj, "xtol"        , 1e-12     )
                set_attribute(model_traj, "xtol_iters"  , 3         )
                set_attribute(model_traj, "bar_initmu"  , 1e-12     )
                set_start_value.(n1_ratio  , n1_ratio_feas  )       ;
                set_start_value.(n2_ratio  , n2_ratio_feas  )       ;
                for n_x in 1:Nx 
                    if typeof(Vars[1,n_x])==VariableRef 
                        set_start_value.(Vars[:,n_x], copy(Vars_feas[:,n_x])); 
                    end 
                end
                set_start_value.( object_value , maximum([Vars_feas[K_mnv*optim_iter,end] for optim_iter in 1:Optim_iter]) )
                @objective(model_traj   ,  Min   , 0.0 )            ;

                global num_csv_files  = length(readdir(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//0_Optim_csvs"))
                global save_filename  = BASE_SAVE_FOLDER*"//1_Mesh_Refinement//0_Optim_csvs//"*string(Optim_iter)*"_"*phase*"_"*string(config_i)*"_"*string(mesh_loop)*"_"*string(round_count)*"_feas_"*string(num_csv_files)*".csv"
                set_attribute(model_traj, "out_csvname", save_filename)    
                sleep(1)
                optimize!(model_traj)
                sleep(1)
            end

            # ---------- Compare & save solution ------------
            if primal_status(model_traj)==MOI.FEASIBLE_POINT
                global JL_obj = value.(maximum([Vars[K_mnv*optim_iter,end] for optim_iter in 1:Optim_iter]))[1];
                if (mesh_feasible_found == false) global msh_loop_obj = [JL_obj]                   ;
                else                              global msh_loop_obj = push!(msh_loop_obj,JL_obj) ;   end

                global any_feasible_found       = true  ; # Indicates if any feasible solution has been found in entire max. phase
                global local_feasible_found     = true  ; # Indicates if any feasible solution has been found for given initial point
                global mesh_feasible_found      = true  ; # Indicates if any feasible solution has been found for given initial point and mesh
                
                if copy(JL_obj)<=minimum(msh_loop_obj)
                    println("Best value found")
                    # --------------- Save best solution ------------------
                    global Vars_res          = copy(value.(Vars))         ;
                    global Derivs_res        = copy(value.(Derivs))       ;
                    global ibr_conf_num_res  = copy(value.(ibr_conf_num)) ;
                    global y_theta_res       = copy(value.(y_theta))      ;
                    global n1_ratio_res      = copy(value.(n1_ratio))     ;
                    global n2_ratio_res      = copy(value.(n2_ratio))     ;
                    global JL_obj_res        = copy(value.(JL_obj))       ;
                    # --------------- Save for plotting -------------------
                    global p_ref_res       = copy(value.(p_ref))                          ;
                    global q_ref_res       = copy(value.(q_ref))                          ;
                    global Pload_traj_res  = copy(value.(Pload_traj))                     ;
                    global Qload_traj_res  = copy(value.(Qload_traj))                     ;
                    global Pibr_res        = copy(value.(p_ibr))                          ;
                    global Qibr_res        = copy(value.(q_ibr))                          ;
                    global w_dq_res        = copy(value.(w_dq))                           ;
                    global ixmag_res       = copy(value.(sqrt.(ix_d.^2 .+ ix_q.^2)))      ;
                    global vmag_res        = copy(value.(sqrt.(vibr_d.^2 .+ vibr_q.^2)))  ;
                    global score_res       = copy(value.(score))                          ;
                    #**** Save solution for interpolation in next mesh ****
                    for optim_iter in 1:Optim_iter 
                        if (optim_iter == 1)
                            global t_p_itp   = copy([     t_p_sec[1+K_mnv*(optim_iter-1):K_mnv*optim_iter]  ]  )  ;
                            global Vars_itp  = copy([    Vars_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:]]  )  ;
                        else
                            global t_p_itp   = push!(t_p_itp  ,      t_p_sec[1+K_mnv*(optim_iter-1):K_mnv*optim_iter]    )  ;
                            global Vars_itp  = push!(Vars_itp ,     Vars_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:]  )  ;
                        end
                    end

                end
                println("("*phase*". model) Objective value: "    , JL_obj ) ;  
            else
                println("("*phase*". model) Unfeasible.")
            end

            # ******** Mesh-Stopping criteria ********
            if (mesh_feasible_found==true)
                global msh_loop_stop_criteria = true
            elseif (round_count>=2)
                global msh_loop_stop_criteria = true
            end
        end
        
        # ******** Config-Stopping criteria ********
        if     (mesh_feasible_found == false)
            println("No feasible solution found...")
        elseif (mesh_loop==2) 
            include("Julia_1_3_EMT.jl")    ;
        end

        if (mesh_loop>=2)
            global continue_mesh_refinement = false;
            println("Done")
            println("eta_1: ", value.(n1_ratio)," - eta_2: ", value.(n2_ratio))
        else 
            println("Continuing mesh refinement")
        end
        println(repeat("-",100))
    end
end
