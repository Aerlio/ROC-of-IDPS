Random.seed!(1234*Optim_iter)

global phase            = "Max"                     ;
global H                = copy(H_max)               ; # Max phase time horizon
global ref_upd_time     = copy(ref_upd_time_max)    ; # Max phase AGC.upd.time
global separate_max     = true                      ; # true : Decomposes MINLP into a NLP based on load_id and jump_sense fixed integers
global add_ED_to_model  = true                      ; # true : Incorporates ED+AGC as constriants to the problem 

# --------- Loading previous solution ---------
if (Optim_iter == 1)
    global ( n1_ratio , n2_ratio , ibr_conf_num) = ( 0.01*ones(Nibr) , 0.01*ones(Nibr) , 7) ; # Initial guess: ALL_GFMs, n1 = 0.01, and n2 = 0.01
else
    global ( _ , _ , _ , _ , _ , _ , _ , _ , _ , _ , _ , ibr_conf_num , y_theta , n1_ratio , n2_ratio , _ ) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*string(Optim_iter-1)*"_Min_mnv_1");
end
load_ibrs(ibr_conf_num)

# ------------  Clean saving folder before starting (ask to do manually) --------------- 
println("Current number of files in saving folder: ", length(readdir(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//0_Optim_csvs")))

#********* Loop Initialization 1/2*********
global global_best_score        = 0     ; # Best score accross all starting points
global stop_criteria            = false ; # Criteria to stop looking for initial points (based on patience)
global patience_for_end         = 10    ; # Patience to stop looking for initial points
global patience_slv             = 0     ; # Patience counter
global ms_count                 = 0     ; # Counts n° of initial points tried

global any_feasible_found       = false ; # Indicates if any feasible solution has been found in entire max. phase
global local_feasible_found     = false ; # Indicates if any feasible solution has been found for given initial point
global mesh_feasible_found      = false ; # Indicates if any feasible solution has been found for given initial point and mesh

while (stop_criteria==false)
    global ms_count  =  ms_count + 1 ;

    # **************************** Looking for positive or negative PQ jump *********************************
    if (ms_count+Optim_iter-1)%2 == 1   global (load_id,jump_sense) = (1,"up")
    else                                global (load_id,jump_sense) = (1,"dn") end 

    # ****************** Load_warmstarting points (positive or negative PQ jump dependant) ******************
    global filename = BASE_SAVE_FOLDER*"//0_Warmstart_points//WS_Points_for_load_id_"*string(load_id)*"_and_jump_sense_"*jump_sense*".mat" ;
    global results = matopen(filename)
    global variables_of_interest = ["WS_p_ref_ed","WS_p_ref_agc","WS_q_ref_ed","WS_Pload_ini","WS_Pload_fin","WS_Qload_ini","WS_Qload_fin","WS_Vars_ed"]
    global (WS_p_ref_ed,WS_p_ref_agc,WS_q_ref_ed,WS_Pload_ini,WS_Pload_fin,WS_Qload_ini,WS_Qload_fin,WS_Vars_ed) = [read(results,i) for i in variables_of_interest]
    close(results)

    # ***** Only start counting patience if a feasible solution has been found for given initial point *******
    if (mesh_feasible_found==true)
        global patience_slv = patience_slv + 1 ;
    end

    # *************************** Multi-search stopping criteria *********************************
    if patience_slv < patience_for_end
        println(repeat("*",100))
        println(repeat("*",100))
        println(repeat("*",44)," ", ms_count, " Big loop ",repeat("*",44))
        println(repeat("*",40)," Patience ", patience_slv, " out of ", patience_for_end ,repeat("*",40))
        println(repeat("*",100))
        println(repeat("*",100)) 
        println(load_id, " - ", jump_sense)
        global continue_mesh_refinement = true  ;
    else
        global continue_mesh_refinement = false ;
        global stop_criteria            = true  ;
    end

    #********* Loop Initialization 2/2*********
    global mesh_loop            = 0     ; # Mesh counter
    global local_feasible_found = false ; # Indicates if any feasible solution has been found for given initial point
    while (continue_mesh_refinement==true)
        global mesh_loop = mesh_loop + 1 ;
        println("Max Mesh loop: ",mesh_loop)            
        #----------- Refine mesh ------------                                               
        global break_times        = custom_mesh(mesh_loop, phase , Optim_iter , ref_upd_time , H )
        global include_saturation = true ; # <- Indicates if saturation is being included or not 

        #*************  Hermitte-Simpson collocation & precalculations ************ #
        global (t_ini,t_fin) = ( 0.0 , break_times[end] )                                                                    ;
        global mid_points    = ((break_times[1:end-1]+break_times[2:end])/2)[(break_times[2:end]-break_times[1:end-1]).!=0]  ;
        global t_p_sec       = sort(cat(break_times, mid_points, dims=1))                                                    ;
        global theta_DQ      = w_b*w_0*t_p_sec                                                                               ;
        global K             = size(t_p_sec,1)                                                                               ;
        global K_upd         = zeros(K)                                                                                      ;
        for k in 1:(K-1)
            if   (t_p_sec[k] <= ref_upd_time) && (t_p_sec[k+1] <= ref_upd_time) global K_upd[k]   = 1            ; 
            else                                                                global K_upd[k]   = 2            ; end
            if k==(K-1)                                                         global K_upd[end] = K_upd[end-1] ; end 
        end
        global K_upd = Int.(K_upd)

        #***************************************************************************************
        #************************************  Create model ************************************
        #***************************************************************************************

        println("Writing Optim. Program...")
        global model_traj      = Model(KNITRO.Optimizer) ;
        # ----------------- Load exog. inputs -----------------
        global Pload_ini = @variable(  model_traj, Pload_ini[n_load=1:Nload] );
        global Qload_ini = @variable(  model_traj, Qload_ini[n_load=1:Nload] );        
        global (Pload_fin,Qload_fin)=(Vector{VariableRef}[],Vector{VariableRef}[])
        for n_load in 1:Nload
            if (n_load == load_id)
                global Pload_fin = cat( Pload_fin ,   @variable( model_traj , [1:1] ) ,dims=1);
                global Qload_fin = cat( Qload_fin ,   @variable( model_traj , [1:1] ) ,dims=1);
            else
                global Pload_fin = cat( Pload_fin , @expression( model_traj , Pload_ini[n_load] ) ,dims=1);
                global Qload_fin = cat( Qload_fin , @expression( model_traj , Qload_ini[n_load] ) ,dims=1);
            end
        end
        global Pload_traj =  @expression(model_traj, Pload_traj[k=1:K,n_load=1:Nload], Pload_fin[n_load] )  ;
        global Qload_traj =  @expression(model_traj, Qload_traj[k=1:K,n_load=1:Nload], Qload_fin[n_load] )  ;

        #----------- Reference variables (Tunable) ------------
        global p_ref_ed     =  @variable(model_traj,   p_ref_ed[n_ibr=1:Nibr]                                                                   )  ;
        global p_ref_agc    =  @variable(model_traj,   p_ref_agc[n_ibr=1:Nibr]                                                                  )  ;
        global q_ref_ed     =  @variable(model_traj,   q_ref_ed[n_ibr=1:Nibr]                                                                   )  ;
        global p_ref        =  @expression(model_traj, p_ref[k=1:K,n_ibr=1:Nibr], [ p_ref_ed[n_ibr] ; p_ref_agc[n_ibr]][Int(K_upd[k]==2)+1]     )  ;
        global q_ref        =  @expression(model_traj, q_ref[k=1:K,n_ibr=1:Nibr],   q_ref_ed[n_ibr]                                             )  ;

        # ----------------- Build model ------------------------
        global model_traj = model_generator(model_traj, true) 

        # ********************************************* Constraints ****************************************

        @constraint(model_traj,[n_load=1:Nload], Pload_ini[n_load] <= Pload_max[n_load]); # (ini)
        @constraint(model_traj,[n_load=1:Nload], Pload_ini[n_load] >= Pload_min[n_load]);
        @constraint(model_traj,[n_load=1:Nload], Qload_ini[n_load] <= Qload_max[n_load]);
        @constraint(model_traj,[n_load=1:Nload], Qload_ini[n_load] >= Qload_min[n_load]);
        if jump_sense == "up"
            @constraint(model_traj, Pload_fin[load_id] - Pload_ini[load_id] <= Pdif_max[load_id]  );
            @constraint(model_traj, Qload_fin[load_id] - Qload_ini[load_id] <= Qdif_max[load_id]  );
            @constraint(model_traj, Pload_fin[load_id] - Pload_ini[load_id] >= 0                  );
            @constraint(model_traj, Qload_fin[load_id] - Qload_ini[load_id] >= 0                  );
            @constraint(model_traj, Pload_fin[load_id]                      <= Pload_max[load_id] );
            @constraint(model_traj, Pload_fin[load_id]                      >= Pload_min[load_id] );
            @constraint(model_traj, Qload_fin[load_id]                      <= Qload_max[load_id] );
            @constraint(model_traj, Qload_fin[load_id]                      >= Qload_min[load_id] );
        else
            @constraint(model_traj, Pload_fin[load_id] - Pload_ini[load_id] <= 0                  );
            @constraint(model_traj, Qload_fin[load_id] - Qload_ini[load_id] <= 0                  );
            @constraint(model_traj, Pload_fin[load_id] - Pload_ini[load_id] >= -Pdif_max[load_id] );
            @constraint(model_traj, Qload_fin[load_id] - Qload_ini[load_id] >= -Qdif_max[load_id] );
            @constraint(model_traj, Pload_fin[load_id]                      <= Pload_max[load_id] );
            @constraint(model_traj, Pload_fin[load_id]                      >= Pload_min[load_id] );
            @constraint(model_traj, Qload_fin[load_id]                      <= Qload_max[load_id] );
            @constraint(model_traj, Qload_fin[load_id]                      >= Qload_min[load_id] );
        end

        # ******** INITIAL STATE interlink with ED ********
        include("Julia_0_1_ED.jl")      ;

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
        global object_value = @expression(model_traj, object_value, score[end])
        @constraint(model_traj, object_value >= 0.0                      );
        @constraint(model_traj, object_value <= 100.0                    );
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

            if (local_feasible_found == false) # Load warm-start from files generated from Julia_0_2_Initial_warmstarts.jl (positive or negative PQ jump dependant)
                global lenght_WS    = length(WS_p_ref_ed)
                global ws_p_ref_ed  = WS_p_ref_ed[rand(1:1:lenght_WS)]  ;
                global ws_q_ref_ed  = WS_q_ref_ed[rand(1:1:lenght_WS)]  ;
                global ws_p_ref_agc = WS_p_ref_agc[rand(1:1:lenght_WS)] ;
                global ws_Pload_ini = WS_Pload_ini[rand(1:1:lenght_WS)] ;
                global ws_Qload_ini = WS_Qload_ini[rand(1:1:lenght_WS)] ;
                global ws_Pload_fin = WS_Pload_fin[rand(1:1:lenght_WS)] ;
                global ws_Qload_fin = WS_Qload_fin[rand(1:1:lenght_WS)] ;
                global ws_Vars      = repeat(WS_Vars_ed[rand(1:1:lenght_WS)],K,1)        ;
                global ws_Vars[:,end-Nibr:end-1] = ws_Vars[:,end-Nibr:end-1] .+ theta_DQ ;
            else                # Load warm-start from interpolating solution from previous mesh_loop
                global ws_p_ref_ed   = copy(p_ref_ed_res)     ;
                global ws_q_ref_ed   = copy(q_ref_ed_res)     ;
                global ws_p_ref_agc  = copy(p_ref_agc_res)    ;
                global ws_Pload_ini  = copy(Pload_ini_res)    ;
                global ws_Pload_fin  = copy(Pload_fin_res)    ;
                global ws_Qload_ini  = copy(Qload_ini_res)    ;
                global ws_Qload_fin  = copy(Qload_fin_res)    ;
                global ws_Vars       = zeros(K,Nx)            ;
                for n_x in 1:Nx
                    global ws_Vars[:,n_x] = linear_interpolation(Interpolations.deduplicate_knots!(t_p_res), Vars_res[:,n_x], extrapolation_bc=Linear())(t_p_sec);  
                end
            end
            
            # --------------------- Add noise (if round count = 2) ----------------------
            if (round_count>1)
                println("Adding noise to starting point")
                global noise_level = 0.1 ;
                global ws_Vars[:,end-Nibr:end-1] = ws_Vars[:,end-Nibr:end-1] .- theta_DQ                            ;
                global ws_Vars = ws_Vars .+ noise_level .* randn(size(ws_Vars)) .* maximum(abs.(ws_Vars),dims=1)    ;
                global ws_Vars[:,end-Nibr:end-1] = ws_Vars[:,end-Nibr:end-1] .+ theta_DQ                            ;
                global ws_p_ref_ed  = ws_p_ref_ed  .+ noise_level .* randn(size(ws_p_ref_ed))  .* ws_p_ref_ed   ;
                global ws_q_ref_ed  = ws_q_ref_ed  .+ noise_level .* randn(size(ws_q_ref_ed))  .* ws_q_ref_ed   ;
                global ws_p_ref_agc = ws_p_ref_agc .+ noise_level .* randn(size(ws_p_ref_agc)) .* ws_p_ref_agc  ;
                global ws_Pload_ini = ws_Pload_ini .+ noise_level .* randn(size(ws_Pload_ini)) .* ws_Pload_ini  ;
                global ws_Pload_fin = ws_Pload_fin .+ noise_level .* randn(size(ws_Pload_fin)) .* ws_Pload_fin  ;
                global ws_Qload_ini = ws_Qload_ini .+ noise_level .* randn(size(ws_Qload_ini)) .* ws_Qload_ini  ;
                global ws_Qload_fin = ws_Qload_fin .+ noise_level .* randn(size(ws_Qload_fin)) .* ws_Qload_fin  ;
            end

            # ********** Set warmstart **********
            set_start_value.(p_ref_ed  , ws_p_ref_ed  ) ;
            set_start_value.(q_ref_ed  , ws_q_ref_ed  ) ;
            set_start_value.(p_ref_agc , ws_p_ref_agc ) ;
            set_start_value.(Pload_ini , ws_Pload_ini ) ;
            set_start_value.(Qload_ini , ws_Qload_ini ) ;
            set_start_value.(Pload_fin[load_id] , ws_Pload_fin[load_id] ) ;
            set_start_value.(Qload_fin[load_id] , ws_Qload_fin[load_id] ) ;
            for n_x in 1:Nx if typeof(Vars[1,n_x])==VariableRef set_start_value.(Vars[:,n_x], copy(ws_Vars[:,n_x])); end end

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
            @objective(model_traj, Max, object_value ) ;

            global num_csv_files  = length(readdir(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//0_Optim_csvs"))
            global save_filename  = BASE_SAVE_FOLDER*"//1_Mesh_Refinement//0_Optim_csvs//"*string(Optim_iter)*"_"*phase*"_"*string(ms_count)*"_"*string(mesh_loop)*"_"*string(round_count)*"_max_"*string(num_csv_files)*".csv"
            set_attribute(model_traj, "out_csvname", save_filename)
            sleep(0.1)
            optimize!(model_traj)
            sleep(0.1)

            # ---------- Feasibility stage (look for closest feasibility point) ------------
            if (primal_status(model_traj)!=MOI.FEASIBLE_POINT)
                println("Non-feasible... starting feasibility stage")
                global p_ref_ed_feas    =   copy(value.(p_ref_ed))  ;
                global q_ref_ed_feas    =   copy(value.(q_ref_ed))  ;
                global p_ref_agc_feas   =   copy(value.(p_ref_agc)) ;
                global Pload_ini_feas   =   copy(value.(Pload_ini)) ;
                global Pload_fin_feas   =   copy(value.(Pload_fin)) ;
                global Qload_ini_feas   =   copy(value.(Qload_ini)) ;
                global Qload_fin_feas   =   copy(value.(Qload_fin)) ;
                global Vars_feas        =   copy(value.(Vars))      ;
                set_attribute(model_traj, "maxit"       , 50      )
                set_attribute(model_traj, "bar_initmu"  , 1e-12   )
                set_attribute(model_traj, "xtol"        , 1e-12   )
                set_attribute(model_traj, "xtol_iters"  , 3       )
                set_start_value.(p_ref_ed  , p_ref_ed_feas  )       ;
                set_start_value.(q_ref_ed  , q_ref_ed_feas  )       ;
                set_start_value.(p_ref_agc , p_ref_agc_feas )       ;
                set_start_value.(Pload_ini , Pload_ini_feas )       ;
                set_start_value.(Qload_ini , Qload_ini_feas )       ;
                set_start_value.(Pload_fin[load_id] , Pload_fin_feas[load_id] ) ;
                set_start_value.(Qload_fin[load_id] , Qload_fin_feas[load_id] ) ;
                for n_x in 1:Nx if typeof(Vars[1,n_x])==VariableRef set_start_value.(Vars[:,n_x], copy(Vars_feas[:,n_x])); end end
                @objective(model_traj   ,  Max   , 0.0 )            ;
                
                global num_csv_files  = length(readdir(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//0_Optim_csvs"))
                global save_filename  = BASE_SAVE_FOLDER*"//1_Mesh_Refinement//0_Optim_csvs//"*string(Optim_iter)*"_"*phase*"_"*string(ms_count)*"_"*string(mesh_loop)*"_"*string(round_count)*"_feas_"*string(num_csv_files)*".csv"
                set_attribute(model_traj, "out_csvname", save_filename)
                sleep(0.1)
                optimize!(model_traj)
                sleep(0.1)
            end

            # ---------- Compare & save solution ------------
            if primal_status(model_traj)==MOI.FEASIBLE_POINT
                global JL_obj = value.(object_value)
                if (mesh_feasible_found == false) global msh_loop_obj = [JL_obj]                   ;
                else                              global msh_loop_obj = push!(msh_loop_obj,JL_obj) ;   end

                global any_feasible_found       = true  ; # Indicates if any feasible solution has been found in entire max. phase
                global local_feasible_found     = true  ; # Indicates if any feasible solution has been found for given initial point
                global mesh_feasible_found      = true  ; # Indicates if any feasible solution has been found for given initial point and mesh
                if copy(JL_obj)>=maximum(msh_loop_obj)
                    println("Best value found")
                    # -------------- Save best solution ------------------
                    global Vars_res          = copy(value.(Vars))         ;
                    global Derivs_res        = copy(value.(Derivs))       ;
                    global p_ref_ed_res      = copy(value.(p_ref_ed))     ;
                    global q_ref_ed_res      = copy(value.(q_ref_ed))     ;
                    global p_ref_agc_res     = copy(value.(p_ref_agc))    ;
                    global Pload_ini_res     = copy(value.(Pload_ini))    ;
                    global Qload_ini_res     = copy(value.(Qload_ini))    ;
                    global Pload_fin_res     = copy(value.(Pload_fin))    ;
                    global Qload_fin_res     = copy(value.(Qload_fin))    ;
                    global JL_obj_res        = copy(JL_obj)               ;
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
                    #**** Save mesh for interpolation in next mesh ****
                    global t_p_res   = copy(t_p_sec)  ;
                end
                println("("*phase*". model) Objective value: "    , JL_obj ) ;  
            else
                println("("*phase*". model) Unfeasible.")
            end

            # ******** Mesh-Stopping criteria ********
            if (mesh_loop==1) && (mesh_feasible_found==true)
                global msh_loop_stop_criteria = true
            elseif (round_count>=2)
                global msh_loop_stop_criteria = true
            end
        end

        # ******** WS-Stopping criteria ********
        if (mesh_feasible_found == false)
            println("No feasible solution found...")
        elseif (mesh_loop==2)
            include("Julia_1_3_EMT.jl")  ;
        end
        
        if (mesh_loop>=2)
            global continue_mesh_refinement = false;
            println("Done")
        else 
            println("Continuing mesh refinement")
        end
        println(repeat("-",100))
    end
end
