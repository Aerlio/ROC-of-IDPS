include("Julia_0_0_General_Calls.jl") ;
global add_ED_to_model   = false      ; # false: Solves a modified ED+AGC optim. problem 
global separate_max      = true       ; # true : Decomposes MINLP into a NLP based on load_id and jump_sense fixed integers
global Num_WS_points     = 200        ; # Number of WS points to generate per NLP

# ------------  Clean saving folder before starting --------------- 
global WS_SAVE_FOLDER = BASE_SAVE_FOLDER*"//0_Warmstart_points";
if length(readdir(WS_SAVE_FOLDER*"//Optim_csvs")) != 0 
    [rm(WS_SAVE_FOLDER*"//Optim_csvs//"*file) for file in readdir(WS_SAVE_FOLDER*"//Optim_csvs")]
end

for slv in 1:2
    if     (slv==1)
        global load_id      = 1      ;
        global jump_sense   = "up"   ;
    elseif (slv==2)
        global load_id      = 1      ;
        global jump_sense   = "dn"   ;
    end

    if (jump_sense == "up") global agc_sign =  1 ;
    else                    global agc_sign = -1 ; end
    for WS_points in 1:Num_WS_points
        println(repeat("-",100))
        println("N° warmstart point being generated: ", WS_points)
        global ws_ix     = -10 .+ 20 .*rand(Nibr)
        global ws_agc    = -10 .+ 20 .*rand(Nibr)
        global Pload_ini = [rand(Pload_min[n_load,1]:0.05:Pload_max[n_load,1]) for n_load in 1:Nload] ;
        global Qload_ini = [rand(Qload_min[n_load,1]:0.05:Qload_max[n_load,1]) for n_load in 1:Nload] ;
        global Pload_fin = copy(Pload_ini);
        global Qload_fin = copy(Qload_ini);
        global Pload_fin[load_id] = max(min(Pload_fin[load_id] + agc_sign*rand()*Pdif_max[load_id] , Pload_max[load_id]) , Pload_min[load_id]) ;
        global Qload_fin[load_id] = max(min(Qload_fin[load_id] + agc_sign*rand()*Qdif_max[load_id] , Qload_max[load_id]) , Qload_min[load_id]) ;

        sleep(0.5)
        include("Julia_0_1_ED.jl") ; # solve modified ED+AGC problem
        sleep(0.5)

        println("Jump sense: ",jump_sense, " | Load jumping: ", load_id)
        println("P(ED): " ,round.(value.(p_ref_ed)   , digits=3))
        println("P(AGC): ",round.(value.(p_ref_agc)  , digits=3))
        println("Q(ED): " ,round.(value.(q_ref_ed)   , digits=3))
        println("P(ini): ",round.(value.(Pload_ini)  , digits=3))
        println("P(fin): ",round.(value.(Pload_fin)  , digits=3))
        println("Q(ini): ",round.(value.(Qload_ini)  , digits=3))
        println("Q(fin): ",round.(value.(Qload_fin)  , digits=3))
        println(repeat("-",100))
        if WS_points==1
            global WS_p_ref_ed  = [ value.(p_ref_ed)  ]
            global WS_p_ref_agc = [ value.(p_ref_agc) ]
            global WS_q_ref_ed  = [ value.(q_ref_ed)  ]
            global WS_Pload_ini = [ value.(Pload_ini) ]
            global WS_Pload_fin = [ value.(Pload_fin) ]
            global WS_Qload_ini = [ value.(Qload_ini) ]
            global WS_Qload_fin = [ value.(Qload_fin) ]
            global WS_Vars_ed   = [ value.(Vars_ed)   ]
        else
            global WS_p_ref_ed  = push!( WS_p_ref_ed  , value.(p_ref_ed)  )
            global WS_p_ref_agc = push!( WS_p_ref_agc , value.(p_ref_agc) )
            global WS_q_ref_ed  = push!( WS_q_ref_ed  , value.(q_ref_ed)  )
            global WS_Pload_ini = push!( WS_Pload_ini , value.(Pload_ini) )
            global WS_Pload_fin = push!( WS_Pload_fin , value.(Pload_fin) )
            global WS_Qload_ini = push!( WS_Qload_ini , value.(Qload_ini) )
            global WS_Qload_fin = push!( WS_Qload_fin , value.(Qload_fin) )
            global WS_Vars_ed   = push!( WS_Vars_ed   , value.(Vars_ed)   )
        end
    end

    # ------------- Save Results for simulation -----------
    global filename = WS_SAVE_FOLDER*"//WS_Points_for_load_id_"*string(load_id)*"_and_jump_sense_"*jump_sense*".mat" ;
    global results = matopen(filename, "w")
    write(results, "WS_p_ref_ed"    , WS_p_ref_ed   );
    write(results, "WS_p_ref_agc"   , WS_p_ref_agc  );
    write(results, "WS_q_ref_ed"    , WS_q_ref_ed   );
    write(results, "WS_Pload_ini"   , WS_Pload_ini  );
    write(results, "WS_Pload_fin"   , WS_Pload_fin  );
    write(results, "WS_Qload_ini"   , WS_Qload_ini  );
    write(results, "WS_Qload_fin"   , WS_Qload_fin  );
    write(results, "WS_Vars_ed"     , WS_Vars_ed    );
    close(results)
end