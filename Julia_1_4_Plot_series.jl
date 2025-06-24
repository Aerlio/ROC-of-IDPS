println("Plotting...")
global colors_by_gen   = palette(:seaborn_bright)[1:3]'     ;
global colors_by_load  = palette(:seaborn_bright)[7:end-1]' ;

# -------------------------------------------------------
# -------------------- Plotting  ------------------------
# -------------------------------------------------------

if (phase=="Min")
    global start_ML = (findall(==(-H/10), Time_ML[:,1]))
    for optim_iter in 1:Optim_iter
        if optim_iter==1
            if optim_iter==Optim_iter global end_ML=[length(Time_ML[:,1])]  
            else                      global end_ML=[start_ML[optim_iter+1]-1]  end
        else
            if optim_iter==Optim_iter global end_ML=push!(end_ML,length(Time_ML[:,1])        ) 
            else                      global end_ML=push!(end_ML,start_ML[optim_iter+1]-1  ) end
        end
    end
end

global l = []
for optim_iter in 1:Optim_iter
    if phase=="Max"
        # ************ JULIA **************
        global p_ref_JL = copy(p_ref_res)         ; 
        global q_ref_JL = copy(q_ref_res)         ; 
        global Pload_JL = copy(Pload_traj_res)    ; 
        global Qload_JL = copy(Qload_traj_res)    ; 
        global Pibr_JL  = copy(Pibr_res)          ;
        global Qibr_JL  = copy(Qibr_res)          ;
        global w_dq_JL  = copy(w_dq_res)          ;
        global ixmag_JL = copy(ixmag_res)         ;
        global vmag_JL  = copy(vmag_res)          ;
        global score_JL = copy(score_res)         ;

        # ************ MATLAB **************
        global p_ref_ML = Vars_ML[:,1:3]          ; 
        global q_ref_ML = Vars_ML[:,4:6]          ; 
        global Pload_ML = Vars_ML[:,7:9]          ; 
        global Qload_ML = Vars_ML[:,10:12]        ; 
        global Pibr_ML  = Vars_ML[:,13:15]        ; 
        global Qibr_ML  = Vars_ML[:,16:18]        ; 
        global w_dq_ML  = Vars_ML[:,19:21]        ; 
        global ixmag_ML = Vars_ML[:,22:24]        ; 
        global vmag_ML  = Vars_ML[:,25:27]        ; 
        global score_ML = Vars_ML[:,end]          ; 
    else
        # ************ JULIA **************
        global p_ref_JL = copy(p_ref_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ; 
        global q_ref_JL = copy(q_ref_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ; 
        global Pload_JL = copy(Pload_traj_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])    ; 
        global Qload_JL = copy(Qload_traj_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])    ; 
        global Pibr_JL  = copy(Pibr_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])          ;
        global Qibr_JL  = copy(Qibr_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])          ;
        global w_dq_JL  = copy(w_dq_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])          ;
        global ixmag_JL = copy(ixmag_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ;
        global vmag_JL  = copy(vmag_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])          ;
        global score_JL = copy(score_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ;

        # ************ MATLAB **************
        global p_ref_ML = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , 1:3]    ; 
        global q_ref_ML = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , 4:6]    ; 
        global Pload_ML = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , 7:9]    ; 
        global Qload_ML = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , 10:12]  ; 
        global Pibr_ML  = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , 13:15]  ; 
        global Qibr_ML  = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , 16:18]  ; 
        global w_dq_ML  = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , 19:21]  ; 
        global ixmag_ML = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , 22:24]  ; 
        global vmag_ML  = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , 25:27]  ; 
        global score_ML = Vars_ML[start_ML[optim_iter]:end_ML[optim_iter] , end]    ;
    end

    if (optim_iter==1) || (phase=="Max")
        global title_off_or_on = 2 # 1 is off, two is on
    else
        global title_off_or_on = 1 # 1 is off, two is on
    end
    
    if (optim_iter==Optim_iter) || (phase=="Max")
        global label_off_or_on = 2 # 1 is off, two is on
    else
        global label_off_or_on = 2 # 1 is off, two is on
    end

    if phase=="Max"
        global t_plot_JL = copy(t_p_sec)
        global t_plot_ML = copy(Time_ML)
    else
        global t_plot_JL = t_p_sec[1+K_mnv*(optim_iter-1):K_mnv*optim_iter]
        global t_plot_ML = copy(Time_ML[start_ML[optim_iter]:end_ML[optim_iter] , :])
    end
    # ************** Plotting comparison **************
    l_PQload = plot(    t_plot_ML, Pload_ML , label=[nothing,permutedims(["Pload("*string(i)*")" for i in 1:Nload])][label_off_or_on]     , xticks=[nothing,([0,ref_upd_time,H],["↑\n","AGC","H"])][label_off_or_on] , title=["","Pload/Qload (Exog. Input)"][title_off_or_on]  , color=colors_by_load, linewidth=[2,1.5,1]')
    l_PQload = plot!(   t_plot_ML, Qload_ML , label=[nothing,permutedims(["Qload("*string(i)*")" for i in 1:Nload])][label_off_or_on]     , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=colors_by_load, linewidth=[2,1.5,1]' , linestyle=:dashdot, legend_columns=3);
    l_PQload = scatter!( t_plot_ML[1:1], Qload_ML[1:1,1],label=" ", alpha=0.0);

    l_Pibr =  plot(   t_plot_ML, p_ref_ML , label=[nothing,permutedims(["Pref("*string(i)*")" for i in 1:Nibr])][label_off_or_on]         , xticks=[nothing,([0,ref_upd_time,H],["↑\n","AGC","H"])][label_off_or_on] , title=["","IBRs Active power"][title_off_or_on]          , color=colors_by_gen , linewidth=1.5 , linestyle=:dash, legend_columns=3);
    l_Pibr = plot!(   t_plot_ML, Pibr_ML  , label=[nothing,permutedims(["Pibr("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on]   , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=colors_by_gen , linewidth=1.5 );
    l_Pibr = scatter!(t_plot_JL, Pibr_JL  , label=[nothing,permutedims(["Pibr("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on]                                                                                                                                                  , color=colors_by_gen , markerstrokecolor=colors_by_gen, marker=:circle, alpha=0.2, markersize=2)

    l_Qibr = plot(    t_plot_ML, q_ref_ML , label=[nothing,permutedims(["Qref("*string(i)*")" for i in 1:Nibr])][label_off_or_on]         , xticks=[nothing,([0,ref_upd_time,H],["↑\n","AGC","H"])][label_off_or_on],  title=["","IBRs Reactive power"][title_off_or_on]        , color=colors_by_gen , linewidth=1.5 , linestyle=:dash, legend_columns=3);
    l_Qibr = plot!(   t_plot_ML, Qibr_ML  , label=[nothing,permutedims(["Qibr("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on]   , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=colors_by_gen , linewidth=1.5 );
    l_Qibr = scatter!(t_plot_JL, Qibr_JL  , label=[nothing,permutedims(["Qibr("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on]                                                                                                                                                  , color=colors_by_gen , markerstrokecolor=colors_by_gen, marker=:circle, alpha=0.2, markersize=2)

    l_Ix   =   plot(  t_plot_ML, ixmag_ML , label=[nothing,permutedims(["|I|("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on]    , xticks=[nothing,([0,ref_upd_time,H],["↑\n","AGC","H"])][label_off_or_on],  title=["","IBRs Current Magnitude"][title_off_or_on]     , color=colors_by_gen, linewidth=1.5 )
    l_Ix   = scatter!(t_plot_JL, ixmag_JL , label=[nothing,permutedims(["|I|("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on]    , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=colors_by_gen, markerstrokecolor=colors_by_gen, marker=:circle, alpha=0.2, markersize=2)
    l_Ix   = hline!(  ix_ibr_max[:,1]'  , label=[nothing,permutedims(["|I|max("*string(i)*")" for i in 1:Nibr])][label_off_or_on]                                                                                                                                                      , color=colors_by_gen, linewidth=1.2 , linestyle=:dash, legend_columns=3)

    l_V =  plot(  t_plot_ML , vmag_ML , label=[nothing,permutedims(["|V|("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on]     , xticks=[nothing,([0,ref_upd_time,H],["↑\n","AGC","H"])][label_off_or_on], title=["","IBRs Voltage Magnitude"][title_off_or_on]      , color=colors_by_gen , linewidth=1.5 )
    l_V = scatter!(t_plot_JL, vmag_JL , label=[nothing,permutedims(["|V|("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on]                                                                                                                                                    , color=colors_by_gen , markerstrokecolor=colors_by_gen, marker=:circle, alpha=0.2, markersize=2)
    l_V = hline!( [1]                 , label=[nothing,"V_nom"][label_off_or_on]                                                       , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=:black        , linewidth=1.2 , linestyle=:dash, alpha=0.2,legend_columns=3)

    l_freq = plot(   t_plot_ML , 50 .*w_dq_ML   , label=[nothing,permutedims(["f("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on], xticks=[nothing,([0,ref_upd_time,H],["↑\n","AGC","H"])][label_off_or_on], title=["","IBRs Frequency"][title_off_or_on], color=colors_by_gen , linewidth=1.5 );
    l_freq = scatter!(t_plot_JL, 50 .*w_dq_JL   , label=[nothing,permutedims(["f("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on], legend=[nothing,:outertop][label_off_or_on]                                                                             , color=colors_by_gen , markerstrokecolor=colors_by_gen, marker=:circle, legend_columns=3, alpha=0.2, markersize=2)
    l_freq = hline!( [50.5]    , label="f_max"                                                                                            , legend=[nothing,:outertop][label_off_or_on]                                                                            , color=:black, linewidth=1.2 , linestyle=:dash, alpha=0.2)
    l_freq = hline!( [49.5]    , label="f_min"                                                                                                                                                                                                                      , color=:black, linewidth=1.2 , linestyle=:dash, alpha=0.2)

    l_score =  plot(   t_plot_ML, score_ML       , label=[nothing,"Score [EMT]"][label_off_or_on]                                          , xticks=[nothing,([0,ref_upd_time,H],["↑\n","AGC","H"])][label_off_or_on], title=["","Score"][title_off_or_on]   , color=:red , linewidth =1.5  )
    l_score = scatter!(t_plot_JL, score_JL       , label=[nothing,"Score [Opt]"][label_off_or_on]                                          , legend=[nothing,:outertop][label_off_or_on]                                                                           , color=:red , markerstrokecolor=:red, marker=:circle, alpha=0.2, markersize=2)
    l_score = scatter!(t_plot_JL[1:1], score_JL[1:1]       , alpha= 0.0, label=" ")

    global l = push!(l,l_PQload,l_Pibr,l_Qibr,l_Ix,l_V,l_freq,l_score)
    if (phase=="Max")
        break
    end
end

if      (Optim_iter == 1) || (phase=="Max")
    global plot_emt = plot(   (li for li in l)..., 
                                layout=grid(1,7), size=(2500,300), legendfontsize=8,
                                xlabelfontsize=9,
                                top_margin = 4Plots.mm, bottom_margin = 4Plots.mm, left_margin=1Plots.mm, right_margin=1Plots.mm, xlim=(-H/10,H) )
else 
    global plot_emt = plot(   (li for li in l)..., 
                                layout=grid(Optim_iter,7), size=(2500,300*Optim_iter), legendfontsize=8,
                                xlabelfontsize=9,
                                top_margin = 4Plots.mm, bottom_margin = 4Plots.mm, left_margin=3Plots.mm, right_margin=3Plots.mm, xlim=(-H/10,H) )
end

if (phase=="Max")
    savefig(plot_emt, BASE_SAVE_FOLDER*"//1_Mesh_Refinement//1_plots//"*string(Optim_iter)*"_"*phase*"_"*string(ms_count)*".svg")
else
    savefig(plot_emt, BASE_SAVE_FOLDER*"//1_Mesh_Refinement//1_plots//"*string(Optim_iter)*"_"*phase*"_"*string(ID_configs[ibr_conf_num])*".svg")
end
