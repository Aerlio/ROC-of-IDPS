println("Plotting...")

# -------------------------------------------------------
# -------------------- Plotting  ------------------------
# -------------------------------------------------------
using Colors

global colors_by_gen_opt = parse.(Colorant, ["blue2","red1","green3"])';  
global colors_by_gen_emt = parse.(Colorant, ["cyan","orange","greenyellow"])';
global colors_by_load_opt = parse.(Colorant, ["blueviolet","gray15","goldenrod1"])';  
global colors_by_load_emt = parse.(Colorant, ["magenta1","gray50","yellow"])';


if (phase=="Min")
    global start_EMT = (findall(==(-H/10), ts_EMT[:,1]))
    for optim_iter in 1:Optim_iter
        if optim_iter==1
            if optim_iter==Optim_iter global end_EMT=[length(ts_EMT[:,1])]  
            else                      global end_EMT=[start_EMT[optim_iter+1]-1]  end
        else
            if optim_iter==Optim_iter global end_EMT=push!(end_EMT,length(ts_EMT[:,1])        ) 
            else                      global end_EMT=push!(end_EMT,start_EMT[optim_iter+1]-1  ) end
        end
    end
end

global l = []
for optim_iter in 1:Optim_iter
    
    if phase=="Max"
        # ************ JULIA **************
        global p_ref_OPT = copy(p_ref_res)         ; 
        global q_ref_OPT = copy(q_ref_res)         ; 
        global Pdist_OPT = copy(p_dist_res)        ; 
        global Qdist_OPT = copy(q_dist_res)        ; 
        global Pibr_OPT  = copy(Pibr_res)          ;
        global Qibr_OPT  = copy(Qibr_res)          ;
        global Pload_OPT = copy(pload_res)         ; 
        global Qload_OPT = copy(qload_res)         ;
        global w_dq_OPT  = copy(w_dq_res)          ;
        global ixmag_OPT = copy(ixmag_res)         ;
        global vmag_OPT  = copy(vmag_res)          ;
        global score_OPT = copy(score_res)         ;

        # ************ MATLAB **************
        global p_ref_EMT = Vars_EMT[:,1:3]          ; 
        global q_ref_EMT = Vars_EMT[:,4:6]          ; 
        global Pload_EMT = Vars_EMT[:,7:9]          ; 
        global Qload_EMT = Vars_EMT[:,10:12]        ; 
        global Pibr_EMT  = Vars_EMT[:,13:15]        ; 
        global Qibr_EMT  = Vars_EMT[:,16:18]        ; 
        global w_dq_EMT  = Vars_EMT[:,19:21]        ; 
        global ixmag_EMT = Vars_EMT[:,22:24]        ; 
        global vmag_EMT  = Vars_EMT[:,25:27]        ; 
        global score_EMT = Vars_EMT[:,end]          ; 

        # To plot PLOAD_INI and QLOAD_INI for Max phase. This is to recover missing saved info in "for_plotting" files using "results" files instead
        global Opt_results = matopen(BASE_SAVE_FOLDER*"//Optimization//"*call_name*".mat")
        global (PLOAD_INI,QLOAD_INI) = [read(Opt_results,i) for i in ["PLOAD_INI","QLOAD_INI"]];
        close(Opt_results)
    else
        # ************ JULIA **************
        global p_ref_OPT = copy(p_ref_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ; 
        global q_ref_OPT = copy(q_ref_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ; 
        global Pdist_OPT = copy(p_dist_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])        ; 
        global Qdist_OPT = copy(q_dist_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])        ;
        global Pibr_OPT  = copy(Pibr_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])          ;
        global Qibr_OPT  = copy(Qibr_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])          ;
        global Pload_OPT = copy(pload_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ; 
        global Qload_OPT = copy(qload_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ;
        global w_dq_OPT  = copy(w_dq_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])          ;
        global ixmag_OPT = copy(ixmag_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ;
        global vmag_OPT  = copy(vmag_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])          ;
        global score_OPT = copy(score_res[1+K_mnv*(optim_iter-1):K_mnv*optim_iter,:])         ;

        # ************ MATLAB **************
        global p_ref_EMT = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , 1:3]    ; 
        global q_ref_EMT = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , 4:6]    ; 
        global Pload_EMT = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , 7:9]    ; 
        global Qload_EMT = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , 10:12]  ; 
        global Pibr_EMT  = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , 13:15]  ; 
        global Qibr_EMT  = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , 16:18]  ; 
        global w_dq_EMT  = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , 19:21]  ; 
        global ixmag_EMT = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , 22:24]  ; 
        global vmag_EMT  = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , 25:27]  ; 
        global score_EMT = Vars_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , end]    ;
        
        # To plot PLOAD_INI and QLOAD_INI for Max phase. This is to recover missing saved info in "for_plotting" files using "results" files instead
        global Opt_results = matopen(BASE_SAVE_FOLDER*"//Optimization//"*call_name[1:6]*"mnv_"*string(optim_iter)*"_"*call_name[7:end]*".mat")
        global (PLOAD_INI,QLOAD_INI) = [read(Opt_results,i) for i in ["PLOAD_INI","QLOAD_INI"]];
        close(Opt_results)
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
        global t_plot_OPT = copy(t_p_sec)
        global t_plot_EMT = copy(ts_EMT)
    else
        global t_plot_OPT = t_p_sec[1+K_mnv*(optim_iter-1):K_mnv*optim_iter]
        global t_plot_EMT = copy(ts_EMT[start_EMT[optim_iter]:end_EMT[optim_iter] , :])
    end
    # ************** Plotting comparison **************
    
    lw_EMT = 2
    alpha_EMT=0.7
    markersize_OPT=2.5
    alpha_OPT = 0.7;
    l_Pload =  plot(   cat(t_plot_OPT[1]-H/10,t_plot_OPT[1],t_plot_OPT,dims=1), cat(repeat(PLOAD_INI',2),Pdist_OPT,dims=1)  , label=[nothing,permutedims(["Pdist("*string(i)*")" for i in 1:Nload])][label_off_or_on]         , xticks=[nothing,([0,ref_upd_time,H],["","",""])][label_off_or_on] , title=["","Pload"][title_off_or_on] , color=colors_by_load_opt , linewidth=1.5 , linestyle=:solid, legend_columns=3);
    l_Pload = scatter!(t_plot_OPT, Pload_OPT  , label=[nothing,permutedims(["Pload("*string(i)*") [Opt]" for i in 1:Nload])][label_off_or_on]                                                                                                                                                  , color=colors_by_load_opt , markerstrokecolor=colors_by_load_opt, marker=:circle,  alpha=alpha_OPT,markersize=markersize_OPT)
    l_Pload =  plot!(   t_plot_EMT, Pload_EMT  , alpha=alpha_EMT, label=[nothing,permutedims(["Pload("*string(i)*") [EMT]" for i in 1:Nload])][label_off_or_on]   , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=colors_by_load_emt , linewidth=lw_EMT );
    
    l_Qload =  plot(   cat(t_plot_OPT[1]-H/10,t_plot_OPT[1],t_plot_OPT,dims=1), cat(repeat(QLOAD_INI',2),Qdist_OPT,dims=1)  , label=[nothing,permutedims(["Qdist("*string(i)*")" for i in 1:Nload])][label_off_or_on]         , xticks=[nothing,([0,ref_upd_time,H],["","",""])][label_off_or_on] , title=["","Qload"][title_off_or_on] , color=colors_by_load_opt , linewidth=1.5 , linestyle=:solid, legend_columns=3);
    l_Qload = scatter!(t_plot_OPT, Qload_OPT  , label=[nothing,permutedims(["Qload("*string(i)*") [Opt]" for i in 1:Nload])][label_off_or_on]                                                                                                                                                  , color=colors_by_load_opt , markerstrokecolor=colors_by_load_opt, marker=:circle,  alpha=alpha_OPT,markersize=markersize_OPT)
    l_Qload = plot!(   t_plot_EMT, Qload_EMT  , alpha=alpha_EMT, label=[nothing,permutedims(["Qload("*string(i)*") [EMT]" for i in 1:Nload])][label_off_or_on]   , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=colors_by_load_emt , linewidth=lw_EMT );
    
    l_Pibr =  plot(   cat(t_plot_OPT[1]-H/10,t_plot_OPT[1],t_plot_OPT,dims=1), cat(p_ref_OPT[1:2,:],p_ref_OPT,dims=1) , label=[nothing,permutedims(["Pref("*string(i)*")" for i in 1:Nibr])][label_off_or_on]         , xticks=[nothing,([0,ref_upd_time,H],["","",""])][label_off_or_on] , title=["","Active power"][title_off_or_on]          , color=colors_by_gen_opt , linewidth=1.5 , linestyle=:solid, legend_columns=3);
    l_Pibr = scatter!(t_plot_OPT, Pibr_OPT  , label=[nothing,permutedims(["Pibr("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on]                                                                                                                                                  , color=colors_by_gen_opt , markerstrokecolor=colors_by_gen_opt, marker=:circle,  alpha=alpha_OPT,markersize=markersize_OPT)
    l_Pibr = plot!(   t_plot_EMT, Pibr_EMT  , alpha=alpha_EMT, label=[nothing,permutedims(["Pibr("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on]   , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=colors_by_gen_emt , linewidth=lw_EMT );
    
    l_Qibr = plot(    cat(t_plot_OPT[1]-H/10,t_plot_OPT[1],t_plot_OPT,dims=1), cat(q_ref_OPT[1:2,:],q_ref_OPT,dims=1) , label=[nothing,permutedims(["Qref("*string(i)*")" for i in 1:Nibr])][label_off_or_on]         , xticks=[nothing,([0,ref_upd_time,H],["","",""])][label_off_or_on],  title=["","Reactive power"][title_off_or_on]        , color=colors_by_gen_opt , linewidth=1.5 , linestyle=:solid, legend_columns=3);
    l_Qibr = scatter!(t_plot_OPT, Qibr_OPT  , label=[nothing,permutedims(["Qibr("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on]                                                                                                                                                  , color=colors_by_gen_opt , markerstrokecolor=colors_by_gen_opt, marker=:circle,  alpha=alpha_OPT,markersize=markersize_OPT)
    l_Qibr = plot!(   t_plot_EMT, Qibr_EMT  , alpha=alpha_EMT, label=[nothing,permutedims(["Qibr("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on]   , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=colors_by_gen_emt , linewidth=lw_EMT );
    
    l_Ix   = scatter(t_plot_OPT, ixmag_OPT , label=[nothing,permutedims(["|I|("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on]    , legend=[nothing,:outertop][label_off_or_on]                                                                                               , color=colors_by_gen_opt, markerstrokecolor=colors_by_gen_opt, marker=:circle,  alpha=alpha_OPT,markersize=markersize_OPT)
    l_Ix   = hline!(  ix_ibr_max[:,1]'  , label=[nothing,permutedims(["|I|max("*string(i)*")" for i in 1:Nibr])][label_off_or_on]                                                                                                                                                      , color=colors_by_gen_opt, linewidth=1.2 , linestyle=:dash, legend_columns=3)
    l_Ix   =   plot!(  t_plot_EMT, ixmag_EMT , alpha=alpha_EMT , label=[nothing,permutedims(["|I|("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on]    , xticks=[nothing,([0,ref_upd_time,H],["","",""])][label_off_or_on],  title=["","Current Magnitude"][title_off_or_on]     , color=colors_by_gen_emt, linewidth=lw_EMT )
    
    l_V = scatter(t_plot_OPT, vmag_OPT , label=[nothing,permutedims(["|V|("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on]                                                                                                                                                    , color=colors_by_gen_opt , markerstrokecolor=colors_by_gen_opt, marker=:circle,  alpha=alpha_OPT,markersize=markersize_OPT)
    l_V = hline!( [1] , label=[nothing,"V_nom"][label_off_or_on] , legend=[nothing,:outertop][label_off_or_on] , color=:grey , linewidth=1.2 , linestyle=:dash, alpha=0.5,legend_columns=3)
    l_V =  plot!(  t_plot_EMT , vmag_EMT  , alpha=alpha_EMT, label=[nothing,permutedims(["|V|("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on]     , xticks=[nothing,([0,ref_upd_time,H],["","",""])][label_off_or_on], title=["","Voltage Magnitude"][title_off_or_on]      , color=colors_by_gen_emt , linewidth=lw_EMT )
    
    l_freq = hline( [50]    , label="f_nom" , legend=[nothing,:outertop][label_off_or_on] , color=:grey, linewidth=1.2 , linestyle=:dash, alpha=0.5)
    l_freq = scatter!(t_plot_OPT, 50 .*w_dq_OPT   , label=[nothing,permutedims(["f("*string(i)*") [Opt]" for i in 1:Nibr])][label_off_or_on], legend=[nothing,:outertop][label_off_or_on]                                                                             , color=colors_by_gen_opt , markerstrokecolor=colors_by_gen_opt, marker=:circle, legend_columns=3,  alpha=alpha_OPT,markersize=markersize_OPT)
    l_freq = plot!(   t_plot_EMT , 50 .*w_dq_EMT  , alpha=alpha_EMT  , label=[nothing,permutedims(["f("*string(i)*") [EMT]" for i in 1:Nibr])][label_off_or_on], xticks=[nothing,([0,ref_upd_time,H],["","",""])][label_off_or_on], title=["","Frequency"][title_off_or_on], color=colors_by_gen_emt , linewidth=lw_EMT );
    
    l_score =  plot(   t_plot_EMT, score_EMT      , alpha=alpha_EMT  , label=[nothing,"Score [EMT]"][label_off_or_on]                                          , xticks=[nothing,([0,ref_upd_time,H],["","",""])][label_off_or_on], title=["","Score"][title_off_or_on]   , color=:red , linewidth =1.5  )
    l_score = scatter!(t_plot_OPT, score_OPT       , label=[nothing,"Score [Opt]"][label_off_or_on]                                          , legend=[nothing,:outertop][label_off_or_on]                                                                           , color=:red , markerstrokecolor=:red, marker=:circle,  alpha=alpha_OPT,markersize=markersize_OPT)
    
    global l = push!(l,l_Pload,l_Qload,l_Pibr,l_Qibr,l_Ix,l_V,l_freq,l_score)
    if (phase=="Max")
        break
    end
end

if      (Optim_iter == 1) || (phase=="Max")
    global plot_emt = plot(   (li for li in l)..., 
                                layout=grid(1,8), size=(1500,150),title="", legend=false,tick_direction = :in,
                                top_margin = 0Plots.mm, bottom_margin = 0Plots.mm, left_margin=0Plots.mm, right_margin=0Plots.mm, xlim=(-H/10,H) )
else 
    global plot_emt = plot(   (li for li in l)..., 
                                layout=grid(Optim_iter,8), size=(1500,150*Optim_iter),title="", legend=false,tick_direction = :in,
                                top_margin = 0Plots.mm, bottom_margin = 0Plots.mm, left_margin=0Plots.mm, right_margin=0Plots.mm, xlim=(-H/10,H) )
end

savefig(plot_emt, BASE_SAVE_FOLDER*"//Plots//"*call_name*".svg")