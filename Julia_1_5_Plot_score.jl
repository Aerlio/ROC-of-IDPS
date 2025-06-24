global V_n1_ratio   = 0.01*ones(Nibr)   ;
global V_n2_ratio   = 0.01*ones(Nibr)   ;

for phase_read in ["Max","Min"]
    for optim_iter in 1:Optim_iter
        if phase_read=="Max"
            global (mnv_break_times   , mnv_t_p           , mnv_Vars          , mnv_Derivs        ,
                    mnv_Pload_ini     , mnv_Pload_fin     , mnv_Qload_ini     , mnv_Qload_fin     , 
                    mnv_p_ref_ed      , mnv_q_ref_ed      , mnv_p_ref_agc     , 
                    mnv_ibr_conf_num  , mnv_y_theta       , mnv_n1_ratio , mnv_n2_ratio , mnv_score) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*string(optim_iter)*"_Max");            
            
            if optim_iter==1 global V_max_scores = [mnv_score/H_max]
            else             global V_max_scores = push!(V_max_scores,mnv_score/H_max) ; end
        else
            
            global (mnv_break_times   , mnv_t_p           , mnv_Vars          , mnv_Derivs        ,
                    mnv_Pload_ini     , mnv_Pload_fin     , mnv_Qload_ini     , mnv_Qload_fin     , 
                    mnv_p_ref_ed      , mnv_q_ref_ed      , mnv_p_ref_agc     , 
                    mnv_ibr_conf_num  , mnv_y_theta       , mnv_n1_ratio , mnv_n2_ratio , mnv_score) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*string(optim_iter)*"_Min_mnv_"*string(optim_iter));            
            
            if optim_iter==1 global V_min_scores = [mnv_score/H_min]
            else             global V_min_scores = push!(V_min_scores,mnv_score/H_min) ; end

            global V_n1_ratio   = cat(V_n1_ratio,mnv_n1_ratio,dims=2)
            global V_n2_ratio   = cat(V_n2_ratio,mnv_n2_ratio,dims=2)

        end
    end
end

# ******************** Plotting ********************

global colors_by_gen   = palette(:seaborn_bright)[1:3]'     ;
global colors_by_load  = palette(:seaborn_bright)[7:end-1]' ;
global l = @layout([A{0.001h}; [B C D]])
global title = plot(title = "GFM/GFL Configuration: "*string(mnv_ibr_conf_num), grid = false, showaxis = false)
global plot_minmax =   plot(1:Optim_iter,log.(1 .+V_max_scores) , ylabel="log(obj. f°)", color=:red    , label="Max. block", title="Min-Max Objective f° Evolution", marker=:circle)
global plot_minmax =  plot!(1:Optim_iter,log.(1 .+V_min_scores) , ylabel="log(obj. f°)", color=:green  , label="Min. block", xlabel="Optim. Iteration", marker=:circle, legend_columns=2, xticks=(collect(1:Optim_iter),string.(Int.(collect(1:Optim_iter)))))
global plot_n1_params = plot(0:Optim_iter,V_n1_ratio' , color=colors_by_gen , ylabel="p.u.",xlabel="Optim. Iteration" , label=permutedims(["IBR(1)","IBR(2)","IBR(3)"]),title="IBRs n1_ratio Evolution", marker=:circle, legend_columns=3, xticks=(collect(0:Optim_iter),string.(Int.(collect(0:Optim_iter)))))
global plot_n2_params = plot(0:Optim_iter,V_n2_ratio' , color=colors_by_gen , ylabel="p.u.",xlabel="Optim. Iteration" , label=permutedims(["IBR(1)","IBR(2)","IBR(3)"]),title="IBRs n2_ratio Evolution", marker=:circle, legend_columns=3, xticks=(collect(0:Optim_iter),string.(Int.(collect(0:Optim_iter)))))

global plot_final=plot(title, plot_minmax,plot_n1_params,plot_n2_params, 
                        legend=:outertop ,layout=l,size=(1200,400), top_margin = 5Plots.mm, bottom_margin = 10Plots.mm, left_margin = 15Plots.mm)

savefig(plot_final,BASE_SAVE_FOLDER*"//1_Mesh_Refinement//1_plots//0_final.svg")
savefig(plot_final,BASE_SAVE_FOLDER*"//1_Mesh_Refinement//2_best_plots//0_final.svg")