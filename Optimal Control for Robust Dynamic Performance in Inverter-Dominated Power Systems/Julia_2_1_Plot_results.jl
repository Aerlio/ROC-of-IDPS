include("Julia_0_0_General_Calls.jl") ;

global selected_files = readdir(BASE_SAVE_FOLDER*"//Optimization//");

# ---------------------------------------------------------
# ----------------- Read MAX phase results ----------------
# ---------------------------------------------------------

for optim_iter in 1:Max_Scenarios
    if optim_iter == 1
        global NLP_time_max     = zeros( Max_Scenarios , NLP_max );
        global EMT_time_max     = zeros( Max_Scenarios , NLP_max );
        global NLP_max_scores   = zeros( Max_Scenarios , NLP_max );
        global EMT_max_scores   = zeros( Max_Scenarios , NLP_max );
    end
    global max_files_optim_iter = filter(f -> startswith(f, string(optim_iter)*"_Max_"), selected_files)
    global max_files_optim_iter = filter(f -> endswith(f, "_for_sim.mat"), max_files_optim_iter)
    global max_filenames_optim_iter = [max_file[1:end-12] for max_file in max_files_optim_iter]
    global max_filenames_optim_iter = sort(max_filenames_optim_iter, by=length)

    for i in 1:NLP_max
        global (_ ,_,_,_,_,_,_,_,_,_,_,_,_,_,_, OBJ_NLP, TIME_NLP, OBJ_EMT, TIME_EMT ) = read_results(BASE_SAVE_FOLDER*"//Optimization//"*max_filenames_optim_iter[i]);
        global NLP_time_max[optim_iter,i] = copy(TIME_NLP)
        global EMT_time_max[optim_iter,i] = copy(TIME_EMT)
        global NLP_max_scores[optim_iter,i] = copy(OBJ_NLP)
        global EMT_max_scores[optim_iter,i] = copy(OBJ_EMT)
    end
end

global worst_EMT_scores =  maximum(EMT_max_scores,dims=2)[:,1]

# ---------------------------------------------------------
# ---------------- Read MIN phase results -------------
# ---------------------------------------------------------

for optim_iter in 2:Max_Scenarios
    if optim_iter == 2
        global Conf_nums = [length(NLP_min)]
        global N1_ratio  = zeros(Max_Scenarios,Nibr)
        global N2_ratio  = zeros(Max_Scenarios,Nibr)
        global N1_ratio[1,:] = 0.01*ones(Nibr)
        global N2_ratio[1,:] = 0.01*ones(Nibr)
        global best_EMT_scores                  = []
        global all_EMT_min_scores               = []
        global all_NLP_min_scores               = []
        
        global EMT_min_score_All_GFMs           = []
        global EMT_min_score_All_GFLs_but_Gen1  = []
        global EMT_min_score_All_GFLs_but_Gen2  = []
        global EMT_min_score_All_GFLs_but_Gen3  = []
        global EMT_min_score_All_GFMs_but_Gen1  = []
        global EMT_min_score_All_GFMs_but_Gen2  = []
        global EMT_min_score_All_GFMs_but_Gen3  = []

        global NLP_min_score_All_GFMs           = []
        global NLP_min_score_All_GFLs_but_Gen1  = []
        global NLP_min_score_All_GFLs_but_Gen2  = []
        global NLP_min_score_All_GFLs_but_Gen3  = []
        global NLP_min_score_All_GFMs_but_Gen1  = []
        global NLP_min_score_All_GFMs_but_Gen2  = []
        global NLP_min_score_All_GFMs_but_Gen3  = []
        
        global EMT_min_time_All_GFMs            = []
        global EMT_min_time_All_GFLs_but_Gen1   = []
        global EMT_min_time_All_GFLs_but_Gen2   = []
        global EMT_min_time_All_GFLs_but_Gen3   = []
        global EMT_min_time_All_GFMs_but_Gen1   = []
        global EMT_min_time_All_GFMs_but_Gen2   = []
        global EMT_min_time_All_GFMs_but_Gen3   = []

        global NLP_min_time_All_GFMs            = []
        global NLP_min_time_All_GFLs_but_Gen1   = []
        global NLP_min_time_All_GFLs_but_Gen2   = []
        global NLP_min_time_All_GFLs_but_Gen3   = []
        global NLP_min_time_All_GFMs_but_Gen1   = []
        global NLP_min_time_All_GFMs_but_Gen2   = []
        global NLP_min_time_All_GFMs_but_Gen3   = []
    end
    global min_files_optim_iter = filter(f -> startswith(f, string(optim_iter-1)*"_Min_"), selected_files)
    global min_files_optim_iter = filter(f -> endswith(f, "_for_sim.mat"), min_files_optim_iter)
    global min_filenames_optim_iter = [min_file[1:end-12] for min_file in min_files_optim_iter]
    global min_filenames_optim_iter = sort(min_filenames_optim_iter, by=length)

    println(min_filenames_optim_iter)
    global All_GFMs_file = filter(s -> occursin("All_GFMs-", s), min_filenames_optim_iter.*"-")
    global All_GFLs_but_Gen1_file = filter(s -> occursin("All_GFLs_but_Gen1-", s), min_filenames_optim_iter.*"-")
    global All_GFLs_but_Gen2_file = filter(s -> occursin("All_GFLs_but_Gen2-", s), min_filenames_optim_iter.*"-")
    global All_GFLs_but_Gen3_file = filter(s -> occursin("All_GFLs_but_Gen3-", s), min_filenames_optim_iter.*"-")
    global All_GFMs_but_Gen1_file = filter(s -> occursin("All_GFMs_but_Gen1-", s), min_filenames_optim_iter.*"-")
    global All_GFMs_but_Gen2_file = filter(s -> occursin("All_GFMs_but_Gen2-", s), min_filenames_optim_iter.*"-")
    global All_GFMs_but_Gen3_file = filter(s -> occursin("All_GFMs_but_Gen3-", s), min_filenames_optim_iter.*"-")
    
    global (_,_,_,_,_,_,_,_,_,_,_,ibr_conf_num,_,n1_ratio,n2_ratio, OBJ_NLP, _ , OBJ_EMT, _) = read_results(BASE_SAVE_FOLDER*"//Optimization//Selected_Results//"*string(optim_iter-1)*"_Min_mnv_1");

    global Conf_nums = push!(Conf_nums,ibr_conf_num)
    global N1_ratio[optim_iter,:] = copy(n1_ratio)
    global N2_ratio[optim_iter,:] = copy(n2_ratio)
    global best_EMT_scores = push!(best_EMT_scores,OBJ_EMT)

    if size(All_GFMs_file,1)>0           global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,OBJ_NLP_7, TIME_NLP_7, OBJ_EMT_7, TIME_EMT_7 ) = read_results(BASE_SAVE_FOLDER*"//Optimization//"*All_GFMs_file[1][1:end-1]);          else global OBJ_EMT_7 = 4321 ;  end 
    if size(All_GFLs_but_Gen1_file,1)>0  global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,OBJ_NLP_6, TIME_NLP_6, OBJ_EMT_6, TIME_EMT_6 ) = read_results(BASE_SAVE_FOLDER*"//Optimization//"*All_GFLs_but_Gen1_file[1][1:end-1]); else global OBJ_EMT_6 = 4321 ;  end 
    if size(All_GFLs_but_Gen2_file,1)>0  global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,OBJ_NLP_5, TIME_NLP_5, OBJ_EMT_5, TIME_EMT_5 ) = read_results(BASE_SAVE_FOLDER*"//Optimization//"*All_GFLs_but_Gen2_file[1][1:end-1]); else global OBJ_EMT_5 = 4321 ;  end 
    if size(All_GFLs_but_Gen3_file,1)>0  global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,OBJ_NLP_4, TIME_NLP_4, OBJ_EMT_4, TIME_EMT_4 ) = read_results(BASE_SAVE_FOLDER*"//Optimization//"*All_GFLs_but_Gen3_file[1][1:end-1]); else global OBJ_EMT_4 = 4321 ;  end 
    if size(All_GFMs_but_Gen1_file,1)>0  global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,OBJ_NLP_3, TIME_NLP_3, OBJ_EMT_3, TIME_EMT_3 ) = read_results(BASE_SAVE_FOLDER*"//Optimization//"*All_GFMs_but_Gen1_file[1][1:end-1]); else global OBJ_EMT_3 = 4321 ;  end 
    if size(All_GFMs_but_Gen2_file,1)>0  global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,OBJ_NLP_2, TIME_NLP_2, OBJ_EMT_2, TIME_EMT_2 ) = read_results(BASE_SAVE_FOLDER*"//Optimization//"*All_GFMs_but_Gen2_file[1][1:end-1]); else global OBJ_EMT_2 = 4321 ;  end 
    if size(All_GFMs_but_Gen3_file,1)>0  global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,OBJ_NLP_1, TIME_NLP_1, OBJ_EMT_1, TIME_EMT_1 ) = read_results(BASE_SAVE_FOLDER*"//Optimization//"*All_GFMs_but_Gen3_file[1][1:end-1]); else global OBJ_EMT_1 = 4321 ;  end 
    
    if OBJ_EMT_7!=4321 
        global EMT_min_score_All_GFMs          = push!( EMT_min_score_All_GFMs          , OBJ_EMT_7)
        global NLP_min_score_All_GFMs          = push!( NLP_min_score_All_GFMs          , OBJ_NLP_7)
        global EMT_min_time_All_GFMs           = push!( EMT_min_time_All_GFMs           , TIME_EMT_7)
        global NLP_min_time_All_GFMs           = push!( NLP_min_time_All_GFMs           , TIME_NLP_7)
        global (all_EMT_min_scores, all_NLP_min_scores) = (push!(all_EMT_min_scores,OBJ_EMT_7), push!(all_NLP_min_scores,OBJ_NLP_7))
    end
    if OBJ_EMT_6!=4321 
        global EMT_min_score_All_GFLs_but_Gen1 = push!( EMT_min_score_All_GFLs_but_Gen1 , OBJ_EMT_6)
        global NLP_min_score_All_GFLs_but_Gen1 = push!( NLP_min_score_All_GFLs_but_Gen1 , OBJ_NLP_6)
        global EMT_min_time_All_GFLs_but_Gen1  = push!( EMT_min_time_All_GFLs_but_Gen1  , TIME_EMT_6)
        global NLP_min_time_All_GFLs_but_Gen1  = push!( NLP_min_time_All_GFLs_but_Gen1  , TIME_NLP_6)
        global (all_EMT_min_scores, all_NLP_min_scores) = (push!(all_EMT_min_scores,OBJ_EMT_6), push!(all_NLP_min_scores,OBJ_NLP_6))
    end       
    if OBJ_EMT_5!=4321 
        global EMT_min_score_All_GFLs_but_Gen2 = push!( EMT_min_score_All_GFLs_but_Gen2 , OBJ_EMT_5)
        global NLP_min_score_All_GFLs_but_Gen2 = push!( NLP_min_score_All_GFLs_but_Gen2 , OBJ_NLP_5)
        global EMT_min_time_All_GFLs_but_Gen2  = push!( EMT_min_time_All_GFLs_but_Gen2  , TIME_EMT_5)
        global NLP_min_time_All_GFLs_but_Gen2  = push!( NLP_min_time_All_GFLs_but_Gen2  , TIME_NLP_5)
        global (all_EMT_min_scores, all_NLP_min_scores) = (push!(all_EMT_min_scores,OBJ_EMT_5), push!(all_NLP_min_scores,OBJ_NLP_5))
    end       
    if OBJ_EMT_4!=4321 
        global EMT_min_score_All_GFLs_but_Gen3 = push!( EMT_min_score_All_GFLs_but_Gen3 , OBJ_EMT_4)
        global NLP_min_score_All_GFLs_but_Gen3 = push!( NLP_min_score_All_GFLs_but_Gen3 , OBJ_NLP_4)
        global EMT_min_time_All_GFLs_but_Gen3  = push!( EMT_min_time_All_GFLs_but_Gen3  , TIME_EMT_4)
        global NLP_min_time_All_GFLs_but_Gen3  = push!( NLP_min_time_All_GFLs_but_Gen3  , TIME_NLP_4)
        global (all_EMT_min_scores, all_NLP_min_scores) = (push!(all_EMT_min_scores,OBJ_EMT_4), push!(all_NLP_min_scores,OBJ_NLP_4))
    end       
    if OBJ_EMT_3!=4321 
        global EMT_min_score_All_GFMs_but_Gen1 = push!( EMT_min_score_All_GFMs_but_Gen1 , OBJ_EMT_3)
        global NLP_min_score_All_GFMs_but_Gen1 = push!( NLP_min_score_All_GFMs_but_Gen1 , OBJ_NLP_3)
        global EMT_min_time_All_GFMs_but_Gen1  = push!( EMT_min_time_All_GFMs_but_Gen1  , TIME_EMT_3)
        global NLP_min_time_All_GFMs_but_Gen1  = push!( NLP_min_time_All_GFMs_but_Gen1  , TIME_NLP_3)
        global (all_EMT_min_scores, all_NLP_min_scores) = (push!(all_EMT_min_scores,OBJ_EMT_3), push!(all_NLP_min_scores,OBJ_NLP_3))
    end      
    if OBJ_EMT_2!=4321 
        global EMT_min_score_All_GFMs_but_Gen2 = push!( EMT_min_score_All_GFMs_but_Gen2 , OBJ_EMT_2)
        global NLP_min_score_All_GFMs_but_Gen2 = push!( NLP_min_score_All_GFMs_but_Gen2 , OBJ_NLP_2)
        global EMT_min_time_All_GFMs_but_Gen2  = push!( EMT_min_time_All_GFMs_but_Gen2  , TIME_EMT_2)
        global NLP_min_time_All_GFMs_but_Gen2  = push!( NLP_min_time_All_GFMs_but_Gen2  , TIME_NLP_2)
        global (all_EMT_min_scores, all_NLP_min_scores) = (push!(all_EMT_min_scores,OBJ_EMT_2), push!(all_NLP_min_scores,OBJ_NLP_2))
    end     
    if OBJ_EMT_1!=4321 
        global EMT_min_score_All_GFMs_but_Gen3 = push!( EMT_min_score_All_GFMs_but_Gen3 , OBJ_EMT_1)
        global NLP_min_score_All_GFMs_but_Gen3 = push!( NLP_min_score_All_GFMs_but_Gen3 , OBJ_NLP_1)
        global EMT_min_time_All_GFMs_but_Gen3  = push!( EMT_min_time_All_GFMs_but_Gen3  , TIME_EMT_1)
        global NLP_min_time_All_GFMs_but_Gen3  = push!( NLP_min_time_All_GFMs_but_Gen3  , TIME_NLP_1)
        global (all_EMT_min_scores, all_NLP_min_scores) = (push!(all_EMT_min_scores,OBJ_EMT_1), push!(all_NLP_min_scores,OBJ_NLP_1))
    end     
end


config_strings=[    "GFM: [2,3] / GFL: [1]", 
                    "GFM: [1] / GFL: [2,3]",
                    "GFM: [1,3] / GFL: [2]", 
                    "GFM: [2] / GFL: [1,3]",
                    "GFM: [1,2] / GFL: [3]", 
                    "GFM: [3] / GFL: [1,2]", 
                    "GFM: [1,2,3]"]
config_symbols=[    :rect, :dtriangle, :diamond, :circle, :utriangle, :pentagon, :star5]
global colors_by_gen   = palette(:seaborn_bright)[1:3]' ;

# ---------------------------------------------------------
# --------------- Local Reduction Main Results ------------
# ---------------------------------------------------------

l = @layout([a b c])

x_min_max = cat(collect(Iterators.flatten(zip(1:Max_Scenarios, 1:(Max_Scenarios-1)))),Max_Scenarios, dims=1)
y_min_max = log1p.(cat(collect(Iterators.flatten(zip(worst_EMT_scores, best_EMT_scores./H_min))),worst_EMT_scores[end], dims=1))

plot_minmax=plot(1:Max_Scenarios     , log1p.(worst_EMT_scores)        , linewidth=4, color=:red    , label="Max-phase ROCP solution", yticks=(log1p.([0.1,0.5,1,2,3,6,12]), string.([0.1,0.5,1,2,3,6,12])) ) 
plot_minmax=plot!(1:(Max_Scenarios-1) , log1p.(best_EMT_scores./H_min)  , linewidth=4, color=:green  , label="Min-phase ROCP solution", legend=:topright, title="Local Reduction Results")
plot_minmax=plot!(x_min_max, y_min_max , linewidth=3 , linestyle=:dash, ylabel="Objective function, Ψ", xlabel="Local reduction iteration", marker=:x, markersize=4 ,  color=:navy, label="Local Reduction trajectory")


plot_max = plot(1:Max_Scenarios, log1p.(worst_EMT_scores) , linewidth=4 , ylabel="Objective function, Ψ", xlabel="Local reduction iteration", markersize=11 ,  color=:red  , label="Selected solution")
plot_max = scatter!(1:Max_Scenarios, log1p.(EMT_max_scores) , alpha=0.6   , markersize=11   ,  color=:red , marker=:circle , label="" , yticks=(log1p.([0.1,0.5,1,2,3,6,12]), string.([0.1,0.5,1,2,3,6,12])))
plot_max = scatter!(1:1, log1p.( EMT_max_scores[1:1]        ) , alpha=0.6   , markersize=11   ,  color=:red , marker=:circle, legend=:topright , label="Solutions across initializations", title="Max-phase ROCP")

x=1
plot_min = plot( 1:(Max_Scenarios-1)                       , log1p.( best_EMT_scores./H_min ), ylim=(0.08,0.4) , yticks=(log1p.([0.1,0.2,0.3,0.4]), string.([0.1,0.2,0.3,0.4])) , linewidth=4,  color=:green   , ylabel="Objective function, Ψ"  ,legend=:topleft, xlabel="Local reduction iteration", label="Selected solution", title="Min-phase ROCP")
plot_min = plot!(1:length(EMT_min_score_All_GFMs)          , log1p.( EMT_min_score_All_GFMs./H_min )          ,alpha=0.8, linewidth=1.4, color=x+1,  label="GFM [1,2,3]", marker=:star5, markersize=11)
plot_min = plot!(1:length(EMT_min_score_All_GFMs_but_Gen1) , log1p.( EMT_min_score_All_GFMs_but_Gen1./H_min ) ,alpha=0.8, linewidth=1.4, color=x+2,  label="GFM [2,3] / GFL [1]", marker=:rect, markersize=11)
plot_min = plot!(1:length(EMT_min_score_All_GFMs_but_Gen2) , log1p.( EMT_min_score_All_GFMs_but_Gen2./H_min ) ,alpha=0.8, linewidth=1.4, color=x+3,  label="GFM [1,3] / GFL [2]", marker=:diamond, markersize=11)
plot_min = plot!(1:length(EMT_min_score_All_GFMs_but_Gen3) , log1p.( EMT_min_score_All_GFMs_but_Gen3./H_min ) ,alpha=0.8, linewidth=1.4, color=x+4,  label="GFM [1,2] / GFL [3]", marker=:utriangle, markersize=11)
plot_min = plot!(1:length(EMT_min_score_All_GFLs_but_Gen1) , log1p.( EMT_min_score_All_GFLs_but_Gen1./H_min ) ,alpha=0.8, linewidth=1.4, color=x+5,  label="GFM [1] / GFL [2,3]", marker=:dtriangle, markersize=11)
plot_min = plot!(1:length(EMT_min_score_All_GFLs_but_Gen2) , log1p.( EMT_min_score_All_GFLs_but_Gen2./H_min ) ,alpha=0.8, linewidth=1.4, color=x+6,  label="GFM [2] / GFL [1,3]", marker=:circle, markersize=11)
plot_min = plot!(1:length(EMT_min_score_All_GFLs_but_Gen3) , log1p.( EMT_min_score_All_GFLs_but_Gen3./H_min ) ,alpha=0.8, linewidth=1.4, color=x+7,  label="GFM [3] / GFL [1,2]", marker=:pentagon, markersize=11, legend_columns=1)

plot_final = plot(plot_minmax,plot_max,plot_min, layout=l,size=(2000,500), bottom_margin = 15Plots.mm, top_margin = 5Plots.mm,  left_margin = 15Plots.mm ,
                       xtickfontsize=14,ytickfontsize=14 ,  #legend_background_color=:transparent,
                       xlabelfontsize=16,ylabelfontsize=16, titlefontsize=20,legendfontsize=14)

savefig(plot_final,BASE_SAVE_FOLDER*"//Plots//Selected_Results_Plots//0_Final_Dynamic_Perf_Results.svg")

# ---------------------------------------------------------
# -------------------- Control Parameters -----------------
# ---------------------------------------------------------


l = @layout([a b])
eta_1 = plot(0:(Max_Scenarios-1), N1_ratio , linewidth=2,color=colors_by_gen, label=permutedims(["IBR 1","IBR 2","IBR 3"]), ylim=(-0.001,0.025), title="K^P")
eta_1 = scatter!(0:(Max_Scenarios-1), N1_ratio , color=colors_by_gen, markershape=config_symbols[Conf_nums], markersize=14, alpha=1, label="", xlabel="Local reduction iteration", ylabel="Param. value")
eta_2 = plot(0:(Max_Scenarios-1), N2_ratio , linewidth=2,color=colors_by_gen, ylim=(-0.005,0.065),legend=false, title="K^Q")
eta_2 = scatter!(0:(Max_Scenarios-1), N2_ratio , color=colors_by_gen, markershape=config_symbols[Conf_nums], markersize=14, alpha=1, label="", xlabel="Local reduction iteration", ylabel="Param. value")

plot_eta=plot(eta_1,eta_2, layout=l,size=(1600,500), bottom_margin = 15Plots.mm, top_margin = 5Plots.mm,  left_margin = 15Plots.mm ,
                       xtickfontsize=14,ytickfontsize=14 , legend_background_color=:transparent,
                       xlabelfontsize=16,ylabelfontsize=16, titlefontsize=20,legendfontsize=14)

                       
savefig(plot_eta,BASE_SAVE_FOLDER*"//Plots//Selected_Results_Plots//0_Final_Control_Params.svg")


# ---------------------------------------------------------
# -------------------- Dynamic Assessment -----------------
# ---------------------------------------------------------

dyn_assess_max=vec(100*abs.(NLP_max_scores.-EMT_max_scores)./EMT_max_scores);
dyn_assess_min=vec(100*abs.(all_NLP_min_scores .- all_EMT_min_scores)./all_EMT_min_scores);
dyn_assess = cat(dyn_assess_max,dyn_assess_min,dims=1)
mean(dyn_assess), std(dyn_assess)

using StatsPlots
plot_dyn_assess = violin(["Max-phase ROCP"] , log1p.(dyn_assess_max); label=false, color=:red  , alpha=0.4, ylabel="Relative Error (%)", title="Dynamic assessment of mesh refinement method")
plot_dyn_assess = violin!(["Min-phase ROCP"], log1p.(dyn_assess_min); label=false, color=:green , alpha=0.4)
plot_dyn_assess = scatter!(0.5 .* ones(size(dyn_assess_max))  , log1p.(dyn_assess_max); label=false, color=:red  , alpha=0.8)
plot_dyn_assess = scatter!(1.9  .* ones(size(dyn_assess_min)) , log1p.(dyn_assess_min); label=false, color=:green, alpha=0.8)

dyn_assess_plot =plot(plot_dyn_assess, layout=l,size=(800,500), bottom_margin = 15Plots.mm, top_margin = 5Plots.mm,  left_margin = 15Plots.mm ,
                       xtickfontsize=14,ytickfontsize=14 , yticks=(log1p.([0,6,12,25,50,100]), string.([0,6,12,25,50,100])),
                       xlabelfontsize=16,ylabelfontsize=16, titlefontsize=20,legendfontsize=14)

savefig(dyn_assess_plot,BASE_SAVE_FOLDER*"//Plots//Selected_Results_Plots//0_Dynamic_Assessment.svg")

# ---------------------------------------------------------
# -------------------- Optimization time ------------------
# ---------------------------------------------------------

total_min_time=zeros((Max_Scenarios-1))
avg_min_time=zeros((Max_Scenarios-1))
for i in 1:(Max_Scenarios-1)
    global confs_solved = 0
    if size(NLP_min_time_All_GFMs,1)>=i            global confs_solved = confs_solved + 1 ; global total_min_time[i] = total_min_time[i] + NLP_min_time_All_GFMs[i] end
    if size(NLP_min_time_All_GFMs_but_Gen1,1)>=i   global confs_solved = confs_solved + 1 ; global total_min_time[i] = total_min_time[i] + NLP_min_time_All_GFMs_but_Gen1[i] end
    if size(NLP_min_time_All_GFMs_but_Gen2,1)>=i   global confs_solved = confs_solved + 1 ; global total_min_time[i] = total_min_time[i] + NLP_min_time_All_GFMs_but_Gen2[i] end
    if size(NLP_min_time_All_GFMs_but_Gen3,1)>=i   global confs_solved = confs_solved + 1 ; global total_min_time[i] = total_min_time[i] + NLP_min_time_All_GFMs_but_Gen3[i] end
    if size(NLP_min_time_All_GFLs_but_Gen1,1)>=i   global confs_solved = confs_solved + 1 ; global total_min_time[i] = total_min_time[i] + NLP_min_time_All_GFLs_but_Gen1[i] end
    if size(NLP_min_time_All_GFLs_but_Gen2,1)>=i   global confs_solved = confs_solved + 1 ; global total_min_time[i] = total_min_time[i] + NLP_min_time_All_GFLs_but_Gen2[i] end
    if size(NLP_min_time_All_GFLs_but_Gen3,1)>=i   global confs_solved = confs_solved + 1 ; global total_min_time[i] = total_min_time[i] + NLP_min_time_All_GFLs_but_Gen3[i] end
    global avg_min_time[i] = total_min_time[i]/confs_solved
end


l = @layout([a b])
plot_max_time = scatter(1:Max_Scenarios, NLP_time_max./60 , alpha=0.6   , markersize=11   ,  color=:red , label="" , marker=:circle , title="Max-phase ROCP Optimization time")
plot_max_time = scatter!(1:1, [NLP_time_max[1,1]./60 ], alpha=0.6   , markersize=11   ,  color=:red , label="Solutions across multiple initializations" , marker=:circle )
plot_max_time = plot!(1:Max_Scenarios, vec(sum(NLP_time_max,dims=2))./60 , linewidth=4,   color=:red , label="Total optimization time" , ylabel="Time (min)", xlabel="Local reduction iteration",ylim=(-2,30), legend=:topright)
plot_max_time = plot!(1:Max_Scenarios, vec(mean(NLP_time_max,dims=2))./60 , linewidth=2,   color=:red , linestyle=:dash, label="Average optimization time")

plot_min_time = plot(1:length(NLP_min_time_All_GFMs)          , NLP_min_time_All_GFMs./60          ,alpha=0.8, color=x+1, linewidth=1.4 ,  label="GFM [1,2,3]", marker=:star5, markersize=11)
plot_min_time = plot!(1:length(NLP_min_time_All_GFMs_but_Gen1) , NLP_min_time_All_GFMs_but_Gen1./60 ,alpha=0.8, color=x+2, linewidth=1.4,  label="GFM [2,3] / GFL [1]", marker=:rect, markersize=11)
plot_min_time = plot!(1:length(NLP_min_time_All_GFMs_but_Gen2) , NLP_min_time_All_GFMs_but_Gen2./60 ,alpha=0.8, color=x+3, linewidth=1.4,  label="GFM [1,3] / GFL [2]", marker=:diamond, markersize=11)
plot_min_time = plot!(1:length(NLP_min_time_All_GFMs_but_Gen3) , NLP_min_time_All_GFMs_but_Gen3./60 ,alpha=0.8, color=x+4, linewidth=1.4,  label="GFM [1,2] / GFL [3]", marker=:utriangle, markersize=11)
plot_min_time = plot!(1:length(NLP_min_time_All_GFLs_but_Gen1) , NLP_min_time_All_GFLs_but_Gen1./60 ,alpha=0.8, color=x+5, linewidth=1.4,  label="GFM [1] / GFL [2,3]", marker=:dtriangle, markersize=11)
plot_min_time = plot!(1:length(NLP_min_time_All_GFLs_but_Gen2) , NLP_min_time_All_GFLs_but_Gen2./60 ,alpha=0.8, color=x+6, linewidth=1.4,  label="GFM [2] / GFL [1,3]", marker=:circle, markersize=11)
plot_min_time = plot!(1:length(NLP_min_time_All_GFLs_but_Gen3) , NLP_min_time_All_GFLs_but_Gen3./60 ,alpha=0.8, color=x+7, linewidth=1.4,  label="GFM [3] / GFL [1,2]", marker=:pentagon, markersize=11, legend_columns=1)
plot_min_time = plot!(1:(Max_Scenarios-1) , total_min_time./60, linewidth=4, color=:green, label="Total optimization time", ylabel="Time (min)"  ,legend=:topright, xlabel="Local reduction iteration", title="Min-phase ROCP Optimization time", 
                    ylim=(-2,140))
plot_min_time = plot!(1:(Max_Scenarios-1) , avg_min_time./60, linewidth=2, color=:green, linestyle=:dash, legend_columns=2, label="Average optimization time")

plot_times=plot(plot_max_time,plot_min_time, layout=l,size=(1100,500), bottom_margin = 15Plots.mm, top_margin = 5Plots.mm,  left_margin = 15Plots.mm ,
                       xtickfontsize=14,ytickfontsize=14 , 
                       xlabelfontsize=16,ylabelfontsize=16, titlefontsize=20,legendfontsize=14)

total_local_reduction_time = (sum(vec(sum(NLP_time_max,dims=2))./60) + sum(total_min_time./60))/60;
println("Total optimization time across all local reduction iterations: ", round(total_local_reduction_time,digits=2), " hours.")

savefig(plot_times,BASE_SAVE_FOLDER*"//Plots//Selected_Results_Plots//0_Optim_time.svg")

# ---------------------------------------------------------
# ------------------------- EMT time ----------------------
# ---------------------------------------------------------


sim_time_max    = vec(EMT_time_max)           ; # replace first sim with maximum sim value, first sim has julia+matlab compilation time noise (first sim of all the process takes the longer)
sim_time_max[1] = maximum(EMT_time_max[2:end]);
sim_time_min = cat( vec([mean(i) for i in EMT_min_time_All_GFLs_but_Gen1]),vec([mean(i) for i in EMT_min_time_All_GFLs_but_Gen2]),vec([mean(i) for i in EMT_min_time_All_GFLs_but_Gen3]),vec([mean(i) for i in EMT_min_time_All_GFMs_but_Gen1]),vec([mean(i) for i in EMT_min_time_All_GFMs_but_Gen2]),vec([mean(i) for i in EMT_min_time_All_GFMs_but_Gen3]),dims=1)

mean_max_sim_time,std_max_sim_time = mean(sim_time_max) , std(sim_time_max)
mean_min_sim_time,std_min_sim_time = mean(sim_time_min) , std(sim_time_min)

println("Mean and Std of EMT simulation time for Max-phase ROCP: ", round(mean_max_sim_time,digits=2), " ± ", round(std_max_sim_time,digits=2), " seconds.")
println("Mean and Std of EMT simulation time for Min-phase ROCP: ", round(mean_min_sim_time,digits=2), " ± ", round(std_min_sim_time,digits=2), " seconds.")