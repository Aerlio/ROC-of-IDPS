include("Julia_0_0_General_Calls.jl") ;

global Optim_iter = 7
global files = readdir(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//")

# ---------------------------------------------------------
# ---------------- Read MAX phase obj. values -------------
# ---------------------------------------------------------

global worst_scores = zeros(Optim_iter,1);
global max_scores   = zeros(Optim_iter,10);

for optim_iter in 1:Optim_iter
    global max_files_optim_iter = filter(f -> startswith(f, string(optim_iter)*"_Max_"), files)
    global max_files_optim_iter = filter(f -> endswith(f, "_for_sim.mat"), max_files_optim_iter)
    global max_filenames_optim_iter = [max_file[1:end-12] for max_file in max_files_optim_iter]
    global max_filenames_optim_iter = sort(max_filenames_optim_iter, by=length)


    for i in 1:length(max_filenames_optim_iter)
        global (_ ,_,_,_,_,_,_,_,_,_,_,_,_,_,_, max_score) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*max_filenames_optim_iter[i]);
        if i == 1
            global worst_scores[optim_iter] = copy(max_score)
        else
            global max_scores[optim_iter,i-1] = copy(max_score)
        end
    end
end

global min_score                    = []
global min_score_All_GFMs           = []
global min_score_All_GFLs_but_Gen1  = []
global min_score_All_GFLs_but_Gen2  = []
global min_score_All_GFLs_but_Gen3  = []
global min_score_All_GFMs_but_Gen1  = []
global min_score_All_GFMs_but_Gen2  = []
global min_score_All_GFMs_but_Gen3  = []

global Conf_nums = [7]
global N1_ratio  = zeros(Optim_iter,Nibr)
global N2_ratio  = zeros(Optim_iter,Nibr)
global N1_ratio[1,:] = 0.01*ones(Nibr)
global N2_ratio[1,:] = 0.01*ones(Nibr)

for optim_iter in 2:Optim_iter

    global min_files_optim_iter = filter(f -> startswith(f, string(optim_iter-1)*"_Min_"), files)
    global min_files_optim_iter = filter(f -> endswith(f, "_for_sim.mat"), min_files_optim_iter)
    global min_filenames_optim_iter = [min_file[1:end-12] for min_file in min_files_optim_iter]
    global min_filenames_optim_iter = sort(min_filenames_optim_iter, by=length)

    global All_GFMs_file = filter(s -> occursin("All_GFMs-", s), min_filenames_optim_iter.*"-")
    global All_GFLs_but_Gen1_file = filter(s -> occursin("All_GFLs_but_Gen1-", s), min_filenames_optim_iter.*"-")
    global All_GFLs_but_Gen2_file = filter(s -> occursin("All_GFLs_but_Gen2-", s), min_filenames_optim_iter.*"-")
    global All_GFLs_but_Gen3_file = filter(s -> occursin("All_GFLs_but_Gen3-", s), min_filenames_optim_iter.*"-")
    global All_GFMs_but_Gen1_file = filter(s -> occursin("All_GFMs_but_Gen1-", s), min_filenames_optim_iter.*"-")
    global All_GFMs_but_Gen2_file = filter(s -> occursin("All_GFMs_but_Gen2-", s), min_filenames_optim_iter.*"-")
    global All_GFMs_but_Gen3_file = filter(s -> occursin("All_GFMs_but_Gen3-", s), min_filenames_optim_iter.*"-")

    global (_,_,_,_,_,_,_,_,_,_,_,ibr_conf_num,_,n1_ratio,n2_ratio,score_min) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*min_filenames_optim_iter[1])
    if size(All_GFMs_file,1)>0          global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,score_7) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*All_GFMs_file[1][1:end-1]);          else score_7 = nothing ; end 
    if size(All_GFLs_but_Gen1_file,1)>0 global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,score_6) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*All_GFLs_but_Gen1_file[1][1:end-1]); else score_6 = nothing ; end 
    if size(All_GFLs_but_Gen2_file,1)>0 global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,score_5) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*All_GFLs_but_Gen2_file[1][1:end-1]); else score_5 = nothing ; end 
    if size(All_GFLs_but_Gen3_file,1)>0 global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,score_4) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*All_GFLs_but_Gen3_file[1][1:end-1]); else score_4 = nothing ; end 
    if size(All_GFMs_but_Gen1_file,1)>0 global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,score_3) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*All_GFMs_but_Gen1_file[1][1:end-1]); else score_3 = nothing ; end 
    if size(All_GFMs_but_Gen2_file,1)>0 global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,score_2) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*All_GFMs_but_Gen2_file[1][1:end-1]); else score_2 = nothing ; end 
    if size(All_GFMs_but_Gen3_file,1)>0 global (_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,score_1) = read_results(BASE_SAVE_FOLDER*"//1_Mesh_Refinement//3_matlab_files//"*All_GFMs_but_Gen3_file[1][1:end-1]); else score_1 = nothing ; end 
    
    global Conf_nums = push!(Conf_nums,ibr_conf_num)
    global N1_ratio[optim_iter,:] = copy(n1_ratio)
    global N2_ratio[optim_iter,:] = copy(n2_ratio)
    global min_score = push!(min_score,score_min) 
    if score_7!=nothing global min_score_All_GFMs          = push!(min_score_All_GFMs         ,score_7) end
    if score_6!=nothing global min_score_All_GFLs_but_Gen1 = push!(min_score_All_GFLs_but_Gen1,score_6) end       
    if score_5!=nothing global min_score_All_GFLs_but_Gen2 = push!(min_score_All_GFLs_but_Gen2,score_5) end       
    if score_4!=nothing global min_score_All_GFLs_but_Gen3 = push!(min_score_All_GFLs_but_Gen3,score_4) end       
    if score_3!=nothing global min_score_All_GFMs_but_Gen1 = push!(min_score_All_GFMs_but_Gen1,score_3) end      
    if score_2!=nothing global min_score_All_GFMs_but_Gen2 = push!(min_score_All_GFMs_but_Gen2,score_2) end     
    if score_1!=nothing global min_score_All_GFMs_but_Gen3 = push!(min_score_All_GFMs_but_Gen3,score_1) end     
end


config_strings=[    "GFM IBRs: [2,3] / GFL IBRs: [1]", 
                    "GFM IBRs: [1] / GFL IBRs: [2,3]",
                    "GFM IBRs: [1,3] / GFL IBRs: [2]", 
                    "GFM IBRs: [2] / GFL IBRs: [1,3]",
                    "GFM IBRs: [1,2] / GFL IBRs: [3]", 
                    "GFM IBRs: [3] / GFL IBRs: [1,2]", 
                    "GFM IBRs: [1,2,3]"]
config_symbols=[    :rect, :dtriangle, :diamond, :circle, :utriangle, :pentagon, :star5]
global colors_by_gen   = palette(:seaborn_bright)[1:3]' ;

l = @layout([a b])

plot_max = plot(1:Optim_iter, worst_scores                , linewidth=2, linestyle=:dash        , ylabel="Obj. function (EMT)", xlabel="N° scenarios", markersize=11 ,  color=:red  , label="Selected worst-case solution", ylim=(-0.5,14))
plot_max = scatter!(1:Optim_iter, max_scores[:,1]     , alpha=0.2                               , markersize=11 ,  color=:red  , label="Worst-case solutions across warm-starts", title="Maximisation phase", marker=:circle)
plot_max = scatter!(1:Optim_iter, max_scores[:,2:end] , alpha=0.2                               , markersize=11 ,  color=:red  , label="",  legend=:topright                                                    , marker=:circle)

plot_min = plot(1:Optim_iter-1, min_score./H_min          , linewidth=2, linestyle=:dash, color=:green   , ylabel="Obj. function (EMT)"  ,legend=:topright, xlabel="N° scenarios", label="Selected control design solution" , ylim=(0.1,0.4), title="Minimisation phase")
plot_min = plot!(1:length(min_score_All_GFMs)          , alpha=0.5 , min_score_All_GFMs./H_min           ,  label="GFM IBRs: [1,2,3]", marker=:star5, markersize=11)
plot_min = plot!(1:length(min_score_All_GFMs_but_Gen1) , alpha=0.5 , min_score_All_GFMs_but_Gen1./H_min  ,  label="GFM IBRs: [2,3] / GFL IBRs: [1]", marker=:rect, markersize=11)
plot_min = plot!(1:length(min_score_All_GFMs_but_Gen2) , alpha=0.5 , min_score_All_GFMs_but_Gen2./H_min  ,  label="GFM IBRs: [1,3] / GFL IBRs: [2]", marker=:diamond, markersize=11)
plot_min = plot!(1:length(min_score_All_GFMs_but_Gen3) , alpha=0.5 , min_score_All_GFMs_but_Gen3./H_min  ,  label="GFM IBRs: [1,2] / GFL IBRs: [3]", marker=:utriangle, markersize=11)
plot_min = plot!(1:length(min_score_All_GFLs_but_Gen1) , alpha=0.5 , min_score_All_GFLs_but_Gen1./H_min  ,  label="GFM IBRs: [1] / GFL IBRs: [2,3]", marker=:dtriangle, markersize=11)
plot_min = plot!(1:length(min_score_All_GFLs_but_Gen2) , alpha=0.5 , min_score_All_GFLs_but_Gen2./H_min  ,  label="GFM IBRs: [2] / GFL IBRs: [1,3]", marker=:circle, markersize=11)
plot_min = plot!(1:length(min_score_All_GFLs_but_Gen3) , alpha=0.5 , min_score_All_GFLs_but_Gen3./H_min  ,  label="GFM IBRs: [3] / GFL IBRs: [1,2]", marker=:pentagon, markersize=11, legend_columns=2)
plot_final = plot(plot_max,plot_min, layout=l,size=(1600,500), bottom_margin = 10Plots.mm,  left_margin = 10Plots.mm 
                       , xtickfontsize=10,ytickfontsize=10
                       , xlabelfontsize=12,ylabelfontsize=12, titlefontsize=16,legendfontsize=11)

savefig(plot_final,"C://Users//tochoa//OneDrive - Imperial College London//Desktop//PhD ICL//Start_of_2025//PaperImages//3-2-MinMaxRES.svg")

l = @layout([a b])
eta_1 = plot(0:(Optim_iter-1), N1_ratio , color=colors_by_gen, label=permutedims(["IBR 1","IBR 2","IBR 3"]), ylim=(-0.001,0.025), title="η₁")
eta_1 = scatter!(0:(Optim_iter-1), N1_ratio , color=colors_by_gen, markershape=config_symbols[Conf_nums], markersize=8, alpha=0.7, label="", xlabel="N° scenarios", ylabel="Param. value")
eta_2 = plot(0:(Optim_iter-1), N2_ratio , color=colors_by_gen, ylim=(-0.005,0.065),legend=false, title="η₂")
eta_2 = scatter!(0:(Optim_iter-1), N2_ratio , color=colors_by_gen, markershape=config_symbols[Conf_nums], markersize=8, alpha=0.7, label="", xlabel="N° scenarios", ylabel="Param. value")

plot_eta=plot(eta_1,eta_2, layout=l,size=(800,250), bottom_margin = 4Plots.mm,top_margin = 4Plots.mm,  left_margin = 10Plots.mm 
                       , xtickfontsize=10,ytickfontsize=10
                       , xlabelfontsize=12,ylabelfontsize=12, titlefontsize=14,legendfontsize=11)

                       
savefig(plot_eta,"C://Users//tochoa//OneDrive - Imperial College London//Desktop//PhD ICL//Start_of_2025//PaperImages//3-3-ETA_RES.svg")
