run("Matlab_1_Generate_Case_Params.m")
fixed_step_size = 5e-5    ;

file="CSVs\"+PScase+"\IBR_Configs.csv";
id_configs =  table2array(readtable(file,opts)); 
id_configs = id_configs(1:8,2:1+Nibr);
load("Generated_files//Local_Reduction//Optimization//Selected_Results//7_Max_for_sim.mat")
best_n1_ratio = n1_ratio;
best_n2_ratio = n2_ratio;
best_y_theta  = y_theta ;
best_score    = 0 ;

for scen = 1:6
    load("Generated_files//Local_Reduction//Optimization//Selected_Results//"+scen+"_Max_for_sim.mat");

    ref_upd_time    = 0.5;
    H               = 3.0;
    jump_time       = H/10                      ;
    ref_upd_time    = ref_upd_time + jump_time  ;
    H               = H +  jump_time            ;

    y_theta  = best_y_theta;
    n1_ratio = best_n1_ratio;
    n2_ratio = best_n2_ratio;
    
    simStates = Simulink.SimulationInput("Simulation");
    simStates = setModelParameter(simStates, 'StartTime',num2str(0.0), 'StopTime' ,num2str(H), 'SaveFinalState','off','SaveOperatingPoint','off');
    out = sim(simStates);
    disp("Simulation done")

    score_res = get(out.logsout,"score").Values.Data(end);
    if score_res >= best_score
        best_score = score_res;
    end
    disp(best_score)
end

% Grids
grid_t_i    = [0.15, 0.25, 0.5, 0.5];
grid_t_f    = [0.3 , 0.5 , 1.0, 3.0];
grid_ratio  = [0.01, 0.05, 0.10];
grid_config = [8, 2 , 5, 7];

% --- Build (config, ratio) combinations (like your Julia repeats) ---
% Each config repeated for all ratios
cfg_col   = repelem(grid_config(:), numel(grid_ratio), 1);         
ratio_col = repmat(grid_ratio(:), numel(grid_config), 1);          
score_to_fill = [cfg_col,ratio_col,ratio_col,zeros(size(cfg_col,1),1)];

for control_design_i = 1:size(score_to_fill,1)
    for scen = 1:6
        load("Generated_files//Local_Reduction//Optimization//Selected_Results//"+scen+"_Max_for_sim.mat");
    
        ref_upd_time    = 0.5;
        H               = 3.0;
        jump_time       = H/10                      ;
        ref_upd_time    = ref_upd_time + jump_time  ;
        H               = H +  jump_time            ;
    
        y_theta  = id_configs(score_to_fill(control_design_i,1),:)';
        n1_ratio = score_to_fill(control_design_i,2) * ones(Nibr,1);
        n2_ratio = score_to_fill(control_design_i,3) * ones(Nibr,1);
        
        simStates = Simulink.SimulationInput("Simulation");
        simStates = setModelParameter(simStates, 'StartTime',num2str(0.0), 'StopTime' ,num2str(H), 'SaveFinalState','off','SaveOperatingPoint','off');
        out = sim(simStates);
        disp("Simulation done")
    
        score_res = get(out.logsout,"score").Values.Data(end);
        if score_res >= score_to_fill(control_design_i,4)
            score_to_fill(control_design_i,4) = score_res;
        end
        disp(score_to_fill)
    end
end

%%

row_best = {'best','best','best',best_score};
res = [row_best; num2cell(score_to_fill)];

writecell(res, 'Generated_Files\benchmark.csv');