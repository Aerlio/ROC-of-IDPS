run("Matlab_1_Generate_Case_Params.m")

fixed_step_size   = 1e-4    ;
load( "Generated_files//tmp//tmp_1.mat")
jump_time       = H/10                      ;
ref_upd_time    = ref_upd_time + jump_time  ;
H               = H +  jump_time            ;

if phase=="Max"
    pause(0.1)
    simStates = Simulink.SimulationInput("Simulation");
    simStates = setModelParameter(simStates, 'StartTime',num2str(0.0), 'StopTime' ,num2str(H), 'SaveFinalState','off','SaveOperatingPoint','off');
    out = sim(simStates);
    pause(0.1)
    display("Simulation done")
    pause(0.1)
    
    var_names = ["p_ref","q_ref","Pload","Qload","p_ibr","q_ibr","w_dq","ixmag","vmag","score"];
    Time_ML = get(out.logsout,var_names(1)).Values.Time - jump_time;
    for var=var_names
        var_value = get(out.logsout,var).Values.Data             ;
        if size(var_value,1) ~= 1   var_value=squeeze(var_value)';
        else                        var_value=squeeze(var_value) ; end
        if var == var_names(1) Vars_ML = var_value               ;
        else                   Vars_ML = cat(2,Vars_ML,var_value); end
    end
    ML_obj = Vars_ML(end,end)   ;
    Time_ML        = Time_ML(1:10:end,:);
    Vars_ML        = Vars_ML(1:10:end,:);
else
    for optim_iter = 1:Optim_iter
        if optim_iter ~= 1 
            load( "Generated_files//tmp//tmp_"+string(optim_iter)+".mat") 
            jump_time       = H/10                      ;
            ref_upd_time    = ref_upd_time + jump_time  ;
            H               = H +  jump_time            ;
        end

        pause(0.1)
        simStates = Simulink.SimulationInput("Simulation");
        simStates = setModelParameter(simStates, 'StartTime',num2str(0.0), 'StopTime' ,num2str(H), 'SaveFinalState','off','SaveOperatingPoint','off');
        out = sim(simStates);
        pause(0.1)
        display("Simulation done")
        pause(0.1)

        var_names = ["p_ref","q_ref","Pload","Qload","p_ibr","q_ibr","w_dq","ixmag","vmag","score"];
        Time_ML_mnv = get(out.logsout,var_names(1)).Values.Time - jump_time;
        for var=var_names
            var_value = get(out.logsout,var).Values.Data                     ;
            if size(var_value,1) ~= 1   var_value=squeeze(var_value)'        ;
            else                        var_value=squeeze(var_value)         ; end
            if var == var_names(1) Vars_ML_mnv = var_value                   ;
            else                   Vars_ML_mnv = cat(2,Vars_ML_mnv,var_value); end
        end
        ML_obj_mnv = Vars_ML_mnv(end,end)   ;
        
        if optim_iter == 1
            Time_ML = Time_ML_mnv(1:10:end,:);
            Vars_ML = Vars_ML_mnv(1:10:end,:);
            ML_obj = ML_obj_mnv ;
        else
            Time_ML = cat(1,Time_ML,Time_ML_mnv(1:10:end,:));
            Vars_ML = cat(1,Vars_ML,Vars_ML_mnv(1:10:end,:));
            if ML_obj_mnv >= ML_obj
                ML_obj = ML_obj_mnv ;
            end
        end
    end
end

save("Generated_files//tmp//"+phase+"_"+string(Optim_iter)+"_tmp.mat","Time_ML", "Vars_ML", "ML_obj");