run("Matlab_1_Generate_Case_Params.m")

fixed_step_size = 5e-5    ;
load( "Generated_files//tmp//tmp_1.mat")
jump_time       = H/10                      ;
ref_upd_time    = ref_upd_time + jump_time  ;
H               = H +  jump_time            ;

if phase=="Max"
    pause(0.1)
    tic
    simStates = Simulink.SimulationInput("Simulation");
    simStates = setModelParameter(simStates, 'StartTime',num2str(0.0), 'StopTime' ,num2str(H), 'SaveFinalState','off','SaveOperatingPoint','off');
    out = sim(simStates);
    EMT_time = toc*1.0;
    pause(0.1)
    display("Simulation done")
    pause(0.1)
    
    var_names = ["p_ref","q_ref","Pload","Qload","p_ibr","q_ibr","w_dq","ixmag","vmag","score"];
    ts_EMT = get(out.logsout,var_names(1)).Values.Time - jump_time;
    for var=var_names
        var_value = get(out.logsout,var).Values.Data             ;
        if size(var_value,1) ~= 1   var_value=squeeze(var_value)';
        else                        var_value=squeeze(var_value) ; end
        if var == var_names(1) Vars_EMT = var_value               ;
        else                   Vars_EMT = cat(2,Vars_EMT,var_value); end
    end
    EMT_obj = Vars_EMT(end,end)   ;
    ts_EMT        = ts_EMT(1:10:end,:);
    Vars_EMT        = Vars_EMT(1:10:end,:);
else
    for optim_iter = 1:Optim_iter
        if optim_iter ~= 1 
            load( "Generated_files//tmp//tmp_"+string(optim_iter)+".mat") 
            jump_time       = H/10                      ;
            ref_upd_time    = ref_upd_time + jump_time  ;
            H               = H +  jump_time            ;
        end

        pause(0.1)
        tic 
        simStates = Simulink.SimulationInput("Simulation");
        simStates = setModelParameter(simStates, 'StartTime',num2str(0.0), 'StopTime' ,num2str(H), 'SaveFinalState','off','SaveOperatingPoint','off');
        out = sim(simStates);
        EMT_time_mnv = toc*1.0;
        pause(0.1)
        display("Simulation done")
        pause(0.1)

        var_names = ["p_ref","q_ref","Pload","Qload","p_ibr","q_ibr","w_dq","ixmag","vmag","score"];
        ts_EMT_mnv = get(out.logsout,var_names(1)).Values.Time - jump_time;
        for var=var_names
            var_value = get(out.logsout,var).Values.Data                     ;
            if size(var_value,1) ~= 1   var_value=squeeze(var_value)'        ;
            else                        var_value=squeeze(var_value)         ; end
            if var == var_names(1) Vars_EMT_mnv = var_value                   ;
            else                   Vars_EMT_mnv = cat(2,Vars_EMT_mnv,var_value); end
        end
        EMT_obj_mnv = Vars_EMT_mnv(end,end)   ;
        
        if optim_iter == 1
            ts_EMT = ts_EMT_mnv(1:10:end,:);
            Vars_EMT = Vars_EMT_mnv(1:10:end,:);
            EMT_obj = EMT_obj_mnv ;
            EMT_time = EMT_time_mnv;
        else
            ts_EMT = cat(1,ts_EMT,ts_EMT_mnv(1:10:end,:));
            Vars_EMT = cat(1,Vars_EMT,Vars_EMT_mnv(1:10:end,:));
            if EMT_obj_mnv >= EMT_obj
                EMT_obj = EMT_obj_mnv ;
            end
            EMT_time = cat(1,EMT_time,EMT_time_mnv);;
        end
    end
end

save("Generated_files//tmp//"+phase+"_"+string(Optim_iter)+"_tmp.mat","ts_EMT", "Vars_EMT", "EMT_obj", "EMT_time" );