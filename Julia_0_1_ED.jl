if add_ED_to_model == false
        #******************  Create model ****************** #
        global model_traj = Model(KNITRO.Optimizer);        
        set_attribute(model_traj, "strat_warm_start", 1)
        set_attribute(model_traj, "outlev", 0)
        set_attribute(model_traj, "maxit", 5000)
        set_attribute(model_traj, "ftol", 1e-2)
        set_attribute(model_traj, "out_csvinfo", 1)

        global optim_filename = WS_SAVE_FOLDER*"//Optim_csvs//"*string(length(readdir(WS_SAVE_FOLDER*"//Optim_csvs//")))*"_ED_AGC.csv"
        set_attribute(model_traj, "out_csvname", optim_filename);
        
        # Reference variables
        global p_ref_ed     = @variable(  model_traj,   p_ref_ed[n_ibr=1:Nibr]           , start = (-s_ibr_nom[n_ibr] + 2*s_ibr_nom[n_ibr])*rand());
        global q_ref_ed     = @variable(  model_traj,   q_ref_ed[n_ibr=1:Nibr]           , start = (-s_ibr_nom[n_ibr] + 2*s_ibr_nom[n_ibr])*rand());
        global p_ref_agc    = @variable(  model_traj,  p_ref_agc[n_ibr=1:Nibr]           , start = (-s_ibr_nom[n_ibr] + 2*s_ibr_nom[n_ibr])*rand());

        # ----------------------- ED Variables -----------------------------

        # State variables X
        global vsh_D_ed           = @variable(model_traj,vsh_D_ed[n_sh=1:Nsh]           , start=ones(Nsh)[n_sh]);
        global vsh_Q_ed           = @variable(model_traj,vsh_Q_ed[n_sh=1:Nsh]           , start=ones(Nsh)[n_sh]);
        global ie_D_ed            = @variable(model_traj,ie_D_ed[n_e=1:Ne]              , start=ones(Ne)[n_e]);
        global ie_Q_ed            = @variable(model_traj,ie_Q_ed[n_e=1:Ne]              , start=ones(Ne)[n_e]);
        global theta_ibr_dq_ed    = @variable(model_traj,theta_ibr_dq_ed[n_ibr=1:Nibr]  , start=zeros(Nibr)[n_ibr]);  
        global ix_d_ed            = @variable(model_traj,ix_d_ed[n_ibr=1:Nibr]          , start=ones(Nibr)[n_ibr]);
        global ix_q_ed            = @variable(model_traj,ix_q_ed[n_ibr=1:Nibr]          , start=ones(Nibr)[n_ibr]);
        global idc_ed             = @variable(model_traj,idc_ed[n_ibr=1:Nibr]           , start=(ones(Nibr)./ Vdc_n)[n_ibr]);
        
        # Fixed variables
        global vibr_d_ed          = ones(Nibr)  ;
        global vibr_q_ed          = zeros(Nibr) ;
        global vdc_ed             = Vdc_n[:,1]  ;
        global hvdc_ed            = zeros(Nibr) ;
        global hi_d_ed            = zeros(Nibr) ;
        global hi_q_ed            = zeros(Nibr) ;
        global hv_d_ed            = zeros(Nibr) ;
        global hv_q_ed            = zeros(Nibr) ;
        global hpll_ed            = zeros(Nibr) ;
        global hlpf_pll_ed        = zeros(Nibr) ;
        global hpq_d_ed           = @expression(model_traj, ix_d_ed);
        global hpq_q_ed           = @expression(model_traj, ix_q_ed);
else
        # ----------------------- ED Variables -----------------------------
        # State variables X
        global vsh_D_ed           = @expression(model_traj,[n_sh=1:Nsh]         , vsh_D[1,n_sh]               );
        global vsh_Q_ed           = @expression(model_traj,[n_sh=1:Nsh]         , vsh_Q[1,n_sh]               );
        global ie_D_ed            = @expression(model_traj,[n_e=1:Ne]           , ie_D[1,n_e]                 );
        global ie_Q_ed            = @expression(model_traj,[n_e=1:Ne]           , ie_Q[1,n_e]                 );
        global theta_ibr_dq_ed    = @expression(model_traj,[n_ibr=1:Nibr]       , theta_ibr_dq[1,n_ibr]       );  
        global ix_d_ed            = @expression(model_traj,[n_ibr=1:Nibr]       , ix_d[1,n_ibr]               );
        global ix_q_ed            = @expression(model_traj,[n_ibr=1:Nibr]       , ix_q[1,n_ibr]               );
        global idc_ed             = @expression(model_traj,[n_ibr=1:Nibr]       , idc[1,n_ibr]                );
        # Fixed variables
        global vibr_d_ed          = ones(Nibr)  ;
        global vibr_q_ed          = zeros(Nibr) ;
        global vdc_ed             = Vdc_n[:,1]  ;
        global hvdc_ed            = zeros(Nibr) ;
        global hi_d_ed            = zeros(Nibr) ;
        global hi_q_ed            = zeros(Nibr) ;
        global hv_d_ed            = zeros(Nibr) ;
        global hv_q_ed            = zeros(Nibr) ;
        global hpll_ed            = zeros(Nibr) ;
        global hlpf_pll_ed        = zeros(Nibr) ;

        fix.(vibr_d[1,:]  ,  vibr_d_ed          ; force=true);
        fix.(vibr_q[1,:]  ,  vibr_q_ed          ; force=true);
        fix.(vdc[1,:]     ,  vdc_ed             ; force=true);
        fix.(hvdc[1,:]    ,  hvdc_ed            ; force=true);
        fix.(hi_d[1,:]    ,  hi_d_ed            ; force=true);
        fix.(hi_q[1,:]    ,  hi_q_ed            ; force=true);
        for n_ibr in 1:Nibr
                if y_theta[n_ibr]==0 #GFM
                        # fix.(hv_d[1,n_ibr]      ,  hv_d_ed[n_ibr]            ; force=true); # <- our Kv_I's are zeros
                        # fix.(hv_q[1,n_ibr]      ,  hv_q_ed[n_ibr]            ; force=true);
                else #GFL
                        fix.(hpll[1,n_ibr]      ,  hpll_ed[n_ibr]            ; force=true);
                        fix.(hlpf_pll[1,n_ibr]  ,  hlpf_pll_ed[n_ibr]        ; force=true);
                        @constraint(model_traj  , hpq_d[1,n_ibr] == ix_d[1,n_ibr])          ;
                        @constraint(model_traj  , hpq_q[1,n_ibr] == ix_q[1,n_ibr])          ;
                end
        end
end

#******************  Create general model ****************** #

# Expressions
global delta_ibr_ed       = @expression(model_traj,[n_ibr=1:Nibr],theta_ibr_dq_ed[n_ibr]);
global vibr_D_ed          = @expression(model_traj,[n_ibr=1:Nibr],vibr_d_ed[n_ibr]*cos(delta_ibr_ed[n_ibr]));
global vibr_Q_ed          = @expression(model_traj,[n_ibr=1:Nibr],vibr_d_ed[n_ibr]*sin(delta_ibr_ed[n_ibr]));

global (vN_D_ed,vN_Q_ed)     = (Matrix{AffExpr}[],Matrix{AffExpr}[])
global (iL_D_ed,iL_Q_ed)     = (Matrix{AffExpr}[],Matrix{AffExpr}[])
global (iibr_D_ed,iibr_Q_ed) = (Matrix{AffExpr}[],Matrix{AffExpr}[])

for n_N in 1:Nn
        if sum(Bsh[n_N,:])!=0
                global vN_D_ed = cat(vN_D_ed, @expression(model_traj, sum(vsh_D_ed[Bsh[n_N,:].==1]   ) ), dims=1);
                global vN_Q_ed = cat(vN_Q_ed, @expression(model_traj, sum(vsh_Q_ed[Bsh[n_N,:].==1]   ) ), dims=1);
        elseif sum(Bibr[n_N,:])!=0
                global vN_D_ed = cat(vN_D_ed, @expression(model_traj, sum(vibr_D_ed[n_ibr]  for n_ibr in 1:Nibr if sum(Bibr[n_N,n_ibr])!=0) ), dims=1 );
                global vN_Q_ed = cat(vN_Q_ed, @expression(model_traj, sum(vibr_Q_ed[n_ibr]  for n_ibr in 1:Nibr if sum(Bibr[n_N,n_ibr])!=0) ), dims=1 );
        else
                println("Problem: Undefined Nodal Voltage Equation")
        end
end

for n_L in 1:Nload
        if sum(Bload[:,n_L])!=0 
                if a_p[n_L]==1 && a_q[n_L]==1
                        global iL_D_ed = cat(iL_D_ed, @expression(model_traj, (
                                                sum(a_p[n_L]*vN_D_ed[Bload[:,n_L].==1])*Pload_ini[n_L] + a_q[n_L]*sum(vN_Q_ed[Bload[:,n_L].==1])*Qload_ini[n_L]
                                                ) ), dims=1);
                        global iL_Q_ed = cat(iL_Q_ed, @expression(model_traj, (
                                                sum(a_p[n_L]*vN_Q_ed[Bload[:,n_L].==1])*Pload_ini[n_L] - a_q[n_L]*sum(vN_D_ed[Bload[:,n_L].==1])*Qload_ini[n_L]
                                                ) ), dims=1);
                elseif c_p[n_L]==1 && c_q[n_L]==1
                        global iL_D_ed = cat(iL_D_ed, @expression(model_traj, (
                                                (sum(c_p[n_L]*vN_D_ed[Bload[:,n_L].==1])*Pload_ini[n_L] + c_q[n_L]*sum(vN_Q_ed[Bload[:,n_L].==1])*Qload_ini[n_L])/
                                                                (sum(vN_D_ed[Bload[:,n_L].==1])^2+sum(vN_Q_ed[Bload[:,n_L].==1])^2)
                                                ) ), dims=1);
                        global iL_Q_ed = cat(iL_Q_ed, @expression(model_traj, (
                                                (sum(c_p[n_L]*vN_Q_ed[Bload[:,n_L].==1])*Pload_ini[n_L] - c_q[n_L]*sum(vN_D_ed[Bload[:,n_L].==1])*Qload_ini[n_L])/
                                                                (sum(vN_D_ed[Bload[:,n_L].==1])^2+sum(vN_Q_ed[Bload[:,n_L].==1])^2)
                                                ) ), dims=1);
                end
        else
                println("Problem: Undefined Load current Equation")
        end
end

for n_ibr in 1:Nibr
        if sum(Be[Bibr[:,n_ibr].==1,:].!=0)!=0 && sum(Bload[Bibr[:,n_ibr].==1,:])!=0
                global iibr_D_ed = cat(iibr_D_ed, @expression(model_traj, sum((Be*ie_D_ed)[Bibr[:,n_ibr].==1]) + sum((Bload*iL_D_ed)[Bibr[:,n_ibr].==1]) ), dims= 1 );
                global iibr_Q_ed = cat(iibr_Q_ed, @expression(model_traj, sum((Be*ie_Q_ed)[Bibr[:,n_ibr].==1]) + sum((Bload*iL_Q_ed)[Bibr[:,n_ibr].==1]) ), dims= 1 );
        
        elseif sum(Be[Bibr[:,n_ibr].==1,:].!=0)!=0
                global iibr_D_ed = cat(iibr_D_ed, @expression(model_traj, sum((Be*ie_D_ed)[Bibr[:,n_ibr].==1]) ), dims= 1 );
                global iibr_Q_ed = cat(iibr_Q_ed, @expression(model_traj, sum((Be*ie_Q_ed)[Bibr[:,n_ibr].==1]) ), dims= 1 );
        else
                println("Problem: Undefined IBR Current Equation")
        end
end

global iibr_d_ed   = @expression(model_traj,[n_ibr=1:Nibr], iibr_D_ed[n_ibr]*cos(-delta_ibr_ed[n_ibr])-iibr_Q_ed[n_ibr]*sin(-delta_ibr_ed[n_ibr]));
global iibr_q_ed   = @expression(model_traj,[n_ibr=1:Nibr], iibr_D_ed[n_ibr]*sin(-delta_ibr_ed[n_ibr])+iibr_Q_ed[n_ibr]*cos(-delta_ibr_ed[n_ibr]));
global p_ibr_ed    = @expression(model_traj,[n_ibr=1:Nibr], (vibr_d_ed[n_ibr]*iibr_d_ed[n_ibr]));
global q_ibr_ed    = @expression(model_traj,[n_ibr=1:Nibr], (-vibr_d_ed[n_ibr]*iibr_q_ed[n_ibr]));

global w_dq_ed = w_0*ones(Nibr);
global vx_d_ed = @expression(model_traj,[n_ibr=1:Nibr], vibr_d_ed[n_ibr] + R_f[n_ibr]*ix_d_ed[n_ibr] - w_dq_ed[n_ibr]*L_f[n_ibr]*ix_q_ed[n_ibr]);
global vx_q_ed = @expression(model_traj,[n_ibr=1:Nibr],                    R_f[n_ibr]*ix_q_ed[n_ibr] + w_dq_ed[n_ibr]*L_f[n_ibr]*ix_d_ed[n_ibr]);

global pdc_ibr_ed  = @expression(model_traj,[n_ibr=1:Nibr], ((vx_d_ed[n_ibr]*ix_d_ed[n_ibr]+vx_q_ed[n_ibr]*ix_q_ed[n_ibr]))); #vdc*iz = vdc*(m_d*ix_d+m_q+ix_q) =vdc*(1/vdc)(vx_d*ix_d+vx_q+ix_q)

# ------------------ Steady-state df/dt = 0 ---------------------------

# Dynamics Shunt voltages
for n_sh in 1:Nsh 
        if sum(Bload[Bsh[:,n_sh].==1,:])!=0
                @constraint(model_traj, w_b*(1 ./C_sh[n_sh])*( -sum((Be*ie_D_ed)[Bsh[:,n_sh].==1]) - sum((Bload*iL_D_ed)[Bsh[:,n_sh].==1]) )  + w_b*vsh_Q_ed[n_sh]==0);
                @constraint(model_traj, w_b*(1 ./C_sh[n_sh])*( -sum((Be*ie_Q_ed)[Bsh[:,n_sh].==1]) - sum((Bload*iL_Q_ed)[Bsh[:,n_sh].==1]) )  - w_b*vsh_D_ed[n_sh]==0);
        else
                @constraint(model_traj, w_b*(1 ./C_sh[n_sh])*( -sum((Be*ie_D_ed)[Bsh[:,n_sh].==1]) )  + w_b*vsh_Q_ed[n_sh]==0);
                @constraint(model_traj, w_b*(1 ./C_sh[n_sh])*( -sum((Be*ie_Q_ed)[Bsh[:,n_sh].==1]) )  - w_b*vsh_D_ed[n_sh]==0);
        end
end

# Dynamics Line Currents
for n_e in 1:Ne
        if sum(Be[:,n_e].!=0)!=0
                @constraint(model_traj, w_b*(1 ./Le[n_e])*(  transpose(Be[:,n_e][Be[:,n_e].!=0])*vN_D_ed[Be[:,n_e].!=0] - Re[n_e]*ie_D_ed[n_e]  )  + w_b*ie_Q_ed[n_e]  == 0);
                @constraint(model_traj, w_b*(1 ./Le[n_e])*(  transpose(Be[:,n_e][Be[:,n_e].!=0])*vN_Q_ed[Be[:,n_e].!=0] - Re[n_e]*ie_Q_ed[n_e]  )  - w_b*ie_D_ed[n_e]  == 0);
        else 
                println("Problem: Undefined Deriv. Line current Equation")
        end
end

# Dynamics Physical variables: RL-filter
@constraint(model_traj,[n_ibr=1:Nibr],w_b*(1 ./L_f[n_ibr])*(vx_d_ed[n_ibr]-vibr_d_ed[n_ibr]-R_f[n_ibr]*ix_d_ed[n_ibr])
                                                        +w_b*w_dq_ed[n_ibr]*ix_q_ed[n_ibr]   == 0 );
@constraint(model_traj,[n_ibr=1:Nibr],w_b*(1 ./L_f[n_ibr])*(vx_q_ed[n_ibr]                 -R_f[n_ibr]*ix_q_ed[n_ibr])
                                                        -w_b*w_dq_ed[n_ibr]*ix_d_ed[n_ibr]   == 0 );
# Dynamics Physical variables: C-filter
@constraint(model_traj,[n_ibr=1:Nibr],w_b*(1 ./C_f[n_ibr])*(ix_d_ed[n_ibr]-iibr_d_ed[n_ibr])
                                                                                             ==0 );
@constraint(model_traj,[n_ibr=1:Nibr],w_b*(1 ./C_f[n_ibr])*(ix_q_ed[n_ibr]-iibr_q_ed[n_ibr])
                                                        -w_b*w_dq_ed[n_ibr]*vibr_d_ed[n_ibr] ==0 ); 

@constraint(model_traj,[n_ibr=1:Nibr], pdc_ibr_ed[n_ibr] == Vdc_n[n_ibr]*idc_ed[n_ibr]);

if add_ED_to_model == false
        global m_p_ibr_ed   = p_ref_ed;
        global m_q_ibr_ed   = q_ref_ed;
        global score_ed     = zeros(1);
else
        @constraint(model_traj,[n_ibr=1:Nibr],m_p_ibr[1,n_ibr]   == p_ref_ed[n_ibr] );
        @constraint(model_traj,[n_ibr=1:Nibr],m_q_ibr[1,n_ibr]   == q_ref_ed[n_ibr] );
        fix.(score[1] , zeros(1) ; force=true);
end


# -----------------------------------------------------------------------------
# ------------------------- DISPATCH CONSTRAINTS ------------------------------
# -----------------------------------------------------------------------------

@constraint(model_traj,[n_ibr=1:Nibr],       p_ref_ed[n_ibr]^2 + q_ref_ed[n_ibr]^2 <= s_ibr_nom[n_ibr]^2  );
@constraint(model_traj,[n_ibr=1:Nibr],       ix_d_ed[n_ibr]^2  + ix_q_ed[n_ibr]^2  <= s_ibr_nom[n_ibr]^2 );

@constraint(model_traj,[n_N=1:Nn]    ,       vN_D_ed[n_N]^2 + vN_Q_ed[n_N]^2        <= 1.1^2);
@constraint(model_traj,[n_N=1:Nn]    ,       vN_D_ed[n_N]^2 + vN_Q_ed[n_N]^2        >= 0.9^2);
@constraint(model_traj,[n_ibr=1:Nibr],       p_ibr_ed[n_ibr]   == p_ref_ed[n_ibr]    );
@constraint(model_traj,[n_ibr=1:Nibr],       q_ibr_ed[n_ibr]   == q_ref_ed[n_ibr]    );

# ---- Avoid redundant results
@constraint(model_traj, [n_ibr=1:Nibr] , theta_ibr_dq_ed[n_ibr]  >= -pi );
@constraint(model_traj, [n_ibr=1:Nibr] , theta_ibr_dq_ed[n_ibr]  <=  pi );
fix.(theta_ibr_dq_ed[1],0;force=true)

# ------------------------------------------------
# --------------- AGC  re-dispatch ---------------
# ------------------------------------------------

@constraint(model_traj,[n_ibr=1:Nibr], p_ref_agc[n_ibr]^2 + q_ref_ed[n_ibr]^2  <= s_ibr_nom[n_ibr]^2 );
@constraint(model_traj,                 sum(p_ref_agc .- p_ref_ed) == (Pload_fin[load_id] - Pload_ini[load_id]) );

if separate_max == true
        if      (jump_sense=="up")     @constraint(model_traj,[n_ibr=1:Nibr], p_ref_agc[n_ibr] - p_ref_ed[n_ibr] >= 0 );
        elseif  (jump_sense=="dn")     @constraint(model_traj,[n_ibr=1:Nibr], p_ref_agc[n_ibr] - p_ref_ed[n_ibr] <= 0 ); end
else
        for n_ibr_i in 1:(Nibr-1)
                for n_ibr_j in (n_ibr_i+1):Nibr
                        @constraint(model_traj, (p_ref_agc[n_ibr_i] - p_ref_ed[n_ibr_i])*(p_ref_agc[n_ibr_j] - p_ref_ed[n_ibr_j]) >= 0);
                        @constraint(model_traj, (p_ref_agc[n_ibr_i] - p_ref_ed[n_ibr_i])*(p_ref_agc[n_ibr_j] - p_ref_ed[n_ibr_j]) >= 0);
                end
        end
end

if add_ED_to_model == false
        
        # -------------------- Variables container ----------------------
        global Vars_ed = @expression(model_traj, cat( vsh_D_ed',vsh_Q_ed',ie_D_ed',ie_Q_ed',ix_d_ed',ix_q_ed',vibr_d_ed',vibr_q_ed',vdc_ed',idc_ed',m_p_ibr_ed', m_q_ibr_ed', 
                                                      hvdc_ed',hv_d_ed',hv_q_ed',hi_d_ed',hi_q_ed',hpll_ed',hlpf_pll_ed',hpq_d_ed',hpq_q_ed',theta_ibr_dq_ed',score_ed',dims=2));
        @objective(model_traj, Max,       sum(ws_ix.*((ix_d_ed.^2  .+ ix_q_ed.^2)./(s_ibr_nom.^2))) +
                                          sum(ws_agc.*(p_ref_agc./s_ibr_nom)) 
                                        );
        optimize!(model_traj)

        global info_opt=readdlm(optim_filename, ';');
        global info_opt=DataFrame(info_opt[2:end,:],info_opt[1,:]);
        
        println("(ws with ED) knot feasiblity: "      , info_opt.slv_feas[1]      );
end