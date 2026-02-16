# ------------------------------------------------------------------------
# ---------------------------- User settings -----------------------------
# ------------------------------------------------------------------------

global BASE_SAVE_FOLDER = "Generated_files//Local_Reduction";

global Max_Scenarios = 7                 ; # Number of scenarios to generate (stopping criterion)
global NLP_max       = 10                ; # Number of initial points to try 
global NLP_min       = [7,6,5,4,3,2,1]   ; # Configurations to be tested in MIN phase (indices of configs_csv)

# Reference update time and horizon
global ( ref_upd_time_max , H_max  ) = ( 0.5 , 1.0  ) ;
global ( ref_upd_time_min , H_min  ) = ( 0.5 , 3.0  ) ;

# Maximum and minimum parameters
global (n1_ratio_min, n1_ratio_max) =  1e-3, 0.1 ;
global (n2_ratio_min, n2_ratio_max) =  1e-3, 0.1 ;

# Initial control design
global (n1_ratio_ini , n2_ratio_ini , ibr_conf_num_ini) = ( 0.01*ones(3) , 0.01*ones(3) , 7) ;

#-------------------------------------------------------------------
#------------------------- Load Libraries --------------------------
#-------------------------------------------------------------------

println("Hello world...")
using LinearAlgebra, JuMP, DelimitedFiles, MAT, KNITRO, DataFrames, Random, PyCall, Plots, Statistics, Interpolations, Colors, CSV, Distributions
println("Libraries loaded...");

# ----------------------- Connect to Matlab ------------------------
if (@isdefined matlab_engine) == false
        println("Connecting to matlab...")
        global matlab_engine = pyimport("matlab.engine")
        global eng = matlab_engine.connect_matlab()                     # <- Connect
        println("Connected to matlab session")
        eng.eval("run('Matlab_1_Generate_Case_Params.m')", nargout=0)   # <- Run Case Params
        println("Case parameters updated")
else
        println("Connected to matlab session")
        eng.eval("run('Matlab_1_Generate_Case_Params.m')", nargout=0)   # <- Run Case Params
        println("Case parameters updated")
end

#-------------------------------------------------------------------
#------------------ Reading Case parameters Data  ------------------
#-------------------------------------------------------------------

global Nx = 85   ;   # Number of state variables

params = matopen("Generated_files//Case_params.mat");
variables_of_interest=["PScase", 
                        "w_b","w_0","Be","Bibr","Bload","Bsh","Nn","Nsh","Nload","Ne","Nibr",
                        "Re","Le", "C_sh","a_p","b_p","c_p","a_q","b_q","c_q",
                        "R_f","L_f","C_f","Cdc","Vdc_n","ix_ibr_max","s_ibr_nom" ,
                        "w_c","tau_dc","Kdc_P","Kdc_I","Ki_P","Ki_I","Kpll_P","Kpll_I","Kpq_P","Kpq_I","Kv_P","Kv_I","tau_lpf_pll",
                        "Pload_max", "Pload_min", "Qload_max", "Qload_min", "Pdif_max", "Qdif_max"];
global (                 PScase,   
                         w_b , w_0 , Be , Bibr , Bload  , Bsh , Nn , Nsh, Nload  , Ne , Nibr, 
                         Re  , Le , C_sh , a_p , b_p , c_p , a_q , b_q , c_q, 
                         R_f , L_f ,C_f  , Cdc , Vdc_n , ix_ibr_max , s_ibr_nom  ,
                         w_c , tau_dc , Kdc_P , Kdc_I , Ki_P , Ki_I , Kpll_P , Kpll_I , Kpq_P , Kpq_I , Kv_P , Kv_I  , tau_lpf_pll,
                         Pload_max, Pload_min, Qload_max, Qload_min, Pdif_max, Qdif_max)=
                        [read(params,i) for i in variables_of_interest];
close(params)

global Nn,Nsh,Nload,Ne,Nibr=Int(Nn),Int(Nsh),Int(Nload),Int(Ne),Int(Nibr);
global Be,Bibr,Bload,Bsh=Int.(Be),Int.(Bibr),Int.(Bload),Int.(Bsh);
global load_buses = unique(i[1] for i in findall(!iszero, Bload));
global configs_csv = readdlm("CSVs//"*PScase*"//IBR_Configs.csv", ',', Any, '\n')       ;
global ID_configs  = String.(configs_csv[2:end,1])                                      ;  
global IBR_configs = Int.(configs_csv[2:end,2:end])                             ;
function load_ibrs(ibr_conf_num)
        global id_config    = ID_configs[ibr_conf_num]                                  ;
        global y_theta      = IBR_configs[ibr_conf_num,:]                               ;
        println("Working on configuration: ", id_config)
end

function custom_mesh(mesh_loop, phase , Optim_iter , ref_upd_time , H )
        if     (mesh_loop == 1)
                global break_times_mnv            =  cat(       range(0.0               , 1e-1                  ; length =  Int(round((            1e-1         )/1e-2  )+1) )[1:end-1]  ,
                                                                range(1e-1              , ref_upd_time          ; length =  Int(round((   ref_upd_time - 1e-1   )/5e-2  )+1) )           ,
                                                                range(ref_upd_time      , ref_upd_time+1e-1     ; length =  Int(round((            1e-1         )/1e-2  )+1) )[1:end-1]  ,
                                                                range(ref_upd_time+1e-1 , H                     ; length =  Int(round(( H - ref_upd_time - 1e-1 )/5e-2  )+1) )           , dims=1)
        else
                global break_times_mnv            =  cat(       range(0.0               , 1e-5                  ; length =  Int(round((        1e-5             )/1e-6  )+1) )[1:end-1]  ,
                                                                range(1e-5              , 1e-4                  ; length =  Int(round((        1e-4 - 1e-5      )/1e-5  )+1) )[1:end-1]  ,
                                                                range(1e-4              , 1e-3                  ; length =  Int(round((        1e-3 - 1e-4      )/1e-4  )+1) )[1:end-1]  ,
                                                                range(1e-3              , 1e-2                  ; length =  Int(round((        1e-2 - 1e-3      )/1e-3  )+1) )[1:end-1]  ,
                                                                range(1e-2              , 1e-1                  ; length =  Int(round((        1e-1 - 1e-2      )/1e-2  )+1) )[1:end-1]  ,
                                                                range(1e-1              , ref_upd_time          ; length =  Int(round((   ref_upd_time - 1e-1   )/5e-2  )+1) )           ,
                                                                range(ref_upd_time      , ref_upd_time+1e-5     ; length =  Int(round((        1e-5             )/1e-6  )+1) )[1:end-1]  ,
                                                                range(ref_upd_time+1e-5 , ref_upd_time+1e-4     ; length =  Int(round((        1e-4 - 1e-5      )/1e-5  )+1) )[1:end-1]  ,
                                                                range(ref_upd_time+1e-4 , ref_upd_time+1e-3     ; length =  Int(round((        1e-3 - 1e-4      )/1e-4  )+1) )[1:end-1]  ,
                                                                range(ref_upd_time+1e-3 , ref_upd_time+1e-2     ; length =  Int(round((        1e-2 - 1e-3      )/1e-3  )+1) )[1:end-1]  ,
                                                                range(ref_upd_time+1e-2 , ref_upd_time+1e-1     ; length =  Int(round((        1e-1 - 1e-2      )/1e-2  )+1) )[1:end-1]  ,
                                                                range(ref_upd_time+1e-1 , H                     ; length =  Int(round(( H - ref_upd_time - 1e-1 )/5e-2  )+1) )           , dims=1)
        end
        global break_times_mnv    = round.( break_times_mnv , digits = 6 )
        if phase == "Max"
                global break_times =         copy(break_times_mnv)                ;
                println("Number of collocation points: ", length(break_times_mnv))
                return (break_times)
        else
                global break_times = repeat( copy(break_times_mnv) , Optim_iter ) ;
                println("Number of collocation points: ", length(break_times))
                return (break_times_mnv,break_times)
        end
end

#-------------------------------------------------------------------
#------------------------ Model Generator --------------------------
#-------------------------------------------------------------------
function model_generator(model_traj, create_variables)

        #-------------------------------------------------------------------
        #-------------------------- Variables ------------------------------
        #-------------------------------------------------------------------

        if create_variables == true

                global score                    = @variable(model_traj,         score[k=1:K]                        )
                
                # ............Power System............
                # Shunt voltages
                global vsh_D                    = @variable(model_traj,         vsh_D[k=1:K,n_sh=1:Nsh]             )
                global vsh_Q                    = @variable(model_traj,         vsh_Q[k=1:K,n_sh=1:Nsh]             )
                # Line Currents
                global ie_D                     = @variable(model_traj,         ie_D[k=1:K,n_e=1:Ne]                )
                global ie_Q                     = @variable(model_traj,         ie_Q[k=1:K,n_e=1:Ne]                )

                # ............IBRs............
                # Local dq theta
                global theta_ibr_dq             = @variable(model_traj,         theta_ibr_dq[k=1:K,n_ibr=1:Nibr]    )  

                # IBR Physical variables
                global vibr_d                   = @variable(model_traj,         vibr_d[k=1:K,n_ibr=1:Nibr]          )
                global vibr_q                   = @variable(model_traj,         vibr_q[k=1:K,n_ibr=1:Nibr]          )
                global ix_d                     = @variable(model_traj,         ix_d[k=1:K,n_ibr=1:Nibr]            )
                global ix_q                     = @variable(model_traj,         ix_q[k=1:K,n_ibr=1:Nibr]            )
                global vdc                      = @variable(model_traj,         vdc[k=1:K,n_ibr=1:Nibr]             )

                # DC current controller
                global idc                      = @variable(model_traj,         idc[k=1:K,n_ibr=1:Nibr]             )
                global hvdc                     = @variable(model_traj,         hvdc[k=1:K,n_ibr=1:Nibr]            )
                
                # Power measurement filters
                global m_p_ibr                  = @variable(model_traj,         m_p_ibr[k=1:K,n_ibr=1:Nibr]         )
                global m_q_ibr                  = @variable(model_traj,         m_q_ibr[k=1:K,n_ibr=1:Nibr]         )

                # ............IBRs' controllers...........
                # Current Controller
                global hi_d                     = @variable(model_traj,         hi_d[k=1:K,n_ibr=1:Nibr]            )
                global hi_q                     = @variable(model_traj,         hi_q[k=1:K,n_ibr=1:Nibr]            )
                
                # GFM or GFL variables
                global (hv_d,hv_q,hpll,hpq_d,hpq_q,hlpf_pll) = (Vector{VariableRef}[],Vector{VariableRef}[],Vector{VariableRef}[],Vector{VariableRef}[],Vector{VariableRef}[],Vector{VariableRef}[]);
                for n_ibr in 1:Nibr
                        if y_theta[n_ibr]==0 
                                # ----- GFM -----
                                # global hv_d       =   cat(hv_d,      @variable(model_traj,[k=1:K],            )        , dims=[1,2][min(n_ibr,2)]);
                                # global hv_q       =   cat(hv_q,      @variable(model_traj,[k=1:K],            )        , dims=[1,2][min(n_ibr,2)]); # <- our Kv_I's are zero
                                global hv_d       =   cat(hv_d,      zeros(K,1) ,dims=[1,2][min(n_ibr,2)]);
                                global hv_q       =   cat(hv_q,      zeros(K,1) ,dims=[1,2][min(n_ibr,2)]);
                                global hpll       =   cat(hpll,      zeros(K,1) ,dims=[1,2][min(n_ibr,2)]);
                                global hlpf_pll   =   cat(hlpf_pll,  zeros(K,1) ,dims=[1,2][min(n_ibr,2)]);
                                global hpq_d      =   cat(hpq_d,     zeros(K,1) ,dims=[1,2][min(n_ibr,2)]);
                                global hpq_q      =   cat(hpq_q,     zeros(K,1) ,dims=[1,2][min(n_ibr,2)]);
                        elseif y_theta[n_ibr]==1           
                                # ----- GFL -----
                                global hpll       =   cat(hpll,      @variable(model_traj,[k=1:K],            )        , dims=[1,2][min(n_ibr,2)]);
                                global hlpf_pll   =   cat(hlpf_pll,  @variable(model_traj,[k=1:K],            )        , dims=[1,2][min(n_ibr,2)]);
                                global hpq_d      =   cat(hpq_d,     @variable(model_traj,[k=1:K],            )        , dims=[1,2][min(n_ibr,2)]);
                                global hpq_q      =   cat(hpq_q,     @variable(model_traj,[k=1:K],            )        , dims=[1,2][min(n_ibr,2)]);
                                global hv_d       =   cat(hv_d,      zeros(K,1) ,dims=[1,2][min(n_ibr,2)]);
                                global hv_q       =   cat(hv_q,      zeros(K,1) ,dims=[1,2][min(n_ibr,2)]);
                        end
                end
        else
                println("Using globally defined variables in model...")
        end

        #-------------------------------------------------------------------
        #-------------------------- Expressions ----------------------------
        #-------------------------------------------------------------------

        # dq-to-DQ IBR model -> Power System
        global delta_ibr        = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],theta_ibr_dq[k,n_ibr]-theta_DQ[k]);
        global vibr_D           = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],vibr_d[k,n_ibr]*cos(delta_ibr[k,n_ibr])-vibr_q[k,n_ibr]*sin(delta_ibr[k,n_ibr]));
        global vibr_Q           = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],vibr_d[k,n_ibr]*sin(delta_ibr[k,n_ibr])+vibr_q[k,n_ibr]*cos(delta_ibr[k,n_ibr]));

        #............Power System............

        # Nodal Voltages
        global (vN_D,vN_Q) = (Matrix{AffExpr}[],Matrix{AffExpr}[])
        for n_N in 1:Nn
                if sum(Bsh[n_N,:])!=0
                        global vN_D = cat(vN_D, @expression(model_traj,[k=1:K], sum(vsh_D[k,Bsh[n_N,:].==1]   )                                       ), dims=[1,2][min(n_N,2)]);
                        global vN_Q = cat(vN_Q, @expression(model_traj,[k=1:K], sum(vsh_Q[k,Bsh[n_N,:].==1]   )                                       ), dims=[1,2][min(n_N,2)]);
                elseif sum(Bibr[n_N,:])!=0
                        global vN_D = cat(vN_D, @expression(model_traj,[k=1:K], sum(vibr_D[k,n_ibr]  for n_ibr in 1:Nibr if sum(Bibr[n_N,n_ibr])!=0) ), dims=[1,2][min(n_N,2)]);
                        global vN_Q = cat(vN_Q, @expression(model_traj,[k=1:K], sum(vibr_Q[k,n_ibr]  for n_ibr in 1:Nibr if sum(Bibr[n_N,n_ibr])!=0) ), dims=[1,2][min(n_N,2)]);
                else
                        println("Problem: Undefined Nodal Voltage Equation")
                end
        end

        # Load Currents
        global (iL_D,iL_Q) = (Matrix{AffExpr}[],Matrix{AffExpr}[])
        for n_L in 1:Nload
                if sum(Bload[:,n_L])!=0 
                        if (a_p[n_L]==1) && (a_q[n_L]==1)
                                global iL_D = cat(iL_D, @expression(model_traj,[k=1:K],(
                                                        sum(a_p[n_L]*vN_D[k,Bload[:,n_L].==1])*Pload_traj[k,n_L] + a_q[n_L]*sum(vN_Q[k,Bload[:,n_L].==1])*Qload_traj[k,n_L]
                                                        ) ), dims=[1,2][min(n_L,2)]);
                                global iL_Q = cat(iL_Q, @expression(model_traj,[k=1:K],(
                                                        sum(a_p[n_L]*vN_Q[k,Bload[:,n_L].==1])*Pload_traj[k,n_L] - a_q[n_L]*sum(vN_D[k,Bload[:,n_L].==1])*Qload_traj[k,n_L]
                                                        ) ), dims=[1,2][min(n_L,2)]);
                        end
                else
                        println("Problem: Undefined Load current Equation")
                end
        end


        global (iibr_D,iibr_Q) = (Matrix{AffExpr}[],Matrix{AffExpr}[])
        for n_ibr in 1:Nibr
                if sum(Be[Bibr[:,n_ibr].==1,:].!=0)!=0 && sum(Bload[Bibr[:,n_ibr].==1,:])!=0
                        global iibr_D = cat(iibr_D, @expression(model_traj,[k=1:K], sum((Be*ie_D[k,:])[Bibr[:,n_ibr].==1]) + sum((Bload*iL_D[k,:])[Bibr[:,n_ibr].==1]) ), dims=[1,2][min(n_ibr,2)]);
                        global iibr_Q = cat(iibr_Q, @expression(model_traj,[k=1:K], sum((Be*ie_Q[k,:])[Bibr[:,n_ibr].==1]) + sum((Bload*iL_Q[k,:])[Bibr[:,n_ibr].==1]) ), dims=[1,2][min(n_ibr,2)]);
                
                elseif sum(Be[Bibr[:,n_ibr].==1,:].!=0)!=0
                        global iibr_D = cat(iibr_D, @expression(model_traj,[k=1:K], sum((Be*ie_D[k,:])[Bibr[:,n_ibr].==1])                                             ), dims=[1,2][min(n_ibr,2)]);
                        global iibr_Q = cat(iibr_Q, @expression(model_traj,[k=1:K], sum((Be*ie_Q[k,:])[Bibr[:,n_ibr].==1])                                             ), dims=[1,2][min(n_ibr,2)]);
                else
                        println("Problem: Undefined IBR Current Equation")
                end
        end

        # DQ-to-dq
        global iibr_d   = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], iibr_D[k,n_ibr]*cos(-delta_ibr[k,n_ibr])-iibr_Q[k,n_ibr]*sin(-delta_ibr[k,n_ibr]));
        global iibr_q   = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], iibr_D[k,n_ibr]*sin(-delta_ibr[k,n_ibr])+iibr_Q[k,n_ibr]*cos(-delta_ibr[k,n_ibr]));

        #............. IBRs .................

        # Power calculation
        global p_ibr    = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],(vibr_d[k,n_ibr]*iibr_d[k,n_ibr]+vibr_q[k,n_ibr]*iibr_q[k,n_ibr]))
        global q_ibr    = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],(vibr_q[k,n_ibr]*iibr_d[k,n_ibr]-vibr_d[k,n_ibr]*iibr_q[k,n_ibr]))

        # ----------------------------------------------------------------------
        # ------------------------- Control Saturation -------------------------
        # ----------------------------------------------------------------------

        # **************************************** w_dq *************************************
        global w_dq    = Vector{AffExpr}[]
        for n_ibr in 1:Nibr
                # ----- GFM -----
                if y_theta[n_ibr]==0
                        global w_dq    =  cat(  w_dq,    @expression(model_traj,[k=1:K], n1_ratio[n_ibr]*(p_ref[k,n_ibr]-m_p_ibr[k,n_ibr])+w_0),dims=[1,2][min(n_ibr,2)]);
                # ----- GFL -----
                elseif y_theta[n_ibr]==1 
                        global w_dq    =  cat(  w_dq,   @expression(model_traj,[k=1:K],   hlpf_pll[k,n_ibr]+w_0),dims=[1,2][min(n_ibr,2)]);
                end
        end


        # ******************************* v_dq_ref *********************************
        global v_d_ref = Vector{AffExpr}[] # v_q_ref is always zero
        for n_ibr in 1:Nibr
                # ----- GFM -----
                if y_theta[n_ibr]==0
                        global v_d_ref = cat(  v_d_ref, @expression(model_traj,[k=1:K], n2_ratio[n_ibr]*(q_ref[k,n_ibr]-m_q_ibr[k,n_ibr]) + 1 ),dims=[1,2][min(n_ibr,2)]);
                # ----- GFL -----
                elseif y_theta[n_ibr]==1
                        # dummy |V|
                        global v_d_ref  = cat( v_d_ref, @expression(model_traj,[k=1:K], (vibr_d[k,n_ibr]^2 + vibr_q[k,n_ibr]^2)^0.5 ) ,dims=[1,2][min(n_ibr,2)]);
                end
        end

        # ***************************** ix_dq_ref *********************************
        global (ix_d_ref, ix_q_ref) = (Vector{AffExpr}[],Vector{AffExpr}[])
        for n_ibr in 1:Nibr
                # ----- GFM -----
                if y_theta[n_ibr]==0
                        global ix_d_ref = cat(   ix_d_ref, @expression(model_traj,[k=1:K], Kv_P[n_ibr]*(v_d_ref[k,n_ibr]-vibr_d[k,n_ibr])+hv_d[k,n_ibr]+iibr_d[k,n_ibr]
                                                                        -w_dq[k,n_ibr]*C_f[n_ibr]*vibr_q[k,n_ibr]),  dims=[1,2][min(n_ibr,2)]);
                        global ix_q_ref = cat(   ix_q_ref, @expression(model_traj,[k=1:K], Kv_P[n_ibr]*(                -vibr_q[k,n_ibr])+hv_q[k,n_ibr]+iibr_q[k,n_ibr]
                                                                        +w_dq[k,n_ibr]*C_f[n_ibr]*vibr_d[k,n_ibr]),  dims=[1,2][min(n_ibr,2)]);
                # ----- GFL -----
                elseif y_theta[n_ibr]==1
                        global ix_d_ref =  cat(  ix_d_ref, @expression(model_traj,[k=1:K], hpq_d[k,n_ibr] + (+Kpq_P[n_ibr])*(p_ref[k,n_ibr]-m_p_ibr[k,n_ibr]+ (1/n1_ratio[n_ibr])*(w_0-w_dq[k,n_ibr]  )))                                        ,dims=[1,2][min(n_ibr,2)]);
                        global ix_q_ref =  cat(  ix_q_ref, @expression(model_traj,[k=1:K], hpq_q[k,n_ibr] + (-Kpq_P[n_ibr])*(q_ref[k,n_ibr]-m_q_ibr[k,n_ibr]+ (1/n2_ratio[n_ibr])*(1- v_d_ref[k,n_ibr])))                              ,dims=[1,2][min(n_ibr,2)]);
                end
        end

        # *********** Current saturation *****************
        if include_saturation == true
                global Z_ix         = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], ix_ibr_max[n_ibr]/ ((ix_d_ref[k,n_ibr]^2 + ix_q_ref[k,n_ibr]^2 + 1e-8 )^0.5) );
                global alpha_ix     = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], 0.5* ( 1 + Z_ix[k,n_ibr] - (((Z_ix[k,n_ibr]-1)^2)^(0.5)) )  );
                global ix_d_ref_sat = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], alpha_ix[k,n_ibr] * ix_d_ref[k,n_ibr] );
                global ix_q_ref_sat = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], alpha_ix[k,n_ibr] * ix_q_ref[k,n_ibr] );
        else
                global ix_d_ref_sat = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], ix_d_ref[k,n_ibr] );
                global ix_q_ref_sat = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], ix_q_ref[k,n_ibr] );
        end


        # ----------------------------------------------------------
        # ----------------------------------------------------------
        # ----------------------------------------------------------

        # Control variables: Current Controller
        global vx_d_ref = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], Ki_P[n_ibr]*(ix_d_ref_sat[k,n_ibr]-ix_d[k,n_ibr])+hi_d[k,n_ibr]+vibr_d[k,n_ibr]
                                                +R_f[n_ibr]*ix_d[k,n_ibr]
                                                -w_dq[k,n_ibr]*L_f[n_ibr]*ix_q[k,n_ibr]);
        global vx_q_ref = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], Ki_P[n_ibr]*(ix_q_ref_sat[k,n_ibr]-ix_q[k,n_ibr])+hi_q[k,n_ibr]+vibr_q[k,n_ibr]
                                                +R_f[n_ibr]*ix_q[k,n_ibr]
                                                +w_dq[k,n_ibr]*L_f[n_ibr]*ix_d[k,n_ibr]);

        # Control variables: Modulation
        global m_d_ref = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],vx_d_ref[k,n_ibr]/Vdc_n[n_ibr]);
        global m_q_ref = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],vx_q_ref[k,n_ibr]/Vdc_n[n_ibr]);

        # Physical variables
        global iz   = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],m_d_ref[k,n_ibr]*ix_d[k,n_ibr]+m_q_ref[k,n_ibr]*ix_q[k,n_ibr]);
        global vx_d = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],vdc[k,n_ibr]*m_d_ref[k,n_ibr]);
        global vx_q = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],vdc[k,n_ibr]*m_q_ref[k,n_ibr]);

        # DC Power calculation
        global pdc_ibr  = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],(vdc[k,n_ibr]*iz[k,n_ibr]))
        
        #-------------------------------------------------------------------
        #------------------------- Derivatives -----------------------------
        #-------------------------------------------------------------------

        #............Power System............

        # Dynamics Shunt voltages
        global (fvsh_D,fvsh_Q) = (Matrix{AffExpr}[],Matrix{AffExpr}[])
        for n_sh in 1:Nsh # there's no shunt capacitor for IBR's buses, already included in IBR's capacitive filter C_f
                if sum(Bload[Bsh[:,n_sh].==1,:])!=0
                        global fvsh_D = cat(fvsh_D, @expression(model_traj,[k=1:K], (w_b/C_sh[n_sh])*( -sum((Be*ie_D[k, :])[Bsh[:,n_sh].==1]) - sum((Bload*iL_D[k, :])[Bsh[:,n_sh].==1]) )  + w_b*vsh_Q[k,n_sh]), dims=[1,2][min(n_sh,2)])
                        global fvsh_Q = cat(fvsh_Q, @expression(model_traj,[k=1:K], (w_b/C_sh[n_sh])*( -sum((Be*ie_Q[k, :])[Bsh[:,n_sh].==1]) - sum((Bload*iL_Q[k, :])[Bsh[:,n_sh].==1]) )  - w_b*vsh_D[k,n_sh]), dims=[1,2][min(n_sh,2)])
                else
                        global fvsh_D = cat(fvsh_D, @expression(model_traj,[k=1:K], (w_b/C_sh[n_sh])*( -sum((Be*ie_D[k, :])[Bsh[:,n_sh].==1]) )  + w_b*vsh_Q[k,n_sh]), dims=[1,2][min(n_sh,2)])
                        global fvsh_Q = cat(fvsh_Q, @expression(model_traj,[k=1:K], (w_b/C_sh[n_sh])*( -sum((Be*ie_Q[k, :])[Bsh[:,n_sh].==1]) )  - w_b*vsh_D[k,n_sh]), dims=[1,2][min(n_sh,2)])                                          
                end
        end

        # Dynamics Line Currents
        global (fie_D,fie_Q) = (Matrix{AffExpr}[],Matrix{AffExpr}[])
        for n_e in 1:Ne
                if sum(Be[:,n_e].!=0)!=0
                        global fie_D = cat(fie_D, @expression(model_traj,[k=1:K], (w_b/Le[n_e])*(  transpose(Be[:,n_e][Be[:,n_e].!=0])*vN_D[k,Be[:,n_e].!=0] - Re[n_e]*ie_D[k,n_e]  )  + w_b*ie_Q[k,n_e]  ), dims=[1,2][min(n_e,2)])
                        global fie_Q = cat(fie_Q, @expression(model_traj,[k=1:K], (w_b/Le[n_e])*(  transpose(Be[:,n_e][Be[:,n_e].!=0])*vN_Q[k,Be[:,n_e].!=0] - Re[n_e]*ie_Q[k,n_e]  )  - w_b*ie_D[k,n_e]  ), dims=[1,2][min(n_e,2)])
                else 
                        println("Problem: Undefined Deriv. Line current Equation")
                end
        end

        #............. IBRs .................
        # Dynamics Local dq theta
        global ftheta_ibr_dq = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],w_b*w_dq[k,n_ibr]);
        # Dynamics Physical variables: DC-link voltage
        global fvdc = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],(w_b/Cdc[n_ibr])*(idc[k,n_ibr]-iz[k,n_ibr]));
        # Dynamics Physical variables: RL-filter
        global fix_d = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],(w_b/L_f[n_ibr])*(vx_d[k,n_ibr]-vibr_d[k,n_ibr]-R_f[n_ibr]*ix_d[k,n_ibr])
                                                                +w_b*w_dq[k,n_ibr]*ix_q[k,n_ibr]);
        global fix_q = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],(w_b/L_f[n_ibr])*(vx_q[k,n_ibr]-vibr_q[k,n_ibr]-R_f[n_ibr]*ix_q[k,n_ibr])
                                                                -w_b*w_dq[k,n_ibr]*ix_d[k,n_ibr]);
        # Dynamics Physical variables: C-filter
        global fvibr_d = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],(w_b/C_f[n_ibr])*(ix_d[k,n_ibr]-iibr_d[k,n_ibr])
                                                                +w_b*w_dq[k,n_ibr]*vibr_q[k,n_ibr]);
        global fvibr_q = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],(w_b/C_f[n_ibr])*(ix_q[k,n_ibr]-iibr_q[k,n_ibr])
                                                                -w_b*w_dq[k,n_ibr]*vibr_d[k,n_ibr]); 

        # Dynamic control variables
        global (fhpll,fhlpf_pll,fhpq_d,fhpq_q,fhv_d,fhv_q) = (Vector{AffExpr}[],Vector{AffExpr}[],Vector{AffExpr}[],Vector{AffExpr}[],Vector{AffExpr}[],Vector{AffExpr}[])
        
        for n_ibr in 1:Nibr
                # ----- GFM -----   
                if y_theta[n_ibr]==0                      
                        # Dynamics Control variables: Voltage Controller
                        if Kv_I[n_ibr]!=0
                                global fhv_d    = cat(fhv_d,@expression(model_traj,[k=1:K],    Kv_I[n_ibr]*(v_d_ref[k,n_ibr]-vibr_d[k,n_ibr]+ix_d_ref_sat[k,n_ibr]-ix_d_ref[k,n_ibr]))         ,dims=[1,2][min(n_ibr,2)]);
                                global fhv_q    = cat(fhv_q,@expression(model_traj,[k=1:K],    Kv_I[n_ibr]*(-vibr_q[k,n_ibr]+ix_q_ref_sat[k,n_ibr]-ix_q_ref[k,n_ibr]))         ,dims=[1,2][min(n_ibr,2)]);
                        else
                                global fhv_d    = cat(fhv_d,    zeros(K,1)      ,dims=[1,2][min(n_ibr,2)]);
                                global fhv_q    = cat(fhv_q,    zeros(K,1)      ,dims=[1,2][min(n_ibr,2)]);
                         end
                        # dummies
                        global fhpll    = cat(fhpll,            zeros(K,1)      ,dims=[1,2][min(n_ibr,2)]);
                        global fhlpf_pll= cat(fhlpf_pll,        zeros(K,1)      ,dims=[1,2][min(n_ibr,2)]);
                        global fhpq_d   = cat(fhpq_d,           zeros(K,1)      ,dims=[1,2][min(n_ibr,2)]);
                        global fhpq_q   = cat(fhpq_q,           zeros(K,1)      ,dims=[1,2][min(n_ibr,2)]);
                
                # ----- GFL -----
                elseif y_theta[n_ibr]==1           
                        # Dynamics Control variables: PLL model_traj
                        global fhpll     = cat(fhpll,      @expression(model_traj,[k=1:K], Kpll_I[n_ibr]*(-vibr_q[k,n_ibr]))                                                                             ,dims=[1,2][min(n_ibr,2)]);
                        global fhlpf_pll = cat(fhlpf_pll,  @expression(model_traj,[k=1:K], (1/tau_lpf_pll[n_ibr])*(Kpll_P[n_ibr]*(-vibr_q[k,n_ibr])+hpll[k,n_ibr] - hlpf_pll[k,n_ibr])  )                ,dims=[1,2][min(n_ibr,2)]);
                        # Dynamics Control variables: PQ controller
                        global fhpq_d    = cat(fhpq_d, @expression(model_traj,[k=1:K], (+Kpq_I[n_ibr])*(p_ref[k,n_ibr] +(1/n1_ratio[n_ibr])*(w_0-w_dq[k,n_ibr] ) -m_p_ibr[k,n_ibr] + ix_d_ref_sat[k,n_ibr]-ix_d_ref[k,n_ibr]))          ,dims=[1,2][min(n_ibr,2)]);
                        global fhpq_q    = cat(fhpq_q, @expression(model_traj,[k=1:K], (-Kpq_I[n_ibr])*(q_ref[k,n_ibr] +(1/n2_ratio[n_ibr])*(1-v_d_ref[k,n_ibr]) -m_q_ibr[k,n_ibr] + ix_q_ref_sat[k,n_ibr]-ix_q_ref[k,n_ibr]))          ,dims=[1,2][min(n_ibr,2)]);
                        
                        # dummies
                        global fhv_d     = cat(fhv_d, zeros(K,1)                ,dims=[1,2][min(n_ibr,2)]);
                        global fhv_q     = cat(fhv_q, zeros(K,1)                ,dims=[1,2][min(n_ibr,2)]);  
                end
        end


        # Dynamics Control variables: Current Controller
        global fhi_d = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],  Ki_I[n_ibr]*(ix_d_ref_sat[k,n_ibr]-ix_d[k,n_ibr]));
        global fhi_q = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],  Ki_I[n_ibr]*(ix_q_ref_sat[k,n_ibr]-ix_q[k,n_ibr]));
        
        # DC link control
        global fidc = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], (1/tau_dc[n_ibr])*  (Kdc_P[n_ibr]*(Vdc_n[n_ibr]-vdc[k,n_ibr]) + hvdc[k,n_ibr] +
                                                                                        (1/Vdc_n[n_ibr])*(p_ref[k,n_ibr]+pdc_ibr[k,n_ibr]-p_ibr[k,n_ibr])
                                                                                        - idc[k,n_ibr]));
        global fhvdc = @expression(model_traj,[k=1:K,n_ibr=1:Nibr], Kdc_I[n_ibr]*(Vdc_n[n_ibr]-vdc[k,n_ibr]));

        # IBR measurements
        global fm_p_ibr   = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],   w_c[n_ibr]* (p_ibr[k,n_ibr]   - m_p_ibr[k,n_ibr]  )  );
        global fm_q_ibr   = @expression(model_traj,[k=1:K,n_ibr=1:Nibr],   w_c[n_ibr]* (q_ibr[k,n_ibr]   - m_q_ibr[k,n_ibr]  )  );

        #-------------------------- Score variables ----------------------------

        global fscore = @expression(model_traj,[k=1:K],  sum( (p_ref[k,:] - p_ibr[k,:]).^2 + (q_ref[k,:] - q_ibr[k,:]).^2 )  )

        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------

        global Vars   = @expression(model_traj,cat(          vsh_D , vsh_Q , ie_D , ie_Q , ix_d , ix_q , vibr_d , vibr_q , vdc , idc , m_p_ibr , m_q_ibr ,
                                                                hvdc , hv_d , hv_q , hi_d , hi_q , hpll , hlpf_pll , hpq_d , hpq_q , theta_ibr_dq , score , dims=2));

        global Derivs = @expression(model_traj,cat(          fvsh_D, fvsh_Q, fie_D, fie_Q, fix_d, fix_q, fvibr_d, fvibr_q, fvdc, fidc, fm_p_ibr, fm_q_ibr,
                                                                fhvdc, fhv_d, fhv_q, fhi_d, fhi_q, fhpll, fhlpf_pll, fhpq_d, fhpq_q, ftheta_ibr_dq, fscore ,dims=2));
        return(model_traj)
end

global (hpq_ini,hpq_end)=(Nx-3*Nibr,Nx-Nibr-1)            # <- Indicates locations of hpq_d and hpq_q in Vars
global ( ix_ini, ix_end)=(2*Nsh+2*Ne+1,2*Nsh+2*Ne+2*Nibr) # <- Indicates locations of ix_d and ix_q in Vars
#-------------------------------------------------------------------
#--------------------- Recover variables from Vars -----------------
#-------------------------------------------------------------------

function recover_vars(vars_matrix)
        n_1 =  1      ;  n_2=n_1-1+Nsh    ; var_vsh_D          = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nsh    ; var_vsh_Q          = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Ne     ; var_ie_D           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Ne     ; var_ie_Q           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_ix_d           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_ix_q           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_vibr_d         = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_vibr_q         = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_vdc            = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_idc            = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_m_p_ibr        = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_m_q_ibr        = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_hvdc           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_hv_d           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_hv_q           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_hi_d           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_hi_q           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_hpll           = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_hlpf_pll       = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_hpq_d          = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_hpq_q          = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+Nibr   ; var_theta_ibr_dq   = vars_matrix[:,n_1:n_2];
        n_1 =  n_2+1  ;  n_2=n_1-1+1      ; var_score          = vars_matrix[:,n_1:n_2];
        return  (       var_vsh_D,var_vsh_Q,var_ie_D,var_ie_Q,var_ix_d,var_ix_q,var_vibr_d,var_vibr_q,var_vdc,var_idc,var_m_p_ibr,var_m_q_ibr,
                        var_hvdc,var_hv_d,var_hv_q,var_hi_d,var_hi_q,var_hpll,var_hlpf_pll,var_hpq_d,var_hpq_q,var_theta_ibr_dq,var_score)
end

function save_to_matlab(file_name,vars_matrix,y_theta_save,n1_ratio_save,n2_ratio_save,Pload_ini_save,Pload_fin_save,Qload_ini_save,Qload_fin_save,p_ref_ed_save,q_ref_ed_save,p_ref_agc_save) # Recover variables to paste into model to recover Derivs

        (var_vsh_D,var_vsh_Q,var_ie_D,var_ie_Q,var_ix_d,var_ix_q,var_vibr_d,var_vibr_q,var_vdc,var_idc,var_m_p_ibr,var_m_q_ibr,
         var_hvdc,var_hv_d,var_hv_q,var_hi_d,var_hi_q,var_hpll,var_hlpf_pll,var_hpq_d,var_hpq_q,var_theta_ibr_dq,var_score) = recover_vars(value.(vars_matrix));

        global matlab_sim = matopen(file_name, "w")

        # -------- User Settings --------
        write(matlab_sim, "H"             , H                   );
        write(matlab_sim, "ref_upd_time"  , ref_upd_time        );
        write(matlab_sim, "phase"         , phase               );
        write(matlab_sim, "Optim_iter"    , Optim_iter          );
        write(matlab_sim, "ibr_conf_num"  , ibr_conf_num        );

        # -------------- IBR Params -------------------        

        write(matlab_sim, "y_theta"     , y_theta_save          );
        write(matlab_sim, "n1_ratio"    , n1_ratio_save         );
        write(matlab_sim, "n2_ratio"    , n2_ratio_save         );

        # --------- Load and references ---------------
        write(matlab_sim, "Pload_ini"   ,  Pload_ini_save       );
        write(matlab_sim, "Pload_fin"   ,  Pload_fin_save       );
        write(matlab_sim, "Qload_ini"   ,  Qload_ini_save       );
        write(matlab_sim, "Qload_fin"   ,  Qload_fin_save       );
        write(matlab_sim, "p_ref_ed"    ,  p_ref_ed_save        );
        write(matlab_sim, "q_ref_ed"    ,  q_ref_ed_save        );
        write(matlab_sim, "p_ref_agc"   ,  p_ref_agc_save       );
        # ------- Shunt Capacitors ------
        write(matlab_sim, "vsh_D_0", var_vsh_D[1,:])
        write(matlab_sim, "vsh_Q_0", var_vsh_Q[1,:])
        # ------------ Lines ------------
        write(matlab_sim, "ie_D_0", var_ie_D[1,:])
        write(matlab_sim, "ie_Q_0", var_ie_Q[1,:])
        # ------------- IBRs ------------

        # Physical model_traj: Inductor Filter
        write(matlab_sim, "ix_d_0", var_ix_d[1,:])
        write(matlab_sim, "ix_q_0", var_ix_q[1,:])
        # Physical model_traj: Capacitor Filter
        write(matlab_sim, "vibr_d_0", var_vibr_d[1,:])
        write(matlab_sim, "vibr_q_0", var_vibr_q[1,:])
        # Physical model_traj: DC voltage
        write(matlab_sim, "vdc_0", var_vdc[1,:]);
        # DC current input
        write(matlab_sim, "idc_0", var_idc[1,:]);
        # Power measurements
        write(matlab_sim, "m_p_ibr_0", var_m_p_ibr[1,:]);
        write(matlab_sim, "m_q_ibr_0", var_m_q_ibr[1,:]);

        # DC-link Control
        write(matlab_sim, "hvdc_0", var_hvdc[1,:])
        # Control model_traj: Voltage control
        write(matlab_sim, "hv_d_0", var_hv_d[1,:])
        write(matlab_sim, "hv_q_0", var_hv_q[1,:])
        # Control model_traj: Current control
        write(matlab_sim, "hi_d_0", var_hi_d[1,:])
        write(matlab_sim, "hi_q_0", var_hi_q[1,:])
        
        # Control model_traj: PLL
        write(matlab_sim, "hpll_0", var_hpll[1,:])
        write(matlab_sim, "hlpf_pll_0", var_hlpf_pll[1,:])

        # Control model_traj: PQ control
        write(matlab_sim, "hpq_d_0", var_hpq_d[1,:])
        write(matlab_sim, "hpq_q_0", var_hpq_q[1,:])

        # theta IBR
        write(matlab_sim, "theta_ibr_dq_0", var_theta_ibr_dq[1,:])
        
        # Scores
        write(matlab_sim, "score_0", var_score[1,:])

        close(matlab_sim)
end


function save_results(filename, BREAK_TIMES , T_INT, VARS_INT, DERIVS_INT, PLOAD_INI , PLOAD_FIN , QLOAD_INI , QLOAD_FIN , P_REF_ED , Q_REF_ED , P_REF_AGC , 
                                IBR_CONF_NUM, Y_THETA , N1_RATIO , N2_RATIO , Obj_NLP, Time_NLP, Obj_EMT, Time_EMT)
        results = matopen(filename*".mat", "w")
        write(results, "BREAK_TIMES"    , BREAK_TIMES   );
        write(results, "T_INT"          , T_INT         );
        write(results, "VARS_INT"       , VARS_INT      );
        write(results, "DERIVS_INT"     , DERIVS_INT    );
        write(results, "PLOAD_INI"      , PLOAD_INI     );
        write(results, "PLOAD_FIN"      , PLOAD_FIN     );
        write(results, "QLOAD_INI"      , QLOAD_INI     );
        write(results, "QLOAD_FIN"      , QLOAD_FIN     );
        write(results, "P_REF_ED"       , P_REF_ED      );
        write(results, "Q_REF_ED"       , Q_REF_ED      );
        write(results, "P_REF_AGC"      , P_REF_AGC     );
        write(results, "IBR_CONF_NUM"   , IBR_CONF_NUM  );
        write(results, "Y_THETA"        , Y_THETA       );
        write(results, "N1_RATIO"       , N1_RATIO      );
        write(results, "N2_RATIO"       , N2_RATIO      );  
        write(results, "OBJ_NLP"        , Obj_NLP       );
        write(results, "TIME_NLP"       , Time_NLP      );
        write(results, "OBJ_EMT"        , Obj_EMT       );
        write(results, "TIME_EMT"       , Time_EMT      );
        close(results)
        println(filename," - RESULTS SAVED")
        save_to_matlab(filename*"_for_sim.mat",VARS_INT,Y_THETA,N1_RATIO,N2_RATIO,PLOAD_INI,PLOAD_FIN,QLOAD_INI,QLOAD_FIN,P_REF_ED,Q_REF_ED,P_REF_AGC)
        println("For_Sim_"*filename," - RESULTS SAVED")
end

function read_results(filename)
        results = matopen(filename*".mat");
        variables_of_interest=[ "BREAK_TIMES",  "T_INT", "VARS_INT", "DERIVS_INT",
                                "PLOAD_INI", "PLOAD_FIN", "QLOAD_INI", "QLOAD_FIN", "P_REF_ED", "Q_REF_ED", "P_REF_AGC",
                                "IBR_CONF_NUM" , "Y_THETA" , "N1_RATIO","N2_RATIO", 
                                "OBJ_NLP", "TIME_NLP", "OBJ_EMT", "TIME_EMT"]
        ( BREAK_TIMES , T_INT, VARS_INT, DERIVS_INT, PLOAD_INI , PLOAD_FIN , QLOAD_INI , QLOAD_FIN , P_REF_ED , Q_REF_ED , P_REF_AGC , IBR_CONF_NUM , Y_THETA , N1_RATIO , N2_RATIO , OBJ_NLP, TIME_NLP, OBJ_EMT, TIME_EMT) = [read(results,i) for i in variables_of_interest];
        close(results)
        println(filename," - RESULTS READED")
        return (BREAK_TIMES , T_INT, VARS_INT, DERIVS_INT, PLOAD_INI , PLOAD_FIN , QLOAD_INI , QLOAD_FIN , P_REF_ED , Q_REF_ED , P_REF_AGC , IBR_CONF_NUM , Y_THETA , N1_RATIO , N2_RATIO , OBJ_NLP, TIME_NLP, OBJ_EMT, TIME_EMT)
end
