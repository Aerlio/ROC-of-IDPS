include("Julia_0_0_General_Calls.jl") ;

# --------- Generate Warmstart Points (a priori) ---------
global Generate_WS_points      = false ;
if Generate_WS_points include("Julia_0_2_Generate_Warmstarts.jl") ; end

# ---------------- Local Reduction Loop ------------------
global Optim_iter = 1
while Optim_iter <= Max_Scenarios
    include("Julia_1_1_Max.jl")
    if (Optim_iter == Max_Scenarios) break end # stopping criterion
    include("Julia_1_2_Min.jl") 
    global Optim_iter = Optim_iter + 1
end