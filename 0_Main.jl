include("Julia_0_0_General_Calls.jl") ;

# ------------ Min-Max Loop -----------
global Optim_iter = 5
while Optim_iter <= 10
    if (Optim_iter != 5) include("Julia_1_1_Max.jl") end
    include("Julia_1_2_Min.jl")
    include("Julia_1_5_Plot_score.jl")
    global Optim_iter = Optim_iter + 1
end