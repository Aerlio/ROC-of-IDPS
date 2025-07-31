include("Julia_0_0_General_Calls.jl") ;

# ------------ Min-Max Loop -----------
global Optim_iter = 1
while Optim_iter <= 7
    include("Julia_1_1_Max.jl")
    if (Optim_iter != 7) include("Julia_1_2_Min.jl") end
    global Optim_iter = Optim_iter + 1
end
