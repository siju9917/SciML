using Pkg
Pkg.activate(".")
using DifferentialEquations
using PowerSystems
using PowerSimulationsDynamics
using Plots
include("utils.jl")
include("vsm.jl")

#SET SOLVER PARAMTERS 
solver = Rodas4()      
abstol = 1e-6
reltol = 1e-3
tfault =  0.1
tspan = (0.0, 2.0)
steps = 150
tsteps = tspan[1]:((tspan[2]-tspan[1])/steps):tspan[2]

#READ IN A PSID SYSTEM 
sys  = System("sys_init.json")

#COLLECT PARAMETERS FROM VARIOUS COMPONENTS TO USE WHEN SOLVING WITHOUT PSID
pvs = collect(get_components(PeriodicVariableSource, sys))[1]
transformer = collect(get_components(Transformer2W,sys))[1]
p_fixed =  [get_x(transformer) + get_X_th(pvs), get_r(transformer)+ get_R_th(pvs)]
device = get_component(DynamicInverter, sys, "gen1")
sim= Simulation!(
    MassMatrixModel,
    sys,
    pwd(),
    (0.0,2.0),
)
x₀_dict = read_initial_conditions(sim)["gen1"]
x₀ = [value for (key,value) in x₀_dict]
wrappers = sim.inputs.dynamic_injectors
ref_dict = get_setpoints(sim)["gen1"]
p_inv = get_parameters(device)
refs = [ref_dict["V_ref"], ref_dict["ω_ref"],ref_dict["P_ref"], ref_dict["Q_ref"],]
p_ode = vcat(p_inv, refs, p_fixed)


#SIMULATE THE PSID SYSTEM 
execute!(sim,
    solver,
    abstol = abstol,
    reltol = reltol,
    reset_simulation=false, saveat=tsteps);
sol_psid = read_results(sim)
ir = get_state_series(sol_psid, ("gen1",:ir_filter))
p1 = plot(ir, label="vsm real output current (psid)")

#SIMULATE THE SYSTEM THAT IS WRITTEN BY HAND 
Vm, Vθ = Source_to_function_of_time(pvs)
M = MassMatrix(19,0)
vsm_func = ODEFunction(vsm, mass_matrix = M)
vsm_prob = ODEProblem(vsm_func, x₀, tspan, p_ode)
sol = solve(vsm_prob, solver, abstol=abstol, reltol=reltol,  saveat=tsteps )
p2 = plot(sol.t,sol[5,:],"vsm real output current (scratch)")

plot(p1,p2,layout=(2,1))


