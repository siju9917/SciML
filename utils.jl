"""
Source_to_function_of_time(source::PeriodicVariableSource)

Takes in a PeriodicVariableSource from PowerSystems and generates functions of time for voltage magnitude and angle
"""
function Source_to_function_of_time(source::PeriodicVariableSource)
    V_bias = get_internal_voltage_bias(source)
    V_freqs = get_internal_voltage_frequencies(source)
    V_coeffs = get_internal_voltage_coefficients(source)
   function V(t)
       val = V_bias
       for (i,ω) in enumerate(V_freqs)
           val += V_coeffs[i][1]* sin.(ω * t)
           val += V_coeffs[i][2]* cos.(ω * t)
       end
       return val
   end
   θ_bias = get_internal_angle_bias(source)
   θ_freqs = get_internal_angle_frequencies(source)
   θ_coeffs = get_internal_angle_coefficients(source)
  function θ(t)
      val = θ_bias
      for (i,ω) in enumerate(θ_freqs)
          val += θ_coeffs[i][1]* sin.(ω * t)
          val += θ_coeffs[i][2]* cos.(ω * t)
      end
       return val
   end 
   return (V, θ)
end

"""

Makes a Float64 Mass Matrix of ones for the ODEProblem. Takes # of differential and algebraic states

"""
function MassMatrix(n_differential::Integer, n_algebraic::Integer)
    n_states = n_differential + n_algebraic
    M = Float64.(zeros(n_states,n_states))
    for i = 1:n_differential 
      M[i,i] = 1.0
    end
    return M
end


function get_parameters(inv::DynamicInverter)
    #pll
    p = Vector{Float64}(undef,23)
    p[1] = get_ω_lp(inv.freq_estimator)
    p[2] = get_kp_pll(inv.freq_estimator)
    p[3] = get_ki_pll(inv.freq_estimator)
    #outer control
    p[4] = PowerSystems.get_Ta(inv.outer_control.active_power)
    p[5] = PowerSystems.get_kd(inv.outer_control.active_power)
    p[6] = PowerSystems.get_kω(inv.outer_control.active_power)
    p[7] = PowerSystems.get_kq(inv.outer_control.reactive_power)
    p[8] = PowerSystems.get_ωf(inv.outer_control.reactive_power)
    #inner control
    p[9] = get_kpv(inv.inner_control)
    p[10] = get_kiv(inv.inner_control)
    p[11] = get_kffv(inv.inner_control)
    p[12] = get_rv(inv.inner_control)
    p[13] = get_lv(inv.inner_control)
    p[14] = get_kpc(inv.inner_control)
    p[15] = get_kic(inv.inner_control)
    p[16] = get_kffi(inv.inner_control)
    p[17] = get_ωad(inv.inner_control)
    p[18] = get_kad(inv.inner_control)
    #lcl
    p[19] = get_lf(inv.filter)
    p[20] = get_rf(inv.filter)
    p[21] = get_cf(inv.filter)
    p[22] = get_lg(inv.filter)
    p[23] = get_rg(inv.filter)
    return p
end

