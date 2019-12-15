function refl_DT(t,β)
    μa = β[1]
    μsp = β[2]
    n = 1.35
    ρ = 1
    c = 29.9792458 # Speed of light cm/ns
    v = c/n # Speed of light in medium
    D = 1/(3μsp) # Diffusion Coefficient
    Rt = @. (v/(4π*D*v*t)^1.5)*exp(-(ρ^2)/(4D*v*t)-μa*v*t)
    replace!(Rt, NaN => 0)
    m = findmax(Rt)
    Rt = Rt./m[1]
    #Rt = Rt[m[2]-5:end]
    #del = length(t)-length(Rt)
    #Rt = [Rt; zeros(del)]
    #Rt = Rt[end:-1:1,end:-1:1]

end
