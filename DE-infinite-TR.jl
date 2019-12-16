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

function refl_DT1(t,β)

nmed = 1.35
ndet = 1.45
μa = β[1]
μsp = β[2]
ρ = 1


if (nmed == ndet)
	Afac = 1
elseif nmed > ndet
	costhetac = sqrt(1 - (ndet/nmed).^2)
	R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
	Afac = (2/(1-R0) -1 + abs(costhetac.^3))./(1-abs(costhetac.^2))
else #nmed < ndet
	R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
	Afac = 2/(1-R0) - 1
end

z0 = 1/(μa+μsp)
D = 1/(3*(μa + μsp))
zb = 2*Afac*D_med

v = 29.9792345/nmed # speed of light in medium

Rt1 = @. v*exp(-μa*v*t)
Rt1 = @. Rt1/((4*pi*D*v*t)^1.5)
Rt1 = @. Rt1*(exp(-(z0^2 + ρ^2)/(4*D*v*t)) - exp(-((2*zb + z0)^2 + ρ^2)/(4*D*v*t)))

Rt2 = @. 3*exp(-μa*v*t)/(2*((4*pi*D*v)^1.5)*(t^2.5))
Rt2 = @. Rt2*(z0*exp(-((z0^2 + ρ^2)/(4*D*v*t))) + (2*zb + z0)*exp(-((2*zb+z0)^2 + ρ^2)/(4*D*v*t)))

Rt = @. abs(Rt1)+abs(Rt2)
replace!(Rt, NaN => 0)
m = findmax(Rt)
Rt = Rt./m[1]
end
