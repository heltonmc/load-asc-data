function refl_DT(t,β,ρ)

	nmed = 1.35
	ndet = 1.45
	μa = β[1]
    μsp = β[2]
    n = 1.35
    c = 29.9792458 # Speed of light cm/ns
    v = c/n # Speed of light in medium
	D = 1/(3*(μa + μsp)) # Diffusion Coefficient
	z0 = 1/μsp

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

	ze = 2*Afac*D
	zp = z0 + 2*ze

	Rt1 = @. -exp(-(ρ^2)/(4*D*v*t)-μa*v*t)/(2*t^(5/2)*(4*π*D*v)^3/2)
	Rt2 = @. (z0*exp(-(z0^2)/(4*D*v*t)) - zp*exp(-(zp^2)/(4*D*v*t)))
	Rt = @. Rt1*Rt2
    #Rt = @. (v/(4π*D*v*t)^1.5)*exp(-(ρ^2)/(4D*v*t)-μa*v*t)

    replace!(Rt, NaN => 0)
    m = findmax(Rt)
    Rt = Rt./m[1]
end



function refl_DT1(t,β,ρ)

	nmed = 1.35
	ndet = 1.45
	μa = β[1]
	μsp = β[2]
	#ρ = β[3]


	if (nmed == ndet)
		Afac = 1
	elseif nmed > ndet
		costhetac = sqrt(1 - (ndet/nmed).^2)
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		Afac = (2/(1-R0) -1 + abs(costhetac.^3))./(1-abs(costhetac.^2))
	else nmed < ndet
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		Afac = 2/(1-R0) - 1
	end

	z0 = 1/(μa+μsp)
	D = 1/(3*(μa + μsp))
	zb = 2*Afac*D

	v = 29.9792345/nmed # speed of light in medium

	Rt1 = @. v*exp(-μa*v*t)
	Rt1 = @. Rt1/((4*π*D*v*t)^1.5)
	Rt1 = @. Rt1*(exp(-(z0^2 + ρ^2)/(4*D*v*t)) - exp(-((2*zb + z0)^2 + ρ^2)/(4*D*v*t)))

	Rt2 = @. 3*exp(-μa*v*t)/(2*((4*π*D*v)^1.5)*(t^2.5))
	Rt2 = @. Rt2*(z0*exp(-((z0^2 + ρ^2)/(4*D*v*t))) + (2*zb + z0)*exp(-((2*zb+z0)^2 + ρ^2)/(4*D*v*t)))

	Rt = @. abs(Rt1)+abs(Rt2)
	replace!(Rt, NaN => 0)
	m = findmax(Rt)
	Rt = Rt./m[1]
end
