function CoTaS_5var(dim, J1, j2, j3, jc1, jc2; b1 = 0.06, rng = MersenneTwister(1234), Kz = 0.001)

  CoTaS = Crystal("CoTaS.cif"; symprec = 1e-3);
  cryst = subcrystal(CoTaS,"Co");  S = 1.5;
  randomN = Int64(round(100*rand(rng)));
  sys = System(cryst, dim, [SpinInfo(1, S = S, g=2)], :dipole; seed=randomN); 
  # option --> :dipole, :heisenberg, :ising

  B1  = b1 * J1;
  J2  = j2 * J1;
  J3  = j3 * J1;
  Jc1 = jc1 * J1;
  Jc2 = jc2 * J1;
  
  if abs(B1) < 10 * eps()  set_exchange!(sys, J1,  Bond(1, 1, [1, 0, 0]));
  else  set_pair_coupling!(sys, (Si,Sj) -> J1*(Si'*Sj) + B1*(Si'*Sj)^2, Bond(1, 1, [1, 0, 0]));
  end
  set_exchange!(sys, Jc1, Bond(1, 2, [0, 0, 0]));
  set_exchange!(sys, Jc2, Bond(1, 2, [1, 1, 0]));
  set_exchange!(sys, J2,  Bond(1, 1, [1, 2, 0]));
  set_exchange!(sys, J3,  Bond(1, 1, [2, 0, 0]));
  set_onsite_coupling!(sys, S -> Kz*S[3]^2, 1);
  
  return sys, cryst;

end

function system_initialize(sys::System, keyword::String, J1::Float64)
        
  if keyword == "1Q_1"
    k = [ 1/3, 0, 0];  ψ = π/6;  ϕ = 120 * π/180;  S0 = [0, 0, 1.5];
    for x in axes(sys.dipoles,1), y in axes(sys.dipoles,2), z in axes(sys.dipoles,3), b in axes(sys.dipoles,4)
      θ = 2π * dot(k,[x,y,z]) + ϕ * Float64(isodd(b));
      sys.dipoles[x,y,z,b] = RotZ(ψ) * RotX(θ) * S0;
    end;
  elseif keyword == "1Q_2"
    k = [ 0,1/3, 0];  ψ = π/2;  ϕ = -120 * π/180;  S0 = [0, 0, 1.5];
    for x in axes(sys.dipoles,1), y in axes(sys.dipoles,2), z in axes(sys.dipoles,3), b in axes(sys.dipoles,4)
      θ = 2π * dot(k,[x,y,z]) + ϕ * Float64(isodd(b));
      sys.dipoles[x,y,z,b] = RotZ(ψ) * RotX(θ) * S0;
    end;
  elseif keyword == "1Q_3"
    k = [-1/3,1/3, 0];  ψ = 5π/6;  ϕ = -120 * π/180;  S0 = [0, 0, 1.5];
    for x in axes(sys.dipoles,1), y in axes(sys.dipoles,2), z in axes(sys.dipoles,3), b in axes(sys.dipoles,4)
      θ = 2π * dot(k,[x,y,z]) + ϕ * Float64(isodd(b));
      sys.dipoles[x,y,z,b] = RotZ(ψ) * RotX(θ) * S0;
    end;
  elseif keyword == "3Q"
    S1 = RotZ(π/3) * RotX(2π/3) * [0, 0, 1.5];
    S2 = RotZ(π/3) * RotX(4π/3) * [0, 0, 1.5];
    S3 = RotZ(π/3) * RotX(6π/3) * [0, 0, 1.5];
    for z in axes(sys.dipoles,3)
      sys.dipoles[1,3,z,1] = S2;  sys.dipoles[2,3,z,1] = S3;  sys.dipoles[3,3,z,1] = S3;
      sys.dipoles[1,2,z,1] = S2;  sys.dipoles[2,2,z,1] = S1;  sys.dipoles[3,2,z,1] = S2;
      sys.dipoles[1,1,z,1] = S1;  sys.dipoles[2,1,z,1] = S1;  sys.dipoles[3,1,z,1] = S3;
    
      sys.dipoles[1,3,z,2] = -S2;  sys.dipoles[2,3,z,2] = -S3;  sys.dipoles[3,3,z,2] = -S2;
      sys.dipoles[1,2,z,2] = -S1;  sys.dipoles[2,2,z,2] = -S1;  sys.dipoles[3,2,z,2] = -S2;
      sys.dipoles[1,1,z,2] = -S1;  sys.dipoles[2,1,z,2] = -S3;  sys.dipoles[3,1,z,2] = -S3;
    end
  else
    error("Invalid keyword.")
  end
      
  dt = 0.1/(J1 * 1.5 * 1.5);  damping = 0.1;  
  langevin = Langevin(dt; damping, kT = 0.0);
  langevin.kT = 0.1 * meV_per_K;  for _ in 1:100000  step!(sys, langevin)  end
  langevin.kT = 0.0;              for _ in 1:100000  step!(sys, langevin)  end
  
  for _ in 1:20  minimize_energy!(sys; maxiters = 3000);  end

  return sys;
end

function broaden2D(data2D, xgrd, ygrd, stdx, stdy)
  data2Dbr = zeros(Float64, size(data2D));

  for (i,x) in enumerate(xgrd), (j,y) in enumerate(ygrd)
    _, id = findmin(abs.(xgrd .- (x - 2 * stdx)));  xstt = max(1, id);  
    _, id = findmin(abs.(xgrd .- (y - 2 * stdy)));  ystt = max(1, id);
    _, id = findmin(abs.(xgrd .- (x + 2 * stdx)));  xend = min(length(xgrd), id);
    _, id = findmin(abs.(xgrd .- (y + 2 * stdy)));  yend = min(length(ygrd), id);
    for i0 in xstt:xend, j0 in ystt:yend
      x0 = xgrd[i0];  dx = x - x0; ex = dx / stdx;
      y0 = ygrd[j0];  dy = y - y0; ey = dy / stdy;
      data2Dbr[i,j] += data2D[i0,j0] * exp(-0.5 * (ex^2 + ey^2));
    end
  end
  return data2Dbr;
end
