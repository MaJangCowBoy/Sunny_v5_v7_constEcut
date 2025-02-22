function define_system(CoTa3S6,main_keyword,J1,B1,J2,J3,Jc1mat,Jc2,Kz,dims)
  
  sys = System(CoTa3S6, [1 => Moment(s=3/2, g=2)], :dipole);

  if main_keyword == "3Q"
    set_pair_coupling!(sys, (Si, Sj) -> J1*(Si'*Sj) + B1*(Si'*Sj)^2, Bond(1, 1, [1, 0, 0]));
  elseif main_keyword == "1Q"
    set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]));
  else  error("Invalid keyword.");  end

  set_exchange!(sys, J2, Bond(1,1,[1,-1, 0]));
  set_exchange!(sys, J3, Bond(1,1,[2, 0, 0]));
  if main_keyword == "1Q"
    Jc1mat = Jc1mat * [1 0 0; 0 1 0; 0 0 1] + Jc1mat * dmvec([0,0,0.001]);
  end
  set_exchange!(sys, Jc1mat, Bond(1,2,[0, 0, 0]));
  set_exchange!(sys, Jc2, Bond(1,2,[1, 1, 0]));
  set_onsite_coupling!(sys, S -> Kz*S[3]^2, 1);
  sys = repeat_periodically(sys, dims);

  return sys;
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
  
  for _ in 1:50  minimize_energy!(sys; maxiters = 3000);  end

  return sys;
end

function define_qgrid(cryst,axis1,axis2,N1,N2; topbot1 = [-0.5, +0.5], topbot2 = [-0.5, +0.5])

  avctr = [1, 0, 0];  bvctr = [0.5,√3/2,0.0];

  range1 = [];  for x in range(topbot1[1],topbot1[2],N1)  push!(range1, x)  end;
  range1 = Float64.(range1);  norm1 = norm(axis1[1]*avctr + axis1[2]*bvctr);

  range2 = [];  for y in range(topbot2[1],topbot2[2],N2)  push!(range2, y)  end;
  range2 = Float64.(range2);  norm2 = norm(axis2[1]*avctr + axis2[2]*bvctr);

  qgrid = q_space_grid(cryst, axis1, range1, axis2, range2; offset = [0,0,1]);
  return qgrid, range1, range2, norm1, norm2;
end

function define_qline(cryst,axis1,axis2,N1,N2; topbot1 = [-1/6, +1/6], topbot2 = [-1/6, +1/6])

  avctr = [1, 0, 0];  bvctr = [0.5,√3/2,0.0];

  range1 = [];  for x in range(topbot1[1],topbot1[2],N1)  push!(range1, x)  end;
  range1 = Float64.(range1);  norm1 = norm(axis1[1]*avctr + axis1[2]*bvctr);

  range2 = [];  for y in range(topbot2[1],topbot2[2],N2)  push!(range2, y)  end;
  range2 = Float64.(range2);  norm2 = norm(axis2[1]*avctr + axis2[2]*bvctr);

  qline1 = [[1/3, 0, 1] + topbot1[1] * axis1, [1/3, 0, 1] + topbot1[2] * axis1];
  qline2 = [[1/3, 0, 1] + topbot2[1] * axis2, [1/3, 0, 1] + topbot2[2] * axis2];
  qpath1 = q_space_path(cryst, qline1, N1);
  qpath2 = q_space_path(cryst, qline2, N2);
  return qpath1, qpath2, range1, range2, norm1, norm2;
end