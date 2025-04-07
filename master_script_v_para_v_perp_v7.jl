using Sunny, HDF5, ProgressBars, CairoMakie
using LinearAlgebra, Statistics, Rotations, Printf
include("function_bundle_v7.jl");  include("param_v7.jl");

#? number of threads
npar = Threads.nthreads();
#? number of threads

#? Data loading part ?#
data = Matrix{Float64}(undef, 0, 6);
# f = open("jlist.dat", "r");
f = open("LT_minimize_1p3610_trim.dat", "r");
for line in eachline(f)  tok = parse.(Float64,split(line));  global data = [data; tok'];  end
close(f);  #* data stored like, j2 - jc1 - jc2 - Q1 - Q2 - Q3.
#? Data loading part ?#

#? Basic parameters ?#
formfactors = [1 => FormFactor("Co2")];
cryst = Crystal("CoTaS.cif",symprec=1e-3);  CoTa3S6 = subcrystal(cryst, "Co");
J1 = 1.311;  Kz = -0.000;  b1 = [0.00, 0.01, 0.02, 0.04, 0.06];  B1 = J1 .* b1;
#? Basic parameters ?#

Threads.@threads for i in 1:npar

  steps = Int(floor(size(data,1)/npar));
  idxStt = 1 + (i-1) * steps;  idxEnd = i * steps;
  idxEnd = (i == npar) ? size(data,1) : idxEnd;
  
  println("I and worker id: ", i, " ", Threads.threadid());
  println("I will take on data from ", idxStt, " to ", idxEnd);

  for j in idxStt:idxEnd

    #? file existence check ?#
    h5name = @sprintf("data_%.4d.h5",j);
    if isfile(h5name)  println("File exists: ", h5name);  continue;  end
    #? file existence check ?#

    #? parameter allocation ?#
    j2 = data[j,1];  jc1 = data[j,2];  jc2 = data[j,3];
    F = jc1 + jc2;  G = jc1 - jc2 * 0.5;
    j3 = 1/2 * (1 - (F*G - jc2*F/2 + 2*jc2*G )/√(F*F-2*F*G+4*G*G));
    J2 = j2 * J1;  J3 = j3 * J1;  Jc1 = jc1 * J1;  Jc2 = jc2 * J1;
    regularization = 1e-3;
    #? parameter allocation ?#

    epsilon = 0.05;
    qs_para = [[1/3, 0, 1], [1/3, 0, 1] + epsilon * [+1, 0, 0]];
    path_para = q_space_path(cryst, qs_para, 5);
    qs_perp = [[1/3, 0, 1], [1/3, 0, 1] + epsilon/sqrt(3) * [-1,+2, 0]];
    path_perp = q_space_path(cryst, qs_perp, 5);

    ratio1Q = NaN;  ratio3Qs = NaN * zeros(length(b1));

    #? 1Q
    try
      sys1Q = define_system(CoTa3S6,"1Q",J1,B1,J2,J3,Jc1,Jc2,Kz,(3,1,1));
      sys1Q = system_initialize(sys1Q, "1Q_1", J1);
      measure = ssf_perp(sys1Q; formfactors);  
      swt1Q = SpinWaveTheory(sys1Q; measure, regularization);
    
      res1Q_para = dispersion(swt1Q, path_para);  
      res1Q_perp = dispersion(swt1Q, path_perp);
      ratio1Q = res1Q_para[end,end] / res1Q_perp[end,end];
    catch
      println("err in 1Q LSWT for 2D,", " J2 = ", j2, " Jc2 = ", jc2);
    end
    #? 3Qs
    try
      for it = axes(b1,1)
        sys3Q = define_system(CoTa3S6,"3Q",J1,B1[it],J2,J3,Jc1,Jc2,Kz,(3,3,1));
        sys3Q = system_initialize(sys3Q, "3Q", J1);
        measure = ssf_perp(sys3Q; formfactors);
        swt3Q = SpinWaveTheory(sys3Q; measure, regularization);
        
        res3Q_para = dispersion(swt3Q, path_para);
        res3Q_perp = dispersion(swt3Q, path_perp);
        ratio3Qs[it] = res3Q_para[end,end] / res3Q_perp[end,end];
      end
    catch
      println("err in 3Q LSWT for 2D,", " J2 = ", j2, " Jc2 = ", jc2);
    end

    #? data saving part ?#
    fid = h5open(h5name,"w");
    write(fid,"j2",j2);  write(fid,"jc1",jc1);  write(fid,"jc2",jc2);  write(fid,"b1",b1);
    write(fid,"ratio1Q",ratio1Q);  write(fid,"ratio3Qs",ratio3Qs);
    close(fid);
  end
end
