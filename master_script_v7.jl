using Sunny, HDF5, ProgressBars, CairoMakie
using LinearAlgebra, Statistics, Rotations, Printf
include("function_bundle_v7.jl");

#? number of threads
npar = Threads.nthreads();
#? number of threads

#? Data loading part ?#
data = Matrix{Float64}(undef, 0, 6);
f = open("LT_minimize_1p3610_trim.dat", "r");
for line in eachline(f)  
  tok = parse.(Float64,split(line));  global data = [data; tok'];
end
close(f);  #* data stored like, j2 - jc1 - jc2 - Q1 - Q2 - Q3.
#? Data loading part ?#

#? Basic parameters ?#
kernel = gaussian(fwhm=0.2);  formfactors = [1 => FormFactor("Co2")];
cryst = Crystal("CoTaS.cif",symprec=1e-3);  CoTa3S6 = subcrystal(cryst, "Co");
J1 = 1.311;  Kz = -0.001;
b1 = [0.00, 0.03, 0.06];  B1 = J1 .* b1;
#? Basic parameters ?#

#? define q-points for LSWT ?#
axis1 = [ 1.0, 0.0, 0.0];  N1 = 240;  axis2 = [-0.5, 1.0, 0.0];  N2 = 240;
qgrid, range1, range2, norm1, norm2 = define_qgrid(cryst,axis1,axis2,N1,N2);
qpath1, qpath2, range1, range2, norm1, norm2 = define_qline(cryst,axis1,axis2,N1,N2);
#? define q-points for LSWT ?#


Threads.@threads for i in 1:npar

  steps = Int(floor(size(data,1)/npar));
  idxStt = 1 + (i-1) * steps;
  idxEnd = i * steps;
  idxEnd = (i == npar) ? size(data,1) : idxEnd;

  println("I and worker id: ", i, " ", Threads.threadid());
  println("I will take on data from ", idxStt, " to ", idxEnd);

  for j in idxStt:idxEnd

    #? parameter allocation ?#
    j2 = data[j,1];  jc1 = data[j,2];  jc2 = data[j,3];
    F = jc1 + jc2;  G = jc1 - jc2 * 0.5;
    j3 = 1/2 * (1 - (F*G - jc2*F/2 + 2*jc2*G )/√(F*F-2*F*G+4*G*G));
    J2 = j2 * J1;  J3 = j3 * J1;  Jc1 = jc1 * J1;  Jc2 = jc2 * J1;

    #? file existence check ?#
    h5name = @sprintf("data_%.4d.h5",j);
    if isfile(h5name)  println("File exists: ", h5name);  continue;  end
    #? file existence check ?#

    energies = [0.5 * (J1+J2), 0.6 * (J1+J2), 0.7 * (J1+J2), 0.8 * (J1+J2), 0.9 * (J1+J2)]; # meV
    regularization = 5e-3;
    #? parameter allocation ?#

    #? define LSWT part ?#
    #* 3Q LSWT part.
    sys_3Q = [];  swt_3Q = [];
    for Bi in B1
      sysi = define_system(CoTa3S6,"3Q",J1,Bi,J2,J3,Jc1,Jc2,Kz,(3,3,1));
      keyword = "3Q";  sysi = system_initialize(sysi, keyword, J1);
      measure = ssf_perp(sysi; formfactors);  swti = SpinWaveTheory(sysi; measure, regularization);
      swt_3Q = [swt_3Q; swti];  sys_3Q = [sys_3Q; sysi];
    end
    #* 1Q LSWT part.
    sys_1Q = [];  swt_1Q = [];
    for key in ["1Q_1", "1Q_2", "1Q_3"]
      sysi = define_system(CoTa3S6,"1Q",J1,B1,J2,J3,Jc1,Jc2,Kz,(3,3,1));
      sysi = system_initialize(sysi, key, J1);
      measure = ssf_perp(sysi; formfactors);  swti = SpinWaveTheory(sysi; measure);
      swt_1Q = [swt_1Q; swti];  sys_1Q = [sys_1Q; sysi];
    end
    #? define LSWT part ?#

    #? 2D: calculate intensities and store the data ?#
    data2D_3Q = NaN * zeros(Float64,length(B1),length(energies),N1,N2);  
    data2D_1Q = NaN * zeros(Float64,length(energies),N1,N2);

    for (k,swt) in enumerate(swt_3Q)
      try
        res = intensities(swt, qgrid; energies, kernel);  
        data2D_3Q[k,:,:,:] = res.data[:,:,:];
        # for (id,en) in enumerate(energies)
        #   res = intensities_static(swt, qgrid; bounds = (en-0.1, en+0.1));
        #   data2D_3Q[k,id,:,:] = res.data[:,:];
        # end
      catch
        println("err in 3Q LSWT for 2D, B1 = ", B1[k], " J2 = ", j2, " Jc2 = ", jc2);
      end
    end

    for (k,swt) in enumerate(swt_1Q)
      try
        res = intensities(swt, qgrid; energies, kernel);  
        if k == 1  data2D_1Q[:,:,:]  = res.data[:,:,:];
        else       data2D_1Q[:,:,:] += res.data[:,:,:];  end
        # for (id,en) in enumerate(energies)
        #   res = intensities_static(swt, qgrid; bounds = (en-0.1, en+0.1));  
        #   if k == 1  data2D_1Q[id,:,:]  = res.data[:,:];
        #   else       data2D_1Q[id,:,:] += res.data[:,:];  end
        # end
      catch
        println("err in 1Q LSWT for 2D,", " J2 = ", j2, " Jc2 = ", jc2);
      end
    end
    #? 2D: calculate intensities and store the data ?#

    #? 1D: calculate intensities and store the data ?#
    data1Da_3Q = NaN * zeros(Float64,length(B1),length(energies),N1);  
    data1Da_1Q = NaN * zeros(Float64,length(energies),N1);
    data1Db_3Q = NaN * zeros(Float64,length(B1),length(energies),N2);  
    data1Db_1Q = NaN * zeros(Float64,length(energies),N2);

    for (k,swt) in enumerate(swt_3Q)
      try
        res1 = intensities(swt, qpath1; energies, kernel);  
        data1Da_3Q[k,:,:] = res1.data[:,:];
        res2 = intensities(swt, qpath2; energies, kernel);  
        data1Db_3Q[k,:,:] = res2.data[:,:];
      catch
        println("err in 3Q LSWT for 1D, B1 = ", B1[k], " J2 = ", j2, " Jc2 = ", jc2);
      end
    end
    for (k,swt) in enumerate(swt_1Q)
      try
        res1 = intensities(swt, qpath1; energies, kernel);  
        if k == 1  data1Da_1Q[:,:]  = res1.data[:,:];
        else       data1Da_1Q[:,:] += res1.data[:,:];  end
        res2 = intensities(swt, qpath2; energies, kernel);  
        if k == 1  data1Db_1Q[:,:]  = res2.data[:,:];
        else       data1Db_1Q[:,:] += res2.data[:,:];  end
      catch
        println("err in 1Q LSWT for 1D,", " J2 = ", j2, " Jc2 = ", jc2);
      end
    end
    #? 1D: calculate intensities and store the data ?#
    
    #? data saving part ?#
    fid = h5open(h5name,"w");
    write(fid,"j2",j2);  write(fid,"jc1",jc1);  write(fid,"jc2",jc2);
    write(fid,"range1",range1);  write(fid,"norm1",norm1);
    write(fid,"range2",range2);  write(fid,"norm2",norm2);
    write(fid,"data1Da_3Q", data1Da_3Q);  write(fid,"data1Db_3Q", data1Db_3Q);
    write(fid,"data1Da_1Q", data1Da_1Q);  write(fid,"data1Db_1Q", data1Db_1Q);
    close(fid);
    #? data saving part ?#
  end
end
