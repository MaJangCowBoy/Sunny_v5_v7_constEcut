using Dates, JLD2, Statistics, Printf, ProgressBars, Rotations
using Sunny, LinearAlgebra, Random, MAT, Interpolations
using JLD2, HDF5, Printf, CairoMakie
include("function_bundle_v5.jl");

# =======================================
#? number of threads and basic params
# =======================================
npar = Threads.nthreads();

seed = Dates.format(now(),"yyyy-mm-dd HH:MM:SS");  seed = parse(Int64,seed[end-1:end]);
rng = MersenneTwister(seed);

TSet = 3;  kT = TSet * Sunny.meV_per_K;
formfactors = [FormFactor("Co2"; g_lande = 2)];

J1  = 1.3111; 
j2  = parse(Float64,ARGS[1]);
j3  = NaN;
jc1 = 1.3610;
jc2 = parse(Float64,ARGS[2]);

vA = jc1 + 1.0 * jc2;  vB = jc1 - 0.5 * jc2;  vC = jc2;
j3 = 0.5 * (1 - (vA*vB-0.5*vA*vC+2*vB*vC)/sqrt(vA*vA-2*vA*vB+4*vB*vB) );
  # j3 is calculated from jc1 and jc2, which are input parameters.

# =======================================
#? define 1Q/3Q system and calculate correlation for the future purpose
# =======================================
#* 3Q / 1Q
dim_small = (3,3,1);  
sys_small_3Q, cryst = CoTaS_5var(dim_small, J1, j2, j3, jc1, jc2; rng);
sys_small_1Q, cryst = CoTaS_5var(dim_small, J1, j2, j3, jc1, jc2; rng, b1 = 0.00);
sys_small_3Q = system_initialize(sys_small_3Q, "3Q", J1);    print_wrapped_intensities(sys_small_3Q);
sys_small_1Q = system_initialize(sys_small_1Q, "1Q_1", J1);  print_wrapped_intensities(sys_small_1Q);

dim = (30, 30, 5);  sys3Qs = [];  sys1Qs = [];

Threads.@threads for id in 1:npar
  tmp3Q, _ = CoTaS_5var(dim, J1, j2, j3, jc1, jc2; rng = MersenneTwister(id));
  tmp1Q, _ = CoTaS_5var(dim, J1, j2, j3, jc1, jc2; rng = MersenneTwister(id), b1 = 0.00);
  for x in axes(tmp3Q.dipoles,1), y in axes(tmp3Q.dipoles,2), 
        z in axes(tmp3Q.dipoles,3), b in axes(tmp3Q.dipoles,4)
    idx = mod1(x,3);  idy = mod1(y,3);
    tmp3Q.dipoles[x,y,z,b] = +1 * sys_small_3Q.dipoles[idx,idy,1,b];
    tmp1Q.dipoles[x,y,z,b] = +1 * sys_small_1Q.dipoles[idx,idy,1,b];
  end;  minimize_energy!(tmp3Q);  minimize_energy!(tmp1Q);
  push!(sys3Qs, tmp3Q);  push!(sys1Qs, tmp1Q);
end

# =======================================
#? calculate the dynamical correlation
# =======================================

dt = 0.02;  ωmax = 10.0;  nω = 501;
sc3Qs = [dynamical_correlations(sys3Qs[id]; dt, nω, ωmax) for id in 1:npar];
sc1Qs = [dynamical_correlations(sys1Qs[id]; dt, nω, ωmax) for id in 1:npar];

damping = 0.1;  ntherm = 4000;
langevin = Langevin(dt; kT, damping);  
Threads.@threads for _ in 1:npar
  for _ in 1:ntherm  step!(sys3Qs[id], langevin);  step!(sys1Qs[id], langevin);  end
  add_sample!(sc3Qs[id], sys3Qs[id]);  add_sample!(sc1Qs[id], sys1Qs[id]);
end

sc3Q = merge_correlations(sc3Qs);  sc1Q = merge_correlations(sc1Qs);
formula = intensity_formula(sc3Q, :perp; formfactors, kT);  Elist = available_energies(sc3Q);

# =======================================
#? calculate the 2D cut
# =======================================

shift  = [ 0.0, 0.0, 1.0];

N = 151;  xGrd = [ x for x in range(-0.5,0.5,N)];  yGrd = [ y for y in range(-0.5,0.5,N)];

basis1s = [ [ 1.0, 0.0, 0.0], [ 0.0, 1.0, 0.0], [-1.0, 1.0, 0.0] ];
basis2s = [ [-0.5, 1.0, 0.0], [ 1.0,-0.5, 0.0], [ 0.5, 0.5, 0.0] ];

sqw_3Q_0 = zeros(Float64,3,N,N,length(Elist));  sqw_1Q_0 = zeros(Float64,3,N,N,length(Elist));

ratios = 0.50:0.05:0.90;
sqw_3Q_1 = zeros(Float64,3,N,N,length(ratios));  sqw_1Q_1 = zeros(Float64,3,N,N,length(ratios));
sqw_3Q_2 = zeros(Float64,3,N,N,length(ratios));  sqw_1Q_2 = zeros(Float64,3,N,N,length(ratios));

Threads.@threads for (q,(basis1, basis2)) in enumerate(zip(basis1s, basis2s))
  qptsL001 = [ shift + x * basis1 + y * basis2  for x in xGrd, y in yGrd ];
  global sqw_3Q_0[q,:,:,:] = intensities_interpolated(sc3Q, qptsL001, formula; interpolation = :linear);
  global sqw_1Q_0[q,:,:,:] = intensities_interpolated(sc1Q, qptsL001, formula; interpolation = :linear);
  for (i,ratio) in enumerate(ratios)
    Ectr1 = ratio * J1 * (1 + j2) - 0.3;  _, idx1 = findmin(abs.(Elist .- Ectr1));
    Ectr2 = ratio * J1 * (1 + j2) + 0.3;  _, idx2 = findmin(abs.(Elist .- Ectr2));
    global sqw_3Q_1[q,:,:,i] = mean(sqw_3Q_0[q,:,:,idx1:idx2], dims=3);
    global sqw_1Q_1[q,:,:,i] = mean(sqw_1Q_0[q,:,:,idx1:idx2], dims=3);
    global sqw_3Q_2[q,:,:,i] = broaden2D(sqw_3Q_1[q,:,:,i], xGrd, yGrd, 0.020, 0.023);
    global sqw_1Q_2[q,:,:,i] = broaden2D(sqw_1Q_1[q,:,:,i], xGrd, yGrd, 0.020, 0.023);
  end
end

# =======================================
#? Check the data with the plot
# =======================================

data2Dcut3Qs = zeros(Float64, length(ratios), N, N);  
data2Dcut1Qs = zeros(Float64, length(ratios), N, N);

for (i,ratio) in enumerate(ratios)
  Ectr1 = ratio * J1 * (1 + j2) - 0.3;  _, idx1 = findmin(abs.(Elist .- Ectr1));
  Ectr2 = ratio * J1 * (1 + j2) + 0.3;  _, idx2 = findmin(abs.(Elist .- Ectr2));
  
  tmp1 = sqw_3Q_2[:,:,:,i];
  tmp2 = reverse(tmp1, dims = 2);  tmp2 = reverse(tmp2, dims = 3);
  tmp3 = mean(tmp1 + tmp2, dims = 1);  tmp3 = tmp3[1,:,:,1];
  
  data2Dcut3Q = tmp3;
  data2Dcut3Q = 0.5 * (data2Dcut3Q + reverse(data2Dcut3Q, dims = 1));
  data2Dcut3Q = 0.5 * (data2Dcut3Q + reverse(data2Dcut3Q, dims = 2));
  data2Dcut3Qs[i,:,:] = data2Dcut3Q;
  
  tmp1 = sqw_1Q_2[:,:,:,i];
  tmp2 = reverse(tmp1, dims = 2);  tmp2 = reverse(tmp2, dims = 3);
  tmp3 = mean(tmp1 + tmp2, dims = 1);  tmp3 = tmp3[1,:,:,1];
  
  data2Dcut1Q = tmp3;
  data2Dcut1Q = 0.5 * (data2Dcut1Q + reverse(data2Dcut1Q, dims = 1));
  data2Dcut1Q = 0.5 * (data2Dcut1Q + reverse(data2Dcut1Q, dims = 2));
  data2Dcut1Qs[i,:,:] = data2Dcut1Q;
  
  figure = CairoMakie.Figure();  ax = CairoMakie.Axis(figure[1,1]; aspect = 2/√3);
  heatmap!(ax, xGrd, yGrd, data2Dcut3Q; colormap = :viridis);
  xlims!(ax,-0.5,0.5);  ylims!(ax,-0.5,0.5);
  filename = @sprintf("data2Dcut3Q_%.3f.png",ratio);  save(filename, figure);
  
  figure = CairoMakie.Figure();  ax = CairoMakie.Axis(figure[1,1]; aspect = 2/√3);
  heatmap!(ax, xGrd, yGrd, data2Dcut1Q; colormap = :viridis);
  xlims!(ax,-0.5,0.5);  ylims!(ax,-0.5,0.5);
  filename = @sprintf("data2Dcut1Q_%.3f.png",ratio);  save(filename, figure);
end

filename = @sprintf("constEcut_%.2f_%.2f.h5",j2,jc2);
h5open(filename, "w") do file
  write(file,"J1", J1);
  write(file,"j2", j2);  write(file,"j3", j3);
  write(file,"jc1", jc1);  write(file,"jc2", jc2);
  write(file, "xGrd", xGrd);  write(file, "yGrd", yGrd);
  write(file, "data2Dcut3Qs", data2Dcut3Qs);  write(file, "data2Dcut1Qs", data2Dcut1Qs);
end;