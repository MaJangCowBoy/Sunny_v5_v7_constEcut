############################################################################

J1 = 1.0;  Jc1 = 1.3610;
filename = @sprintf("LT_minimize_%.4f",Jc1);
filename = replace(filename,"." => "p");  filename = filename*".dat";

J2stt  = 0.25;  J2end  = 0.75;  J2stp  = ( J2end- J2stt)/100;  J2arr  = J2stt:J2stp:J2end;
Jc2stt =-0.25;  Jc2end =+0.25;  Jc2stp = (Jc2end-Jc2stt)/100;  Jc2arr = Jc2stt:Jc2stp:Jc2end;
# J3stt =-0.10;  J3stp = 0.01;  J3end = 1.00;  J3arr = J3stt:J3stp:J3end;


############################################################################

A_1NN_1 = [+1.0; 0.0; 0.0];  A_1NN_2 = [-1.0; 0.0; 0.0];
A_1NN_3 = [ 0.0;+1.0; 0.0];  A_1NN_4 = [ 0.0;-1.0; 0.0];
A_1NN_5 = [+1.0;+1.0; 0.0];  A_1NN_6 = [-1.0;-1.0; 0.0];
dAA_1NN = hcat(A_1NN_1, A_1NN_2, A_1NN_3, A_1NN_4, A_1NN_5, A_1NN_6);

B_1NN_1 = [+1.0; 0.0; 0.0];  B_1NN_2 = [-1.0; 0.0; 0.0];
B_1NN_3 = [ 0.0;+1.0; 0.0];  B_1NN_4 = [ 0.0;-1.0; 0.0];
B_1NN_5 = [+1.0;+1.0; 0.0];  B_1NN_6 = [-1.0;-1.0; 0.0];
dBB_1NN = hcat(B_1NN_1, B_1NN_2, B_1NN_3, B_1NN_4, B_1NN_5, B_1NN_6);

A_2NN_1 = [+2.0;+1.0; 0.0];  A_2NN_2 = [-2.0;-1.0; 0.0];
A_2NN_3 = [+1.0;+2.0; 0.0];  A_2NN_4 = [-1.0;-2.0; 0.0];
A_2NN_5 = [-1.0;+1.0; 0.0];  A_2NN_6 = [+1.0;-1.0; 0.0];
dAA_2NN = hcat(A_2NN_1, A_2NN_2, A_2NN_3, A_2NN_4, A_2NN_5, A_2NN_6);

B_2NN_1 = [+2.0;+1.0; 0.0];  B_2NN_2 = [-2.0;-1.0; 0.0];
B_2NN_3 = [+1.0;+2.0; 0.0];  B_2NN_4 = [-1.0;-2.0; 0.0];
B_2NN_5 = [-1.0;+1.0; 0.0];  B_2NN_6 = [+1.0;-1.0; 0.0];
dBB_2NN = hcat(B_2NN_1, B_2NN_2, B_2NN_3, B_2NN_4, B_2NN_5, B_2NN_6);

A_3NN_1 = [+2.0; 0.0; 0.0];  A_3NN_2 = [-2.0; 0.0; 0.0];
A_3NN_3 = [ 0.0;+2.0; 0.0];  A_3NN_4 = [ 0.0;-2.0; 0.0];
A_3NN_5 = [+2.0;+2.0; 0.0];  A_3NN_6 = [-2.0;-2.0; 0.0];
dAA_3NN = hcat(A_3NN_1, A_3NN_2, A_3NN_3, A_3NN_4, A_3NN_5, A_3NN_6);

B_3NN_1 = [+2.0; 0.0; 0.0];  B_3NN_2 = [-2.0; 0.0; 0.0];
B_3NN_3 = [ 0.0;+2.0; 0.0];  B_3NN_4 = [ 0.0;-2.0; 0.0];
B_3NN_5 = [+2.0;+2.0; 0.0];  B_3NN_6 = [-2.0;-2.0; 0.0];
dBB_3NN = hcat(B_3NN_1, B_3NN_2, B_3NN_3, B_3NN_4, B_3NN_5, B_3NN_6);

A_c1NN_1 = [ 0.0; 0.0; 0.0];  A_c1NN_2 = [ 0.0; 0.0;-1.0];
A_c1NN_3 = [ 1.0; 0.0; 0.0];  A_c1NN_4 = [ 1.0; 0.0;-1.0];
A_c1NN_5 = [ 0.0;-1.0; 0.0];  A_c1NN_6 = [ 0.0;-1.0;-1.0];
dAB_c1NN = hcat(A_c1NN_1, A_c1NN_2, A_c1NN_3, A_c1NN_4, A_c1NN_5, A_c1NN_6);

B_c1NN_1 = [ 0.0; 0.0; 0.0];  B_c1NN_2 = [ 0.0; 0.0;+1.0];
B_c1NN_3 = [-1.0; 0.0; 0.0];  B_c1NN_4 = [-1.0; 0.0;+1.0];
B_c1NN_5 = [ 0.0;+1.0; 0.0];  B_c1NN_6 = [ 0.0;+1.0;+1.0];
dBA_c1NN = hcat(B_c1NN_1, B_c1NN_2, B_c1NN_3, B_c1NN_4, B_c1NN_5, B_c1NN_6);

A_c2NN_1 = [ 1.0; 1.0; 0.0];  A_c2NN_2 = [ 1.0; 1.0;-1.0];
A_c2NN_3 = [ 1.0;-1.0; 0.0];  A_c2NN_4 = [ 1.0;-1.0;-1.0];
A_c2NN_5 = [-1.0;-1.0; 0.0];  A_c2NN_6 = [-1.0;-1.0;-1.0];
dAB_c2NN = hcat(A_c2NN_1, A_c2NN_2, A_c2NN_3, A_c2NN_4, A_c2NN_5, A_c2NN_6);

B_c2NN_1 = [-1.0;-1.0; 0.0];  B_c2NN_2 = [-1.0;-1.0;+1.0];
B_c2NN_3 = [-1.0; 1.0; 0.0];  B_c2NN_4 = [-1.0; 1.0;+1.0];
B_c2NN_5 = [ 1.0; 1.0; 0.0];  B_c2NN_6 = [ 1.0; 1.0;+1.0];
dBA_c2NN = hcat(B_c2NN_1, B_c2NN_2, B_c2NN_3, B_c2NN_4, B_c2NN_5, B_c2NN_6);

############################################################################

function calc_matrix(J2J3JcJd, Q1Q2Q3)
  J1 = 1.0;
  J2  = J2J3JcJd[1];  J3  = J2J3JcJd[2];  
  Jc1 = J2J3JcJd[3];  Jc2 = J2J3JcJd[4];
  Q1 = Q1Q2Q3[1];  Q2 = Q1Q2Q3[2];  Q3 = Q1Q2Q3[3];
  Q = hcat(Q1,Q2,Q3);

  Emat = zeros(ComplexF64,2,2);  
  Emat[1,1] =  J1 * sum( exp.(2π*im*Q* dAA_1NN) ) + 
               J2 * sum( exp.(2π*im*Q* dAA_2NN) ) +
               J3 * sum( exp.(2π*im*Q* dAA_3NN) ) ;
  Emat[1,2] = Jc1 * sum( exp.(2π*im*Q*dAB_c1NN) ) + 
              Jc2 * sum( exp.(2π*im*Q*dAB_c2NN) ) ;
  Emat[2,1] = Jc1 * sum( exp.(2π*im*Q*dBA_c1NN) ) +
              Jc2 * sum( exp.(2π*im*Q*dBA_c2NN) ) ;
  Emat[2,2] =  J1 * sum( exp.(2π*im*Q* dBB_1NN) ) + 
               J2 * sum( exp.(2π*im*Q* dBB_2NN) ) +
               J3 * sum( exp.(2π*im*Q* dBB_3NN) ) ;
  k = eigen(Emat);
  E1 = real(k.values[1]);  E2 = real(k.values[2]);

  return minimum([E1,E2]);
end

############################################################################

function find_minpos(arr, tolerance)
  min_val = minimum(arr);  
  min_pos = [];
  for (i, val) in enumerate(arr)
    if abs(val - min_val) <= tolerance  push!(min_pos, i)  end
  end;  min_pos = convert(Array{Int64},min_pos);

  return min_pos;
end

############################################################################
