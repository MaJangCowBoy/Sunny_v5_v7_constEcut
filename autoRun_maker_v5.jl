using Printf

j2_sweep  = 0.4:0.1:0.8;  jc2_sweep =-0.1:0.1:0.7;

for (i,j2) in enumerate(j2_sweep), (j,jc2) in enumerate(jc2_sweep)
  id = (i-1) * length(jc2_sweep) + j;
  fid = open("autoBatch_$(id).sh","w");
  @printf(fid,"#!/bin/bash\n");
  @printf(fid,"#SBATCH --cpus-per-task=1\n");
  @printf(fid,"#SBATCH --job-name=woonghee_CTS\n");
  @printf(fid,"source ~/.bashrc\n");
  @printf(fid,"julia_wh master_script_v5.jl %f %f 1> so_%d.txt 2> se_%d.txt\n",j2,jc2,id,id);
  close(fid);
end

fid = open("autoRun_v5.sh","w");
@printf(fid,"#!/bin/bash\n")
for (i,j2) in enumerate(j2_sweep), (j,jc2) in enumerate(jc2_sweep)
  id = (i-1) * length(jc2_sweep) + j;
  @printf(fid,"sbatch autoRun_%d.sh\n",id);
end
close(fid)