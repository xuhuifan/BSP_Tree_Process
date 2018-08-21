function truncated_particle = truncate_current_particle(currentParticle, flag_ii)
% truncate the currentParticle to the stage_ii of flag_ii

truncated_particle = currentParticle;
truncated_particle.kdtree((2*flag_ii):end, :) = 0;
truncated_particle.kdtree(truncated_particle.kdtree(:, 2)>(2*flag_ii-1), 2:3) = 0;

truncated_particle.nodeNum = 2*flag_ii-1;


truncated_particle.points((2*flag_ii+3):end, :) = 0;

truncated_particle.pIndex((2*flag_ii):end) = [];

truncated_particle.pointNum = 2*flag_ii+2;

truncated_particle.ll_ratio_seq((flag_ii+1):end) = [];

truncated_particle.remainingBudget((flag_ii+1):end) = [];

truncated_particle.z_label((flag_ii+1):end) = [];







