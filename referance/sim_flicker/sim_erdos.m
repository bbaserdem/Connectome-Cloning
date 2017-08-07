function [ conv, f ] = sim_erdos( N, d )
%SIM_ERDOS Generate Erdos random connectome and find convergance time and
%fn4

[C, f] = init_conn(N, 'Erdos', d);
conv = sim_cc( C );

end

