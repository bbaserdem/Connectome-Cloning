function [ conv, dens ] = conclone_wrapper( N, f )
%CONCLONE_WRAPPER Wrapper that takes in N and average f values for output

[conn, dens] = init_conn( N, 'Erdos', f );
conv = simulation_flicker( conn );

end

