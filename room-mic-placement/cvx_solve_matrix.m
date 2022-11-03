function [G] = cvx_solve_matrix(G_prev, G_next, G_cur)
% this is a semidefinite program because of LMI (linear matrix inequality)

[Nmic, Nsrc] = size(G_cur);
LL = max(G_prev,G_next);
cvx_begin sdp
    variable G(Nmic,Nsrc)
    minimize( norm( G - G_cur, 'fro' ) )
    vec(G) >= vec(LL);
cvx_end
end

