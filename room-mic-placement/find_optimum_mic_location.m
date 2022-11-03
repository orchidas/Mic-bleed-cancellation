function [m_opt] = find_optimum_mic_location(src_loc, time_delays, N_extra_mics, Nsrc, room_dims, close_mic_dist, fs)
%%
% Solve least squares problem to find optimum mic locations
%
if size(src_loc,2) ~= 3
    src_loc = src_loc.';
end

Nsquare =(Nsrc * (Nsrc-1))/2;
tau = zeros(Nsquare, N_extra_mics);
s_plus = zeros(Nsquare,3);
s_minus = zeros(Nsquare,3);
k = 1;
c = 343;
dim_tiled = repmat(room_dims,[N_extra_mics,1]).';

% form TDOA matrix
%time delays has mic along column and source among row
% but tau has mic along row and source along column
for i = 1:Nsrc-1
    for j = i+1:Nsrc
        tau(k,:) = (time_delays(:,i).^2 - time_delays(:,j).^2).';
        s_minus(k,:) = src_loc(j,:) - src_loc(i,:);
        s_plus(k,:) = src_loc(j,:) + src_loc(i,:);
        k = k+1;
    end
end

S_pm = kron(diag(s_plus*s_minus'), ones(1, N_extra_mics));
% m_opt = s_minus\(0.5*(c/fs)*tau + S_pm);

% we need to ensure that the location is not negative and lesser than room
% dimension
cvx_begin sdp
    variable m_opt(3, N_extra_mics);
    minimize norm(s_minus*m_opt - 0.5*((c^2/fs^2)*tau + S_pm), 'fro')
    subject to:
        vec(close_mic_dist*ones(3, N_extra_mics)) <= vec(m_opt) <= vec(dim_tiled);
cvx_end 

end

