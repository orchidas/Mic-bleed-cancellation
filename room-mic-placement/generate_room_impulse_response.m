function [h_ir,H_ir] = generate_room_impulse_response(fs, dims, src_pos, rec_pos, rt60, nsamp,varargin)
%%
% Generate RIR matrix from Nsrc to Nmic with IM 
% Returns
% h_ir - RIR matrix in time-domain
% H_ir - RIR matrix in frequency domain
%%

addpath(genpath('../room_acoustics_simulation/RIR-Generator/.'));   %rir toolbox

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'orientation',[0,0]);
addParameter(p,'mtype','omnidirectional');
addParameter(p,'ref_order',3);
addParameter(p, 'hp_filter', 1);

parse(p,varargin{:});
orientation = p.Results.orientation;
mtype = p.Results.mtype;
order = p.Results.ref_order;
hp_filter = p.Results.hp_filter;

Nsrc = size(src_pos,1);
Nmic = size(rec_pos,1);
dim = size(src_pos,2);
c = 343;
h_ir = zeros(nsamp, Nmic, Nsrc);




%% RIR generator

for k = 1:Nsrc
    h = rir_generator(c, fs, rec_pos, src_pos(k,:), dims, rt60, nsamp, mtype, order, dim, orientation, hp_filter);
    h_ir(:,:,k) = h.';

end

H_ir = fft(h_ir, 2^ceil(log2(nsamp)), 1);

