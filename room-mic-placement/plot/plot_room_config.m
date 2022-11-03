function [ax] = plot_room_config(room_dims, src_pos, rec_pos, varargin)


p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'marker_color','k');
addParameter(p, 'title' ,[]);
parse(p,varargin{:});
col = p.Results.marker_color;
tit = p.Results.title;

%% plot mic and source configuration
    figure;
    Nsrc = size(src_pos,1);
    sp = scatter3(src_pos(:,1), src_pos(:,2), src_pos(:,3));hold on;grid on;
    sp.MarkerEdgeColor = 'b';
    scatter3(rec_pos(1:Nsrc,1), rec_pos(1:Nsrc,2), rec_pos(1:Nsrc,3), 'rx');hold on;
    rp = scatter3(rec_pos(Nsrc+1:end,1), rec_pos(Nsrc+1:end,2), rec_pos(Nsrc+1:end,3), 'x');hold on;
    rp.MarkerFaceColor = col;
    xlim([0 room_dims(1)]);ylim([0 room_dims(2)]);zlim([0 room_dims(3)]);
    hold off;
    drawnow;
    ax = gca;
    title(tit);
    legend('Sources','Mics');

end

