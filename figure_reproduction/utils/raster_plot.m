function raster_plot(spk_t_all,dt,colr,tRange,lWidth)
% RASTER_PLOT Draw a spike raster into the current axes.
%
% YQ 2021. Documentation revised 2026-09-18 for the public repository; the
% executable code below is unchanged from the lab library copy at
% matlab\raster_plot.m.
%
% Inputs
% ------
% spk_t_all : {nRows x 1} cell
%     One cell array element per raster row. Each element is a vector of
%     spike times in *sample units* (for example 1-based frame indices);
%     they are multiplied by dt before plotting. Row or column vectors are
%     both accepted. Row k is drawn spanning y = [k-0.5, k+0.5], so row
%     centers fall on integers 1..nRows.
%
% dt : scalar
%     Seconds per sample. Pass 1 if spk_t_all already holds seconds.
%
% colr : [1 x 3] or [nRows x 3] double, optional
%     RGB raster color. A single row colors every raster row and draws them
%     all into one line object; an nRows x 3 matrix colors each row
%     separately using one line object per row. Defaults to black.
%
% tRange : [1 x 2] double, optional
%     [tMin tMax] in seconds. Spikes outside the range are dropped.
%     Defaults to [-inf inf].
%
% lWidth : scalar, optional
%     Line width. Defaults to 0.5. Only applied when colr is supplied.
%
% Notes
% -----
% Spikes are drawn as a single line object per color, with NaN separators
% between ticks, which is much faster than one plot call per spike.
%
% Because the color branch is selected with nargin rather than isempty,
% call raster_plot(spk_t,dt) for the black default; do not pass an empty
% colr to request it.

% start with an empty line
is_multi_color = false;

if ~exist('tRange','var')
    tRange = [-inf inf];
end
if ~exist('lWidth','var')
    lWidth = .5;
end

if nargin == 2
    l=line(nan,nan,'color','k');
elseif nargin > 2
    if size(colr,1)==1
        l=line(nan,nan,'color',colr,'linewidth',lWidth);
    else
        is_multi_color = true;
        for ii=1:size(colr,1)
            l(ii) = line(nan,nan,'color',colr(ii,:),'linewidth',lWidth);
        end
    end
end

% fast way to fill empty line with raster data 
for i=1:length(spk_t_all)
    time_raster = spk_t_all{i}*dt;
    time_raster(time_raster<tRange(1) | time_raster > tRange(2)) = [];
%     time_raster = time_raster - tRange(1) + dt;
    
    if size(time_raster,2)>1
        time_raster = time_raster';
    end
    
    x_append = zeros(1,length(time_raster)*3);
    y_append = zeros(1,length(time_raster)*3);
    
    x_append(1:3:end) = time_raster;
    x_append(2:3:end) = time_raster;
    x_append(3:3:end) = nan;
    
    y_append(1:3:end) = ones(size(time_raster))*(i-1)+.5;
    y_append(2:3:end) = ones(size(time_raster))*(i)+.5;
    y_append(3:3:end) = nan;
    
    if is_multi_color
        l(i).XData = [l(i).XData x_append nan];
        l(i).YData = [l(i).YData y_append nan];
    else
        l.XData = [l.XData x_append nan];
        l.YData = [l.YData y_append nan];
    end
%     plot([time_raster;time_raster],[ones(size(time_raster))*(i-1);ones(size(time_raster))*(i)]+.5,'k')
%     hold on
end