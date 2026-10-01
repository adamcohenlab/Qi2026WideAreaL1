% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Author: Yitong Qi.
% Origin: the author's MATLAB voltage-imaging library (imaging/voltage).
% --------------------------------------------------------------------------
function [ftprnt,spkTrace] = cellSpkFtprntInit2(mov,initFtprnt,mcTrace,nPkLim)

    trace = zscore(SeeResiduals_vec(mov',initFtprnt,0,1)',[],2);
    traceFilt = trace - movmedian(trace,15,2);
    traceFilt = fft_clean(traceFilt')';
%     trace(1,:) = SeeResiduals_vec(trace(1,:),trace(2:end,:));
    [traceFilt,w] = SeeResiduals_vec(traceFilt(1,:),traceFilt(2:end,:),0);
    trace = trace(1,:) - w*trace(2:end,:);
%     trace = SeeResiduals_vec(trace(1,:),trace(2:end,:));
%     traceOut = trace;
%     trace = mov' * initFtprnt;
%     trace = trace(1,:);
    
%     traceFilt = trace_motion_correction(traceFilt,mcTrace,'aggressive');

    [pk,locs] = findpeaks(traceFilt,'sortstr','descend');
%     [~,nLocs] = findpeaks(-traceFilt,'sortstr','descend');
%     thres = -traceFilt(nLocs(25));
    thres = -prctile(traceFilt,.1);
%     locsOut = locs(1:nPkLim);
%% filter out bad peaks 
%     nWinPk = 5;
%     locs = locs(pk > thres);
%     pkWav = traceFilt(locs + [-nWinPk:nWinPk]');
%     [u,s,v] = svds(double(pkWav));
%     distWavFt = mahal(v(:,[2 3]),v(:,[2 3]));
% 
%     locs = locs(distWavFt < median(distWavFt)*2);
% 
%     spkTrace = zeros(size(trace));
%     spkTrace(locs + [-nWinPk:nWinPk]') = trace(locs + [-nWinPk:nWinPk]');
%%
    spkTrace = zeros(size(trace));
    spkTrace(traceFilt>thres) = trace(traceFilt>thres);
    [~,ftprnt] = SeeResiduals_vec(mov,spkTrace);
    ftprnt = max(ftprnt,0);

    spkTrace = spkTrace * squeeze(mean(ftprnt(:,2),[1 2]));
    ftprnt = ftprnt(:,2)/mean(ftprnt(:,2),[1 2]);
    
end