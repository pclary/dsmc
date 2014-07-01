function [findtimes, targets, phi2, mur, c, xt1, xt2] = runsearch(K, v1, v2, h, x1, x2, mu, nsamplepts, ...
    xt1i, xt2i, umax, algorithm, ntargets, findradius, findtimeconst, maxtime, ...
    lambda, axmain, gifname, axaux, axmu, axcov, aborthandle, stopallfound, ...
    pausebutton, au, tu, spherical)
%RUNSEARCH Runs a search using the specified algorithm
%   Returns the time at which each target was found and a vector of the
%   number of remaining targets at each detection time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/14/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim = [min(x1(:)), max(x1(:))];
ylim = [min(x2(:)), max(x2(:))];

% Precompute Fourier basis functions
fks = cell(K + 1);
hks = zeros(K + 1);
for K1 = 0:K
    for K2 = 0:K
        k1 = K1*pi/(xlim(2)-xlim(1));
        k2 = K2*pi/(ylim(2)-ylim(1));
        hks(K1+1, K2+1) = sqrt(trapz(x2(:,1), trapz(x1(1,:), cos(k1*(x1-xlim(1))).^2.*cos(k2*(x2-ylim(1))).^2, 2), 1));
        fks{K1+1, K2+1} = 1./hks(K1+1, K2+1).*cos(k1*(x1-xlim(1))).*cos(k2*(x2-ylim(1)));
    end
end

% Sample distribution to get particles and targets
[mu1, mu2] = sampledist(mu, x1(1, :), x2(:, 1), nsamplepts);
[mu1tar, mu2tar] = sampledist(mu, x1(1, :), x2(:, 1), ntargets);
foundtargets = [];
findtimes = 0;
findchance = exp(-h/findtimeconst);
phi2 = [];

% Set initial positions
xt1 = xt1i;
xt2 = xt2i;
if strcmp(algorithm, 'Lawnmower')
    xt1c = xt1i;
    xt2c = xt2i;
end

fcount = 0;
if h > 1
    frate = 1/30;
elseif h > 1/30
    frate = h;
elseif h > 1/60
    frate = h/2;
else
    frate = 1/30;
end
stepnum = size(xt1, 1);
t = 0;
co = get(gca,'ColorOrder');

% Used for abort button signalling with the GUI
if nargin >= 23 && ~isempty(aborthandle)
    abort = @() getappdata(aborthandle, 'Abort');
else
    abort = @() 0;
end

% Used for pause button signalling with the GUI
if nargin >= 25 && ~isempty(pausebutton)
    paused = @() get(pausebutton, 'Value');
else
    paused = @() 0;
end

if nargin < 24 || isempty(stopallfound)
    stopallfound = 0;
end

% Step through simulation until a stop condition is met
while t < maxtime && (~stopallfound || numel(foundtargets) < ntargets) && ~abort()
    
    % Step particles, targets, and trajectories forward in time according 
    % to the velocity field
    [mu1, mu2] = rk4step(t, h, mu1, mu2, v1, v2);
    [mu1tarn, mu2tarn] = rk4step(t, h, mu1tar, mu2tar, v1, v2);
    tdisps = sqrt((mu1tarn-mu1tar).^2 + (mu2tarn-mu2tar).^2);
    mu1tar = mu1tarn + tu*rand(size(mu1tar)).*tdisps;
    mu2tar = mu2tarn + tu*rand(size(mu1tar)).*tdisps;
    [xt1, xt2] = rk4step(t, h, xt1, xt2, v1, v2);
    
    % Use search algorithm to determine new agent positions
    if strcmp(algorithm, 'DSMC')
        xt10 = xt1(1:stepnum, :);
        xt20 = xt2(1:stepnum, :);
        muks = logmuks(x1, x2, mu1, mu2, lambda, K, fks, hks);
        [xt1n, xt2n] = dsmcstep(xt10, xt20, muks, h, umax, hks, K, xlim, ylim, au, spherical);
%         for i = 1:numel(xt1i)
%             [xt1n(i), xt2n(i)] = dsmcstep(xt10(:, i), xt20(:, i), muks, h, umax, hks, K, L1, L2, au, tu);
%         end
    elseif strcmp(algorithm, 'Lawnmower')
        xt10 = xt1c(1:stepnum, :);
        xt20 = xt2c(1:stepnum, :);
        [xt1n, xt2n] = lawnmowerstep(xt10, xt20, h, umax, xlim, ylim, 10, spherical);
        xt1c(stepnum+1, :) = xt1n;
        xt2c(stepnum+1, :) = xt2n;
        if ((nargin >= 20 && ~isempty(axaux)) || (nargin >= 21&& ~isempty(axmu))) && ...
                mod(stepnum, 10) == 1
            muks = logmuks(x1, x2, mu1, mu2, lambda, K, fks, hks);
        end
    elseif strcmp(algorithm, 'Random Walk')
        xt10 = xt1(1:stepnum, :);
        xt20 = xt2(1:stepnum, :);
        [xt1n, xt2n] = randomwalkstep(xt10, xt20, h, umax, xlim, ylim, spherical);
        if ((nargin >= 20 && ~isempty(axaux)) || (nargin >= 21&& ~isempty(axmu))) && ...
                mod(stepnum, 10) == 1
            muks = logmuks(x1, x2, mu1, mu2, lambda, K, fks, hks);
        end
    end
    xt1(stepnum+1, :) = xt1n;
    xt2(stepnum+1, :) = xt2n;
    stepnum = stepnum + 1;
    
    t = t + h;
    
    % Keep track of discovered targets
    notfound = ~ismember(1:ntargets, foundtargets);
    dists = min(sqrt(bsxfun(@minus, xt1(stepnum, :), mu1tar).^2 + ...
        bsxfun(@minus, xt2(stepnum, :), mu2tar).^2), [], 2);
    found = find(dists < findradius & notfound' & rand(ntargets, 1) > findchance);
    foundtargets = [foundtargets; found];
    findtimes = [findtimes; t*ones(numel(found), 1)];
    
    % Plot particles and trajectories
    if nargin >= 18 && ~isempty(axmain)
        plot(axmain, mu1, mu2, '.b', 'MarkerSize', 2);
        set(axmain, 'ColorOrder', co(2:end, :), 'NextPlot', 'replacechildren');
        hold(axmain, 'on');
        plot(axmain, mu1tar, mu2tar, 'co', 'MarkerSize', 5, 'LineWidth', 2);
        plot(axmain, mu1tar(foundtargets), mu2tar(foundtargets), 'mo', ...
            'MarkerSize', 5, 'LineWidth', 2);
        
        plot(axmain, [xt1(1, :); xt1(1:stepnum, :)], [xt2(1, :); xt2(1:stepnum, :)]);
        plot(axmain, [xt1(stepnum, :); xt1(stepnum, :)], ...
            [xt2(stepnum, :); xt2(stepnum, :)], '.', 'MarkerSize', 15);
        
        hold(axmain, 'off');
        axis(axmain, 'equal');
        axis(axmain, [xlim, ylim]);
        title(axmain, sprintf('t = %.2f s', t));
    end
    
    % Save trajectory animation
    if nargin >= 19 && ~isempty(gifname) && numel(gifname) > 0 && t >= fcount*frate
        frame=getframe(axmain);
        im=frame2im(frame);
        [imind,map]=rgb2ind(im,256);
        if fcount == 0
            imwrite(imind,map,gifname,'DelayTime',frate,'LoopCount',inf);
        else
            imwrite(imind,map,gifname, 'DelayTime',frate, 'WriteMode', 'append');
        end
        fcount = fcount+1;
    end
    
%     % Plot DSMC surface
%     if strcmp(algorithm, 'DSMC') && nargin >= 20 && ~isempty(axaux)
%         sk = @(K1, K2) ck(K1, K2, xt1(1:stepnum, :), xt2(1:stepnum, :), ...
%             hks(K1+1, K2+1), L1, L2) - muks(K1+1, K2+1);
%         s = zeros(size(mu));
%         for K1 = 0:K
%             for K2 = 0:K
%                 k1 = K1*pi/L1;
%                 k2 = K2*pi/L2;
%                 La = 1/(1 + K1^2 + K2^2)^(3/2);
%                 s = s + fks{K1+1, K2+1}*sk(K1, K2)*La;
%             end
%         end
%         surf(axaux, x1, x2, s);
%     end

    % Plot convergence metric
    if nargin >= 20 && ~isempty(axaux)
        sk = @(K1, K2) ck(K1, K2, xt1(1:stepnum, :), xt2(1:stepnum, :), ...
            hks(K1+1, K2+1), xlim, ylim) - muks(K1+1, K2+1);
        phi = 0;
        for K1 = 0:K
            for K2 = 0:K
                La = 1/(1 + K1^2 + K2^2)^(3/2);
                phi = phi + La*sk(K1, K2)^2;
            end
        end
        phi2 = [phi2; phi];
        loglog(axaux, phi2);
        title(axaux, 'Convergence');
        xlabel(axaux, 'Steps');
        ylabel(axaux, '\phi^2');
    end

    
    % Plot log density of particles
    if nargin >= 21 && ~isempty(axmu)
        mur = zeros(size(mu));
        for K1 = 0:K
            for K2 = 0:K
                mur = mur + fks{K1+1, K2+1}*muks(K1+1, K2+1);
            end
        end
        mur = mur / trapz(x2(:, 1), trapz(x1(1, :), mur, 2));
        surf(axmu, x1, x2, mur, 'EdgeColor', 'none');
        axis(axmu, 'equal');
        axis(axmu, [xlim, ylim]);
        title(axmu, 'Log(Particle distribution)');
        caxis(axmu, sort([0, max(max(mur))]));
    end
    
    % Plot trajectory density
    if nargin >= 22 && ~isempty(axcov)
        c = zeros(size(mu));
        for K1 = 0:K
            for K2 = 0:K
                c = c + fks{K1+1, K2+1}*ck(K1, K2, xt1(1:stepnum, :), ...
                    xt2(1:stepnum, :), hks(K1+1, K2+1), xlim, ylim);
            end
        end
        c = c / trapz(x2(:, 1), trapz(x1(1, :), c, 2));
        surf(axcov, x1, x2, c, 'EdgeColor', 'none');
        axis(axcov, 'equal');
        axis(axcov, [xlim, ylim]);
        title(axcov, 'Coverage density');
        caxis(axcov, [0, max(max(c))]);
    end
    
    drawnow;
    
    % Pause
    while paused()
        pause(0.05);
    end
end

set(axmain, 'ColorOrder', co);

targets = ntargets:-1:ntargets-numel(findtimes)+1;

end

