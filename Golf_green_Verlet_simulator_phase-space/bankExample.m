function bankExample(angleDeg, speed, numBalls, mode, startLoc, isReplay, showTrail)
    if nargin < 6
        isReplay = false;
    end
    if nargin < 7
        showTrail = false;
    end

%% Simulation setup
if ~isReplay
    close all;  % Only close figures for the main simulation, not replays
end
clc;

dt = 0.01;
Cg = 50; Cd = 1; minSpeed = 10;
holeRadius = 6;
name = 'greenbankexample.png.tiff';
greenheight = fliplr(rot90(rot90(double(imread(name)) / 256 * 0.3)));
[gradX, gradY] = gradient(greenheight(:,:,1));
forceFun = @(R, RL) (RL - R)/dt * Cd + ...
    [interp2(-gradX, R(1,:), R(2,:)); interp2(-gradY, R(1,:), R(2,:))] * Cg;
holeLoc = size(greenheight)' / 2;
l = size(greenheight, 1);

%% Angle and speed configuration
if strcmp(mode, "Manual")
    angleRad = angleDeg * pi / 180;
    if numBalls == 1
        angles = angleRad;
        speeds = speed;
    else
        spread = 10 * pi / 180;
        angles = linspace(angleRad + spread/2, angleRad - spread/2, numBalls)';
        speeds = repmat(speed, numBalls, 1);  % Column vector
    end
else
    angleSpread = 15;
    speedSpread = 200;
    angleRange = linspace(angleDeg + angleSpread/2, angleDeg - angleSpread/2, numBalls);
    speedRange = linspace(speed + speedSpread/2, speed - speedSpread/2, numBalls);
    [A, S] = meshgrid(angleRange, speedRange);
    angles = A(:) * pi / 180;
    speeds = S(:);
end

%% Launch setup
if strcmp(mode, 'Ideal Finder')
    numParticles = numBalls * numBalls;
else
    numParticles = length(angles);
end

angleVectors = [cos(angles)'; sin(angles)'];
startConditions = angleVectors .* repmat(speeds', 2, 1);
r = repmat(startLoc, 1, numParticles);
rl = r - startConditions * dt;
rstart = r;
rlstart = rl;

if showTrail
    trailPoints = [];
end

%% Simulation loop
i = 0; launchTicks = 80;
while sum(sum(r ~= rl))
    i = i + 1;
    rn = 2 * r - rl + forceFun(r, rl) * dt^2;
    rl = r; r = rn;

    if i < 100
        r = rstart; rl = rlstart;
    end

    dr = rl - r;
    s = sqrt(dr(1,:).^2 + dr(2,:).^2) / dt;
    haltBallsEdge = (r > l) | (r < 1);
    distToHole = sqrt(sum((rl - holeLoc).^2, 1));
    holeBalls = (distToHole < holeRadius);
    rl(:, holeBalls) = repmat(holeLoc, 1, sum(holeBalls));
    haltBalls = (haltBallsEdge(1,:) | haltBallsEdge(2,:)) | ...
                 (s < minSpeed) | (distToHole < holeRadius);
    r(:, haltBalls) = rl(:, haltBalls);

    if showTrail && numParticles == 1
        trailPoints(:, end+1) = r(:,1);
    end

    if mod(i,10) == 0
        figure(1); clf
        surf(greenheight, 'edgecolor', 'none'); hold on;
        colormap(flipud([181,228,138;160,220,104;124,216,87;81,212,76; ...
                         52,188,67;37,167,60;32,154,61;1,114,56;0,86,19]/255));
        set(gca, 'DataAspectRatio', [1 1 1])

        inHole = distToHole < holeRadius;
        visible = ~inHole;

        if any(visible)
            scatter3(r(1,visible), r(2,visible), ...
                3 + greenheight(sub2ind(size(greenheight), ...
                floor(r(2,visible)), floor(r(1,visible)))), 'wo', 'filled');
        end

        if isReplay && any(inHole)
            sinkZ = -5 * ones(1, sum(inHole));
            scatter3(r(1,inHole), r(2,inHole), sinkZ, 'ko', 'filled');
        end

        if showTrail && exist('trailPoints', 'var') && size(trailPoints,2) > 1
            trailZ = 3 + greenheight(sub2ind(size(greenheight), ...
                floor(trailPoints(2,:)), floor(trailPoints(1,:))));
            plot3(trailPoints(1,:), trailPoints(2,:), trailZ, 'r-', 'LineWidth', 2);
        end

        holeX = [holeLoc(1), holeLoc(1)];
        holeY = [holeLoc(2), holeLoc(2)];
        holeZ = [greenheight(floor(holeLoc(1)),floor(holeLoc(2))), ...
                 greenheight(floor(holeLoc(1)),floor(holeLoc(2))) + 100];
        plot3(holeX, holeY, holeZ, 'k-');
        view([-39 30.85]); xlim([0 512]); ylim([0 512]);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    end
end

%% Output result
finalDistance = distToHole(1);
scaleFeetPerUnit = (4.25 / 12) / holeRadius;
realDistance = finalDistance * scaleFeetPerUnit;

if strcmp(mode, 'Manual') && numBalls == 1
    if isReplay
        fprintf('\nIdeal Putt Found:\n');
        fprintf('Angle: %.4f degrees\n', angles(1) * 180/pi);
        fprintf('Speed: %.4f units\n', speeds(1));
    else
        if finalDistance < holeRadius
            fprintf('The ball went in the hole!\n');
        else
            fprintf('The ball was %.2f feet from the hole.\n', realDistance);
        end
        fprintf('Angle: %.4f degrees\n', angles(1) * 180/pi);
        fprintf('Speed: %.4f units\n', speeds(1));
    end
elseif strcmp(mode, 'Manual')
    centerIdx = ceil(numBalls / 2);
    if finalDistance < holeRadius
        fprintf('The ball went in the hole!\n');
    else
        fprintf('The center ball was %.2f feet from the hole.\n', realDistance);
    end
    fprintf('Angle: %.4f degrees\n', angles(centerIdx) * 180/pi);
    fprintf('Speed: %.4f units\n', speeds(centerIdx));
end

%% Phase space visualization (Ideal Finder)
if strcmp(mode, 'Ideal Finder')
    dists = reshape(distToHole, numBalls, numBalls);
    figure(2)
    imshow(dists, [holeRadius, 256], 'colormap', colormap('jet'))
    hold on; h = imshow(repmat(ones(size(dists)), 1, 1, 3)); hold off
    set(h, 'AlphaData', dists < holeRadius);
    ylabel('Speed index'); xlabel('Angle index');

    [~, bestIdx] = min(distToHole);
    bestAngle = angles(bestIdx);
    bestSpeed = speeds(bestIdx);
    fprintf('\nIdeal Putt Found:\n');
    fprintf('Angle: %.4f degrees\n', bestAngle * 180/pi);
    fprintf('Speed: %.4f units\n', bestSpeed);

    pause(1);
    fprintf('Replaying best shot...\n');
    bankExample(bestAngle * 180/pi, bestSpeed, 1, 'Manual', startLoc, true, true);
end
end
