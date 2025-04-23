function interactiveLauncher()
% interactiveLauncher - GUI for selecting putt parameters and launching simulation

f = uifigure('Name', 'Golf Putt Simulator', 'Position', [100 100 600 650]);
f.UserData.startLoc = [64; 64];  % Default start location stored in UserData

%% Green Preview for Clickable Start Location
ax = uiaxes(f, 'Position', [50 340 512 220]);
greenImage = fliplr(rot90(rot90(double(imread('greenbankexample.png.tiff')) / 256 * 0.3)));
imagesc(ax, greenImage); axis(ax, 'image'); title(ax, 'Click to select start location');

%% Launch Angle Slider
uilabel(f, 'Text', 'Launch Angle (°)', 'Position', [30 300 150 22]);
angleSlider = uislider(f, 'Limits', [0 90], 'Value', 45, 'Position', [30 285 540 3]);

%% Launch Speed Slider
uilabel(f, 'Text', 'Launch Speed', 'Position', [30 240 150 22]);
speedSlider = uislider(f, 'Limits', [100 600], 'Value', 400, 'Position', [30 225 540 3]);

%% Number of Balls Slider
uilabel(f, 'Text', 'Number of Balls', 'Position', [30 180 150 22]);
ballSlider = uislider(f, 'Limits', [1 20], 'MajorTicks', 1:1:20, 'Value', 1, 'Position', [30 165 540 3]);

%% Mode Selector
uilabel(f, 'Text', 'Mode', 'Position', [30 100 150 22]);
modeDropdown = uidropdown(f, 'Items', {'Manual', 'Ideal Finder'}, 'Position', [100 80 120 22]);

%% Run Simulation Button
uibutton(f, 'Text', 'Run Simulation', 'Position', [250 30 140 30], ...
    'ButtonPushedFcn', @(btn, event) runSimCallback());

%% Stimp Test Button
uibutton(f, 'Text', 'Run Stimp Test', 'Position', [420 30 140 30], ...
    'ButtonPushedFcn', @(btn, event) runStimpTestFromGUI(f.UserData.startLoc));

%% Callback Function
function runSimCallback()
    angle = angleSlider.Value;
    speed = speedSlider.Value;
    numBalls = round(ballSlider.Value);
    mode = modeDropdown.Value;
    startLoc = f.UserData.startLoc;

    % If Ideal Finder, inject extra argument to show trail
    if strcmp(mode, 'Ideal Finder')
        bankExample(angle, speed, numBalls, mode, startLoc, false, true);
    else
        bankExample(angle, speed, numBalls, mode, startLoc);
    end
end

%% Enable click-to-set-start-loc
f.WindowButtonDownFcn = @figureClickCallback;

function figureClickCallback(~, ~)
    clickedObj = gco;
    if isa(clickedObj, 'matlab.graphics.axis.Axes') && clickedObj == ax
        cp = ax.CurrentPoint;
        x = round(cp(1,1));
        y = round(cp(1,2));

        % Clip to bounds of image
        x = max(1, min(x, size(greenImage, 2)));
        y = max(1, min(y, size(greenImage, 1)));

        f.UserData.startLoc = [x; y];
        title(ax, sprintf('Start location: [%d, %d]', x, y));

        % Optional marker
        hold(ax, 'on');
        plot(ax, x, y, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    end
end
end

function runStimpTestFromGUI(startLoc)
    % Parameters
    dt = 0.01;
    Cg = 50; Cd = 1; minSpeed = 10;
    holeRadius = 6;
    scaleFeetPerUnit = (4.25 / 12) / holeRadius;

    % Load real green (same as bankExample)
    name = 'greenbankexample.png.tiff';
    greenheight = fliplr(rot90(rot90(double(imread(name)) / 256 * 0.3)));
    [gradX, gradY] = gradient(greenheight);

    forceFun = @(R, RL) (RL - R)/dt * Cd + ...
        [interp2(-gradX, R(1,:), R(2,:)); interp2(-gradY, R(1,:), R(2,:))] * Cg;

    % Launch at 6 ft/s in sim units
    angleDeg = 0;
    angleRad = angleDeg * pi / 180;
    simSpeed = 6.0 / scaleFeetPerUnit;

    r = startLoc;
    rl = r - [cos(angleRad); sin(angleRad)] * simSpeed * dt;

    i = 0;
    while any(r ~= rl)
        i = i + 1;
        rn = 2*r - rl + forceFun(r, rl)*dt^2;
        rl = r; r = rn;

        dr = rl - r;
        s = norm(dr) / dt;

        if s < minSpeed
            r = rl;
            break;
        end
    end

    rollDistUnits = norm(r - startLoc);
    rollDistFeet = rollDistUnits * scaleFeetPerUnit;

    fprintf('\nStimp Test at Selected Location:\n');
    fprintf('Ball rolled: %.2f simulation units\n', rollDistUnits);
    fprintf('Stimp Rating: %.2f feet\n', rollDistFeet);
end
