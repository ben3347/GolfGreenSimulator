function interactiveLauncher()
% interactiveLauncher - GUI for selecting putt parameters and launching simulation

f = uifigure('Name', 'Golf Putt Simulator', 'Position', [100 100 600 650]);

%% Green Preview for Clickable Start Location
ax = uiaxes(f, 'Position', [50 340 512 220]);
greenImage = fliplr(rot90(rot90(double(imread('greenbankexample.png.tiff')) / 256 * 0.3)));
imagesc(ax, greenImage); axis(ax, 'image'); title(ax, 'Click to select start location');

% Default start location
startLoc = [64; 64];

% Callback to capture click
ax.ButtonDownFcn = @(src, event) setStart(event);

function setStart(event)
    coords = round(event.IntersectionPoint(1:2));
    startLoc = [coords(1); coords(2)];
    title(ax, sprintf('Start location: [%d, %d]', coords(1), coords(2)));
end

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

%% Callback Function
function runSimCallback()
    angle = angleSlider.Value;
    speed = speedSlider.Value;
    numBalls = round(ballSlider.Value);
    mode = modeDropdown.Value;

    bankExample(angle, speed, numBalls, mode, startLoc);
end
end
