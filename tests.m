clear;
files = dir("tests");
scripts = regexp({files.name}, ".*\.m$", "match");
scripts = vertcat(scripts{:});
for script = scripts'
    run(script{:});
end
% Displayed if no errors occured inside the scripts
disp("The tests were ran successfully")