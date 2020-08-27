clear;
files = dir("tests");
scripts = regexp({files.name}, ".*\.m$", "match");
scripts = vertcat(scripts{:});
failed = 0;
errors = MException.empty();
for script = scripts'
    try
        run(script{:});
    catch ME
        errors(end+1) = ME;
        failed = failed + 1;
    end
end
% Displayed if no errors occured inside the scripts
if failed == 0
    disp("The tests were ran successfully")
else
    for ME = errors
        disp(ME.message);
        disp("In file: " + ME.stack(end-2).name);
    end
    error(failed + " tests failed");
end