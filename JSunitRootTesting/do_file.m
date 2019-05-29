clear; close all;

% Parameterisations
fid = fopen('ParameterValues.mod','wt');
fprintf(fid, 'thetaN = 10;\n');
fclose(fid);

% For the first call
fid = fopen('LoadOptions.mod','wt');
fprintf(fid, '@#define LoadSteadyState = 0\n');
fprintf(fid, '@#define Simulate = 0\n');
fclose(fid);

dynare DynamicSpatialModel

% Second call and simulation
fid = fopen('LoadSteadyStateOption.mod','wt');
fprintf(fid, '@#define LoadSteadyState = 1\n');
fprintf(fid, '@#define Simulate = 1\n');
fclose(fid);

dynare DynamicSpatialModel mingw
