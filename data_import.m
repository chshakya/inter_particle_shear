%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function structure=data_import(file_loc,cols_to_import)
%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", cols_to_import);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Intervaldata", "Time", "Gap", "Torque", "NormalForce", "DeflectionAngle", "TestTime", "EncoderDeflectionAngle", "Velocity", "ShearRate", "IntervalTime", "RotationalSpeed"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Intervaldata", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Intervaldata", "EmptyFieldRule", "auto");

% Import the data
data_numeric = readtable(file_loc, opts);


%% Clear temporary variables
clear opts

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", cols_to_import);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Intervaldata", "Time", "Gap", "Torque", "NormalForce", "DeflectionAngle", "TestTime", "EncoderDeflectionAngle", "Velocity", "ShearRate", "IntervalTime", "RotationalSpeed"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Intervaldata", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Intervaldata", "EmptyFieldRule", "auto");

% Import the data
data_string = readtable(file_loc, opts);


%% Clear temporary variables
clear opts

%%storedata to variables
linenos_interval=find(data_numeric.Intervaldata=="Number of Intervals:"); %find data on no. of intervals
interval_number=table2array(data_numeric(linenos_interval,2)); % number of intervals
linenos_interval_start=find(data_numeric.Intervaldata=="Interval data:")+3; % line numbers for interval start
linenos_interval_end=find(data_numeric.Intervaldata=="Interval and data points:")-2; % line numbers for interval end

%for last point
n=size(data_numeric);
linenos_interval_end1=[linenos_interval_end;n(1)];

linenos_interval_end=linenos_interval_end1(2:end);
clear linenos_interval_end1 n

for i=1:interval_number
   structure(i).interval_no=i;
   structure(i).time=data_numeric.Time(linenos_interval_start(i):linenos_interval_end(i)); % in sec
   structure(i).gap=data_numeric.Gap(linenos_interval_start(i):linenos_interval_end(i)); % mm
   structure(i).torque=data_numeric.Torque(linenos_interval_start(i):linenos_interval_end(i)); %mN.m
   structure(i).normal_force=data_numeric.NormalForce(linenos_interval_start(i):linenos_interval_end(i)); %N
   structure(i).defangle=data_numeric.DeflectionAngle(linenos_interval_start(i):linenos_interval_end(i)); %mrad
   structure(i).vel_norm=data_numeric.Velocity(linenos_interval_start(i):linenos_interval_end(i)); %[micrometer/s]
   structure(i).interval_time=data_numeric.IntervalTime(linenos_interval_start(i):linenos_interval_end(i)); %sec
   structure(i).rot_speed=data_numeric.RotationalSpeed(linenos_interval_start(i):linenos_interval_end(i))*2*pi()*1000/60; %converted to mrad/sec
end

end