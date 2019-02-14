function raw_data = LTspice2Matlab( filename, varargin )
%
%
% LTspice2Matlab -- Reads an LTspice IV or LTspice XVII .raw waveform file containing data from a Transient
%           Analysis (.tran), AC Analysis (.ac), DC Sweep (.dc), Operating Point (.op), Transfer Function
%           (.tf), FFT (.four), or Noise (.noise) simulation, and converts voltages and currents vs. time
%           (or frequency) into MATLAB variables. It also supports import of stepped simulations.
%           This function can read compressed binary, uncompressed binary, and ASCII file formats. It does
%           not currently support files saved in the Fast Access Format. In the case of compressed binary,
%           the data is automatically uncompressed using fast quadratic point insertion. Note that LTspice
%           uses a lossy compression format (enabled by default) with user adjustable error bounds.
%
%           Use LTspice2Matlab to import LTspice waveforms into MATLAB for additional analysis or to
%           compare with measured data.
%
%           This function has been tested with LTspice IV version 4.01p, and MATLAB versions 6.1 and 7.5.
%           Regression testing has been used to expose the function to a wide range of LTspice settings.
%           Author:     Paul Wagner         2009-04-25
%
%           Not so much testing has been done for the added features (steps, .op, .dc, .tf, and .noise
%           simulations, FFT calculations), but they should work well at least with uncompressed files and
%           without downsampling.
%           Modified:   Peter Feichtinger   2015-04
%
%           Support for reading files written by LTspice XVII (which encodes text in UTF-16) added.
%           Modified:   Peter Feichtinger   2019-02
%
%
%    Calling Convention:
%        RAW_DATA = LTspice2Matlab( FILENAME );                  % Returns all variables found in FILENAME
%                             (or)
%        RAW_DATA = LTspice2Matlab( FILENAME, SELECTED_VARS );   % Returns only selected variables
%               Set SELECTED_VARS to [] to quickly determine the number and names of variables present in
%               FILENAME without actually loading the variables.
%                             (or)
%        RAW_DATA = LTspice2Matlab( FILENAME, SELECTED_VARS, N );
%               Returns variables listed in SELECTED_VARS, with all waveforms downsampled by N. Set N > 1 to
%               load very large data files using less memory, at the price of degraded waveform accuracy and
%               possible aliasing.
%
%    Inputs:  FILENAME is a string containing the name and path of the LTspice .raw file to be converted.
%
%             SELECTED_VARS (optional) is a vector of indexes indicating which variables to extract from the
%                  .raw file. For example, if a .raw file has 14 variables and SELECTED_VARS is [1 8 9], then
%                  the output RAW_DATA.variable_mat will be a 3 x num_data_pnts matrix containing waveforms
%                  for variables 1, 8, and 9 only. Note that SELECTED_VARS does not cover the time (or
%                  frequency) variable (index 0), which is returned separately in RAW_DATA.time_vect (or
%                  RAW_DATA.freq_vect). Extracting only a subset of variables is a way to use less memory
%                  when loading very large simulation files.
%
%                  If this parameter is not specified, then all variables are returned by default. Setting
%                  SELECTED_VARS to 'all' will also cause all variables to be returned.
%
%                  To quickly determine the number and names of variables present in a .raw file, call
%                  LTspice2Matlab with SELECTED_VARS set to []. In this case, all fields in RAW_DATA will be
%                  populated, except .time_vect (or .freq_vect) and .variable_mat, which will both be empty.
%                  Since only the header is read, the function call should execute very quickly, even for
%                  large files.
%
%             N (optional) must be a positive integer >= 1. If N is specified, then SELECTED_VARS must also
%                  be specified. If N is unspecified, it defaults to 1, which does not change the sampling
%                  rate. If this value is 2 or larger, the returned voltage, current, and time data will be
%                  downsampled by keeping every N-th sample in the original data, starting with the first.
%                  Caution: No lowpass filtering is applied prior to downsampling, so aliasing may occur.
%                  Also, in many cases LTspice saves data with a non-constant sampling rate, in which case
%                  downsampling can result in substantial waveform distortion. This option should only be
%                  used if the waveform of interest is initially oversampled.
%
%    Outputs:  RAW_DATA is a Matlab structure containing the following fields:
%                  title:                String containing the title appearing in the .raw file header.
%                  date:                 String containing the date appearing in the .raw file header.
%                  plotname:             String indicating simulation type (e.g. 'Transient Analysis')
%                  conversion_notes:     Description of modifications (if any) done to the data during
%                                        conversion.
%                  num_variables:        Number of variables (does not include the time/frequency variable).
%                  variable_type_list:   A cell of strings indicating the variable type (e.g. 'voltage').
%                  variable_name_list:   A cell of strings indicating the name of each variable.
%                  selected_vars:        A vector of indicies referencing variable_type_list cells,
%                                        corresponding to each row in variable_mat.
%                  num_data_pnts:        Number of data points for each variable.
%                  num_steps:            The number of steps for stepped simulations.
%                  variable_mat:         Double precision matrix with num_variables rows and num_data_pnts
%                                        columns. This matrix contains node voltages (in Volts) and device
%                                        currents (in Amps) for each variable and each time point listed in
%                                        time_vect (or freq_vect).
%                                        For AC Analysis simulations, variable_mat will have complex values
%                                        showing the real and imaginary components of the voltage or current
%                                        at the corresponding frequency. To convert this to log magnitude and
%                                        normalized phase representation used in LTspice plots, use the
%                                        following formulas:
%                                            Log_Magnitude_dB    =  20*log10(abs(variable_mat))
%                                            Norm_Phase_degrees  =  angle(variable_mat) * 180/pi
%                  time_vect:            [Field returned for Transient Analysis only] Double precision row
%                                        vector of time values (in seconds) at each simulation point.
%                    (or)
%                  freq_vect:            [Field returned for AC Analysis only] Double precision row vector of
%                                        frequency values (in Hz) at each simulation point.
%                    (or)
%                  source_vect:          [Field returned for DC Sweeps only] Double precision row vector of
%                                        source values (in Amperes or Volts) at each simulation point.
%                    (or)
%                  param_vect:           [Field returned for stepped DC Operating Point or Transfer Function
%                                        Analysis only] Double precision row vector of values for the stepped
%                                        parameter.
%                  source_name:          [Field returned for DC Sweeps only] The name of the source that is
%                                        swept.
%                    (or)
%                  param_name:           [Field returned for stepped DC Operating Point or Transfer Function
%                                        Analysis only] The name of the stepped parameter.
%
%                  ** For stepped simulations that have been properly recognized as such, an additional
%                     dimension is added to variable_mat and the time, frequency, source or param vector,
%                     respectively.
%
%    ** Currently this function is able to import results from Transient Analysis (.tran), AC Analysis (.ac),
%       DC Sweeps (.dc), Operating Point Analysis (.op), Transfer Function (.tf) and Noise (.noise)
%       simulations, and FFT calculations (.four).
%
%
%    Examples
%    --------
%    These examples assume you've run a .tran simulation in LTspice for a hypothetical file called
%    BASIC_CIRCUIT.ASC, and that an output file called BASIC_CIRCUIT.RAW has been created. It also
%    assumes your current MATLAB directory is pointing to the directory where the .raw file is
%    located (or that you prepended the full path to the input parameter FILENAME).
%
%    To import BASIC_CIRCUIT.RAW into MATLAB and create a labeled plot of a single variable vs. time:
%
%       raw_data = LTspice2Matlab('BASIC_CIRCUIT.RAW');
%       variable_to_plot = 1;   % This example plots the first variable in the data structure.
%       plot( raw_data.time_vect, raw_data.variable_mat(variable_to_plot,:), 'k' );
%       title( sprintf('Waveform %s', raw_data.variable_name_list{variable_to_plot}) );
%       ylabel( raw_data.variable_type_list{variable_to_plot} );
%       xlabel( 'Time (sec)' );
%
%    To superimpose all variables in BASIC_CIRCUIT.RAW on a single graph with a legend:
%
%       raw_data = LTspice2Matlab('BASIC_CIRCUIT.RAW');
%       plot( raw_data.time_vect, raw_data.variable_mat );
%       title( sprintf( 'File:  %s', raw_data.title) );
%       legend( raw_data.variable_name_list );
%       ylabel( 'Voltage (V) or Current (A)' );
%       xlabel( 'Time (sec)' );
%
%     To quickly determine the number and names of variables in BASIC_CIRCUIT.RAW without loading the entire
%     file:
%
%       raw_data = LTspice2Matlab('BASIC_CIRCUIT.RAW', []);
%       disp( sprintf('\n\nThis file contains %.0f variables:\n', raw_data.num_variables) );
%       disp( sprintf('NAME         TYPE\n-------------------------') );
%       disp( [char(raw_data.variable_name_list), char(zeros(raw_data.num_variables,5)), ...
%              char(raw_data.variable_type_list)] );
%

    % Initialize the output structure.
    raw_data = [];

    % Process function arguments
    if nargin == 1
        selected_vars = 'all';
        downsamp_N = 1;
    elseif nargin == 2
        selected_vars = varargin{1};
        if ischar(selected_vars),  selected_vars = lower(selected_vars);  end
        downsamp_N = 1;
    elseif nargin == 3
        selected_vars = varargin{1};
        if ischar(selected_vars),  selected_vars = lower(selected_vars);  end
        downsamp_N = varargin{2};
    else
        error( 'LTspice2Matlab takes only 1, 2, or 3 input parameters.  Type "help LTspice2Matlab" for details' );
    end

    if length(downsamp_N) ~= 1 || ~isnumeric(downsamp_N) || isnan(downsamp_N) || mod(downsamp_N,1) ~= 0.0 || downsamp_N <= 0
        error( 'Optional parameter DOWNSAMP_N must be a positive integer >= 1' );
    end

    % Try to open file
    filename = strtrim(filename);   % Remove leading and trailing spaces from filename.
    fid = fopen(filename, 'rb');
    if length(fid) == 1 && isnumeric(fid) && fid == -1
        % try to append ".raw" to the file name ...
        fid = fopen(sprintf( '%s.raw', filename ), 'rb');
        if length(fid) == 1 && isnumeric(fid) && fid == -1
            error( 'Could not open file "%s"', filename );
        end
    end

    [filename, ~, machineformat] = fopen(fid);

    % Detect text encoding: LTspice IV uses ASCII, LTspice XVII uses UTF-16
    utf16 = false;
    [buf, count] = fread(fid, 2, '*uint8');
    if count == 2 && any(buf == 0)
        utf16 = true;
    end
    clear buf count;
    frewind(fid);

    % Load header tags & information
    % Variables include voltages and currents only.  Does not include the time vector.
    variable_name_list = {};
    variable_type_list = {};
    variable_flag = 0;
    file_format = '';
    while 1
        if ~utf16
            the_line = fgetl(fid);
        else
            % fopen doesn't know UTF-16, read manually using fread instead.
            the_line = [];
            while ~feof(fid)
                [buf, count] = fread(fid, 100, 'uint16=>char', 0, 'ieee-le');
                nl = find(buf == newline, 1);
                if isempty(nl)
                    the_line = [the_line buf']; %#ok<AGROW>
                    continue;
                end

                the_line = [the_line buf(1:nl(1)-1)']; %#ok<AGROW>
                fseek(fid, (-count + nl(1)) * 2, 0);
                break;
            end
            if isempty(the_line) && feof(fid)
                % Do as fgetl does
                the_line = -1;
            end
        end
        if ~ischar(the_line)
            try fclose( fid );  catch, end
            error( 'Format error in LTspice file "%s" ... Unexpected end of file', filename );
        end

        if contains(the_line, 'Binary:')
            file_format = 'binary';
            break;
        elseif contains(the_line, 'Values:')
            file_format = 'ascii';
            break;
        end

        if variable_flag == 0  % Non-variable header section
            if isempty(the_line)
                colon_index = [];
            else
                colon_index = find( the_line == ':' );
            end
            if isempty(colon_index)
                try fclose( fid );  catch, end
                error( 'Format error in LTspice file "%s"', filename );
            end
            var_name = the_line(1:(colon_index(1)-1));
            var_value = strtrim(the_line((colon_index(1)+1):end));

            vn_keep = var_name( var_name ~= ' ' & var_name ~= '.' & var_name ~= char(9) & var_name ~= newline & var_name ~= char(13) );
            var_name = lower(vn_keep);
            if isempty(var_name) || (var_name(1)>='0' && var_name(1)<='9')
                try fclose( fid );  catch, end
                error( 'Format error in LTspice file "%s" ... Bad tag name found', filename );
            end

            if strcmpi( var_name, 'variables' ) || strcmpi( var_name, 'variable' )
                variable_flag = 1;
                continue;
            end
            value_try = str2double(var_value);
            try
                if isnan(value_try)
                    raw_data.(var_name) = var_value;
                else
                    raw_data.(var_name) = value_try;
                end
            catch
                try fclose( fid );  catch, end
                error( 'Format error in LTspice file "%s" ... Bad tag name found', filename );
            end

        else  % Variable header section
            leading_ch_index = find( (the_line(1:end-1) == ' ' | the_line(1:end-1) == char(9)) & (the_line(2:end) ~= ' ' & the_line(2:end) ~= char(9)) );
            if length(leading_ch_index) ~= 3
                try fclose( fid );  catch, end
                error( 'Format error in LTspice file "%s" ... Wrong number of columns in the variable define section', filename );
            end

            part1 = strtrim(the_line( (leading_ch_index(1)+1) : leading_ch_index(2) ));
            part2 = strtrim(the_line( (leading_ch_index(2)+1) : leading_ch_index(3) ));
            part3 = strtrim(the_line( (leading_ch_index(3)+1) : end ));

            if str2double(part1) ~= length(variable_name_list)
                try fclose( fid );  catch, end
                error( 'Format error in LTspice file "%s" ... Inconsistency found in the variable define section', filename );
            end
            variable_name_list{end+1} = part2; %#ok<AGROW>
            variable_type_list{end+1} = part3; %#ok<AGROW>
        end
    end

    % Check raw_data structure for required fields
    expected_tags      = {'title', 'date', 'plotname', 'flags', 'novariables',   'nopoints'    };
    expected_tags_full = {'Title', 'Date', 'Plotname', 'Flags', 'No. Variables', 'No. Points'  };
    for q=1:length(expected_tags)
        if ~isfield( raw_data, lower(expected_tags{q}) )
            try fclose( fid );  catch, end
            error( 'Format error in LTspice file "%s" ... tag "%s" not found', filename, expected_tags_full{q} );
        end
    end

    % Rename number fields
    raw_data.conversion_notes = '';
    raw_data.num_data_pnts  = raw_data.nopoints;         raw_data = rmfield( raw_data, 'nopoints' );
    raw_data.num_variables  = raw_data.novariables - 1;  raw_data = rmfield( raw_data, 'novariables' );
    % "raw_data.num_variables" does not include the time vector (index 0 in the .raw file)

    % Check for stepped simulation
    if contains(raw_data.flags, 'stepped', 'IgnoreCase',true)
        step_flag = true;
    else
        step_flag = false;
    end
    raw_data.num_steps = 1;

    % Remove unneeded fields and store offset
    if isfield( raw_data, 'command' ),          raw_data = rmfield( raw_data, 'command' );  end
    if isfield( raw_data, 'backannotation' ),   raw_data = rmfield( raw_data, 'backannotation' );  end
    if isfield( raw_data, 'offset' )
        general_offset = raw_data.offset;  %(sec)
        raw_data = rmfield( raw_data, 'offset' );
    else
        general_offset = 0.0;
    end

    % cut off the time variable.
    raw_data.variable_name_list = variable_name_list(2:end);
    raw_data.variable_type_list = variable_type_list(2:end);

    % determine simulation type
    if     contains(raw_data.plotname, 'transient analysis', 'IgnoreCase',true),           simulation_type = '.tran';
    elseif contains(raw_data.plotname, 'ac analysis', 'IgnoreCase',true),                  simulation_type = '.ac';
    elseif contains(raw_data.plotname, 'dc transfer characteristic', 'IgnoreCase',true),   simulation_type = '.dc';
    elseif contains(raw_data.plotname, 'operating point', 'IgnoreCase',true),              simulation_type = '.op';
    elseif contains(raw_data.plotname, 'transfer function', 'IgnoreCase',true),            simulation_type = '.tf';
    elseif contains(raw_data.plotname, 'fft of time domain data', 'IgnoreCase',true),      simulation_type = '.four';
    elseif contains(raw_data.plotname, 'noise spectral density', 'IgnoreCase',true),       simulation_type = '.noise';
    else
        try fclose( fid );  catch, end
        error( 'Currently LTspice2Matlab is only able to import results from Transient Analysis (.tran), AC Analysis (.ac), Operating Point (.op), Transfer Function (.tf), DC Sweep (.dc), Noise (.noise) simulations and FFT (.four) calculations.' );
    end

    % Check for the expected formats for every simulation type
    if contains(raw_data.flags, 'fastaccess', 'IgnoreCase',true)
        try fclose( fid );  catch, end
        error( 'LTspice2Matlab cannot convert files saved in the "Fast Access" format.' );
    end
    switch simulation_type
        case '.tran'
            if ~contains(raw_data.flags, 'real', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "real" flag for a Transient Analysis (.tran) simulation.  Unsure how to convert the data' );
            end
            if ~contains(raw_data.flags, 'forward', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "forward" flag for a Transient Analysis (.tran) simulation.  Unsure how to convert the data' );
            end

        case '.ac'
            if ~contains(raw_data.flags, 'complex', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "complex" flag for an AC Analysis (.ac) simulation.  Unsure how to convert the data' );
            end
            if ~contains(raw_data.flags, 'forward', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "forward" flag for an AC Analysis (.ac) simulation.  Unsure how to convert the data' );
            end

        case '.dc'
            if ~contains(raw_data.flags, 'real', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "real" flag for a DC Sweep (.dc) simulation.  Unsure how to convert the data' );
            end
            if ~contains(raw_data.flags, 'forward', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "forward" flag for a DC Sweep (.dc) simulation.  Unsure how to convert the data' );
            end

        case '.op'
            if ~contains(raw_data.flags, 'real', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "real" flag for an Operating Point (.op) simulation.  Unsure how to convert the data' );
            end

        case '.tf'
            if ~contains(raw_data.flags, 'real', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "real" flag for a Transfer function (.tf) simulation.  Unsure how to convert the data' );
            end

        case '.four'
            if ~contains(raw_data.flags, 'complex', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "complex" flag for an FFT calculation.  Unsure how to convert the data' );
            end
            if ~contains(raw_data.flags, 'forward', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "forward" flag for an FFT calculation.  Unsure how to convert the data' );
            end

        case '.noise'
            if ~contains(raw_data.flags, 'real', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "real" flag for a Noise simulation.  Unsure how to convert the data' );
            end
            if ~contains(raw_data.flags, 'forward', 'IgnoreCase',true)
                try fclose( fid );  catch, end
                error( 'Expected to find "forward" flag for a Noise simulation.  Unsure how to convert the data' );
            end
    end

    % Remove flags field
    %if isfield( raw_data, 'flags' ),  raw_data = rmfield( raw_data, 'flags' );  end

    % Look if selected_vars is a string that says 'all' or similar
    if ischar(selected_vars)
        if strcmpi(selected_vars, 'all') || strcmpi(selected_vars, 'everything') || strcmpi(selected_vars, 'complete') || ...
                strcmpi(selected_vars, 'all variables') || strcmpi(selected_vars, 'all vars') || ...
                strcmpi(selected_vars, 'every thing') || strcmpi(selected_vars, 'every')
            selected_vars = 1:raw_data.num_variables;     %Return all variables
        else
            try fclose( fid );  catch, end
            error( 'Bad value for optional input parameter SELECTED_VARS' );
        end
    end

    % selected_vars is empty, only return structure without data
    if size(selected_vars,1) == 0 || size(selected_vars,2) == 0
        raw_data.selected_vars = [];
        raw_data.variable_mat = [];
        raw_data.time_vect = [];
        try fclose( fid );  catch, end
        return;
    end

    % Sanity check selected_vars
    if size(selected_vars,1) > 1 && size(selected_vars,2) > 1
        try fclose( fid );  catch, end
        error( 'SELECTED_VARS must be a row or column vector, not a matrix' );
    end
    if ~isempty(find(selected_vars == 0, 1))
        try fclose( fid );  catch, end
        error( 'The time vector (index 0) is returned separately.\n   Values in input parameter SELECTED_VARS must be positive integers >= 1 and <= NUM_VARIABLES' );
    end
    non_integer_index = find(isnan(selected_vars) | ~isnumeric(selected_vars) | mod( selected_vars, 1 ) ~= 0.0, 1);
    if ~isempty(non_integer_index)
        try fclose( fid );  catch, end
        error( 'Values in input parameter SELECTED_VARS must be positive integers >= 1 and <= NUM_VARIABLES' );
    end
    missing_index = find( ~ismember( selected_vars, 1:raw_data.num_variables ), 1 );
    if ~isempty(missing_index)
        try fclose( fid );  catch, end
        error( 'Error in input parameter SELECTED_VARS ... Out of range value(s) found' );
    end

    selected_vars = unique(selected_vars);   % remove duplicates and sort in ascending order.
    raw_data.selected_vars = selected_vars;

    % Apply donwsampling
    NumPnts = raw_data.num_data_pnts;
    NumPnts_DS = floor(NumPnts / downsamp_N);
    raw_data.num_data_pnts = NumPnts_DS;   % Updated # of points
    NumVars = raw_data.num_variables + 1;


    % READ IN THE ACTUAL WAVEFORM DATA
    if strcmpi(file_format, 'binary')
        binary_start = ftell(fid);      % start of binary data section.

        if strcmpi( simulation_type, '.tran' ) || strcmpi( simulation_type, '.dc' ) || strcmpi( simulation_type, '.noise' )
            % For Transient Analysis simulations, the time data is stored in double precision floating point binary format,
            % and everything else is stored in single precision format. For stepped simulations, a new step is intoduced
            % by a duplicate start entry.
            % For DC Sweep simulations, the source value is stored in double precision format (8 bytes) and the
            % operating point values in real format. A new step is introduced by a duplicate start entry.
            % For Noise simulations, the frequency data is stored in double precision floating point binary format,
            % and gain, output noise and input referred noise are stored in single precision format. It is unknown,
            % how new steps are introduced.

            % Extract the binary data in the fewest possible number of contiguous blocks
            if length(selected_vars) > 1
                g_border = find( [2, diff(selected_vars), 2] ~= 1 );
                block_list = cell(1, length(g_border) - 1);
                for k=1:length(g_border)-1,  block_list{k} = g_border(k):(g_border(k+1)-1);  end
            else
                block_list = {1:length(selected_vars)};
            end

            raw_data.variable_mat = zeros(length(selected_vars), NumPnts_DS);  % Initialize.
            for k=1:length(block_list)
                target_var_index = selected_vars(block_list{k});
                fseek(fid, binary_start + (target_var_index(1)+1)*4, 'bof');
                TVIL = length(target_var_index);
                bytes_skip = (NumVars+1-TVIL)*4  +  (downsamp_N-1)*(NumVars+1)*4;
                precision_str = sprintf('%.0f*float',TVIL);
                raw_data.variable_mat(block_list{k},:) = reshape( fread(fid, NumPnts_DS*TVIL, precision_str, bytes_skip, machineformat), TVIL, NumPnts_DS );
            end

            fseek(fid, binary_start, 'bof');    % rewind to start, then extract the time/source vector.
            raw_data.time_vect = fread( fid, NumPnts_DS, 'double', (NumVars-1)*4 + (downsamp_N-1)*(NumVars+1)*4, machineformat ).';

            % Process stepped data
            if step_flag
                if downsamp_N > 1
                    warning( 'LTspice2Matlab:downsampleSteps', ...
                            'Stepped data found and downsampling enabled. Steps may not be recognized properly.' );
                end

                % Stepped data starts on a duplicate entry in LTspice IV
                steps = find( diff(raw_data.time_vect) == 0.0 );
                if numel(steps) > 0
                    % Remove duplicate entries
                    raw_data.variable_mat(:,steps) = [];
                    raw_data.time_vect(steps) = [];

                    num_steps = numel(steps) + 1;
                else
                    % Stepped data starts with the same time/source value in LTspice XVII
                    num_steps = numel(find( raw_data.time_vect == raw_data.time_vect(1) ));
                end

                raw_data.num_steps = num_steps;
                % Reshape value matrix and time/source vector
                mat_size = size(raw_data.variable_mat);
                raw_data.variable_mat = reshape(raw_data.variable_mat, mat_size(1), mat_size(2) / num_steps, num_steps);
                raw_data.time_vect = reshape(raw_data.time_vect, mat_size(2) / num_steps, num_steps).';
                % Update num_data_pnts
                raw_data.num_data_pnts = mat_size(2) / num_steps;
            end

            % Rename time_vect and add source name for .dc simulations
            if strcmpi( simulation_type, '.dc' )
               raw_data.source_vect = raw_data.time_vect;
               raw_data = rmfield( raw_data, 'time_vect' );
            % Rename time_vect for .noise simulations
            elseif strcmpi( simulation_type, '.noise' )
               raw_data.freq_vect = raw_data.time_vect;
               raw_data = rmfield( raw_data, 'time_vect' );
            end

        elseif strcmpi( simulation_type, '.ac' ) || strcmpi( simulation_type, '.four' )
            % For AC Analysis simulations, the frequency data is stored in double precision format (8 + 8 unused bytes),
            % and the variables are stored as complex double precision arrays (8 bytes real followed by 8 bytes imag).
            % For stepped simulations, every step is simply appended.
            % The same goes for FFT calculations.

            % Extract the binary data in the fewest possible number of contiguous blocks
            if length(selected_vars) > 1
                g_border = find( [2, diff(selected_vars), 2] ~= 1 );
                block_list = cell(1, length(g_border) - 1);
                for k=1:length(g_border)-1,  block_list{k} = g_border(k):(g_border(k+1)-1);  end
            else
                block_list = {1:length(selected_vars)};
            end

            raw_data.variable_mat = zeros(length(selected_vars), NumPnts_DS);  % Initialize.
            if numel(raw_data.variable_mat) ~= 0,  raw_data.variable_mat(1,1) = 0.0 + 1i*0.0;  end  % Allocate memory for complex double.
            for k=1:length(block_list)
                target_var_index = selected_vars(block_list{k});
                fseek(fid, binary_start + target_var_index(1)*16, 'bof');
                TVIL = length(target_var_index);
                bytes_skip = (NumVars-TVIL)*16  +  (downsamp_N-1)*NumVars*16;
                precision_str = sprintf('%.0f*double',TVIL*2);
                temp_buff = reshape(fread(fid, NumPnts_DS*TVIL*2, precision_str, bytes_skip, machineformat), TVIL*2, NumPnts_DS );
                raw_data.variable_mat(block_list{k},:) = temp_buff(1:2:end-1,:) + 1i*temp_buff(2:2:end,:);
                clear temp_buff;
            end

            fseek(fid, binary_start, 'bof');    % rewind to start, then extract the time vector.
            raw_data.freq_vect = fread( fid, NumPnts_DS, 'double', (NumVars-1)*16 + 8 + (downsamp_N-1)*NumVars*16, machineformat ).';

            % Process stepped data
            if step_flag
                % Stepped data starts with low frequency again, should work with downsampled data
                steps = find( diff(raw_data.freq_vect) < 0.0 ) + 1;
                num_steps = numel(steps) + 1;
                raw_data.num_steps = num_steps;
                % Reshape value matrix and time/source vector
                mat_size = size(raw_data.variable_mat);
                raw_data.variable_mat = reshape(raw_data.variable_mat, mat_size(1), mat_size(2) / num_steps, num_steps);
                raw_data.freq_vect = reshape(raw_data.freq_vect, mat_size(2) / num_steps, num_steps).';
                % Update num_data_pnts
                raw_data.num_data_pnts = mat_size(2) / num_steps;
            end

        elseif strcmpi( simulation_type, '.op' ) || strcmpi( simulation_type, '.tf' )
            % For simple Operating Point and Transfer Function simulations, the first variable is stored in double
            % precision format (8 bytes), the remaining variables are stored in single precision format.
            % For stepped Operating point and Transfer Function simulations, the stepping parameter is stored in double
            % precision format (8 bytes) and the variables are stored in single precision format.

            if step_flag
                % Extract the binary data in the fewest possible number of contiguous blocks
                if length(selected_vars) > 1
                    g_border = find( [2, diff(selected_vars), 2] ~= 1 );
                    block_list = cell(1, length(g_border) - 1);
                    for k=1:length(g_border)-1,  block_list{k} = g_border(k):(g_border(k+1)-1);  end
                else
                    block_list = {1:length(selected_vars)};
                end

                raw_data.variable_mat = zeros(length(selected_vars), NumPnts_DS);  % Initialize.
                for k=1:length(block_list)
                    target_var_index = selected_vars(block_list{k});
                    fseek(fid, binary_start + (target_var_index(1)+1)*4, 'bof');
                    TVIL = length(target_var_index);
                    bytes_skip = (NumVars+1-TVIL)*4  +  (downsamp_N-1)*(NumVars+1)*4;
                    precision_str = sprintf('%.0f*float',TVIL);
                    raw_data.variable_mat(block_list{k},:) = reshape( fread(fid, NumPnts_DS*TVIL, precision_str, bytes_skip, machineformat), TVIL, NumPnts_DS );
                end

                fseek(fid, binary_start, 'bof');    % rewind to start, then extract the step parameter vector.
                raw_data.param_vect = fread( fid, NumPnts_DS, 'double', (NumVars-1)*4 + (downsamp_N-1)*(NumVars+1)*4, machineformat ).';
                raw_data.param_name = variable_name_list{1};
                raw_data.num_steps = raw_data.num_data_pnts;
                raw_data.num_data_pnts = 1;

            else
                % Ignore selected_vars and read everything, since we have only one set of variables
                raw_data.variable_mat = zeros(1, NumVars);
                fseek(fid, binary_start, 'bof');
                raw_data.variable_mat(1) = fread(fid, 1, 'double', 0, machineformat);
                raw_data.variable_mat(2:end) = fread(fid, NumVars - 1, sprintf('%.0f*float', NumVars - 1), 0, machineformat);
                % First variable is normal variable, not time or step parameter
                raw_data.selected_vars = 0:(NumVars - 1);
                raw_data.variable_name_list = variable_name_list;
                raw_data.variable_type_list = variable_type_list;
                raw_data.num_variables = NumVars;

            end

        else
            try fclose( fid );  catch, end
            error( 'Simulation type (%s) not currently supported', simulation_type );
        end

        if downsamp_N == 1
            raw_data.conversion_notes = 'Converted from Binary format';
        else
            raw_data.conversion_notes = sprintf( 'Converted from Binary format.  Downsampled from %.0f to %.0f points', NumPnts, NumPnts_DS );
        end


    elseif strcmpi(file_format, 'ascii' )

        if strcmpi( simulation_type, '.tran' ) || strcmpi( simulation_type, '.dc')
            % Format:  point number, time value, var1, var2 ... varN
            raw_data.variable_mat = fscanf( fid, '%g', [raw_data.num_variables+2, raw_data.num_data_pnts] );   %matrix is filled in column order.
            if (size(raw_data.variable_mat,1) ~= raw_data.num_variables+2) || (size(raw_data.variable_mat,2) ~= raw_data.num_data_pnts)
                error( 'Format error in ASCII Transient Analysis or DC Sweep LTspice file "%s" ... Incorrect number of data values read', filename );
            end
            raw_data.time_vect = raw_data.variable_mat(2,1:downsamp_N:end);
            raw_data.variable_mat = raw_data.variable_mat(2+selected_vars,1:downsamp_N:end);

            % Process stepped data
            if step_flag
                if downsamp_N > 1
                   warning( 'LTspice2Matlab:downsampleSteps', 'Stepped data found and downsampling enabled. Steps may not be recognized properly.' );
                end

                % Stepped data starts with same start time/source value
                num_steps = numel(find( raw_data.time_vect == raw_data.time_vect(1) ));
                raw_data.num_steps = num_steps;
                % Reshape value matrix and time/source vector
                mat_size = size(raw_data.variable_mat);
                raw_data.variable_mat = reshape(raw_data.variable_mat, mat_size(1), mat_size(2) / num_steps, num_steps);
                raw_data.time_vect = reshape(raw_data.time_vect, mat_size(2) / num_steps, num_steps).';
                % Update num_data_pnts
                raw_data.num_data_pnts = mat_size(2) / num_steps;
            end

            % Rename time_vect and add source name for .dc simulations
            if strcmpi( simulation_type, '.dc' )
               raw_data.source_vect = raw_data.time_vect;
               raw_data = rmfield( raw_data, 'time_vect' );
               raw_data.source_name = variable_name_list{1};
            end

        elseif strcmpi( simulation_type, '.ac' ) || strcmpi( simulation_type, '.four' )
            % Format:  point number, freq value, 0, var1 real, var1 imag, var2 real, var2 imag ... varN real, varN imag
            all_data = fread( fid, inf, 'uchar' );
            all_data( all_data == ',' ) = sprintf( '\t' );  %Replace commas with tab characters
            raw_data.variable_mat = sscanf( char(all_data), '%g', [3+2*raw_data.num_variables, raw_data.num_data_pnts] );
            clear all_data;
            % raw_data.variable_mat = fscanf( fid, '%g', [3+2*raw_data.num_variables, raw_data.num_data_pnts] );   %matrix is filled in column order.

            if (size(raw_data.variable_mat,1) ~= (3+2*raw_data.num_variables)) || (size(raw_data.variable_mat,2) ~= raw_data.num_data_pnts)
                error( 'Format error in ASCII AC Analysis LTspice file "%s" ... Incorrect number of data values read', filename );
            end
            raw_data.freq_vect = raw_data.variable_mat(2,1:downsamp_N:end);
            raw_data.variable_mat = raw_data.variable_mat(3+selected_vars*2-1,1:downsamp_N:end) + 1i*raw_data.variable_mat(3+selected_vars*2,1:downsamp_N:end);

            % Process stepped data
            if step_flag
                % Stepped data starts with low frequency again, should work with downsampled data
                steps = find( diff(raw_data.freq_vect) < 0.0 ) + 1;
                num_steps = numel(steps) + 1;
                raw_data.num_steps = num_steps;
                % Reshape value matrix and time/source vector
                mat_size = size(raw_data.variable_mat);
                raw_data.variable_mat = reshape(raw_data.variable_mat, mat_size(1), mat_size(2) / num_steps, num_steps);
                raw_data.freq_vect = reshape(raw_data.freq_vect, mat_size(2) / num_steps, num_steps).';
                % Update num_data_pnts
                raw_data.num_data_pnts = mat_size(2) / num_steps;
            end

        elseif strcmpi( simulation_type, '.op' ) || strcmpi( simulation_type, '.tf' )
            % Format: point number, first value, var2, var3 ... varN
            % Format stepped: point number, step parameter, var1, var2 ... varN

            if step_flag
                raw_data.variable_mat = fscanf( fid, '%g', [raw_data.num_variables+2, raw_data.num_data_pnts] );   %matrix is filled in column order.
                if (size(raw_data.variable_mat,1) ~= raw_data.num_variables+2) || (size(raw_data.variable_mat,2) ~= raw_data.num_data_pnts)
                    error( 'Format error in ASCII Operating Point or Transfer Function LTspice file "%s" ... Incorrect number of data values read', filename );
                end
                raw_data.param_vect = raw_data.variable_mat(2,1:downsamp_N:end);
                raw_data.variable_mat = raw_data.variable_mat(2+selected_vars,1:downsamp_N:end);
                raw_data.param_name = variable_name_list{1};
                raw_data.num_steps = raw_data.num_data_pnts;
                raw_data.num_data_pnts = 1;

            else
                raw_data.variable_mat = fscanf( fid, '%g', [raw_data.num_variables+2, raw_data.num_data_pnts] );   %matrix is filled in column order.
                if (size(raw_data.variable_mat,1) ~= raw_data.num_variables+2) || (size(raw_data.variable_mat,2) ~= raw_data.num_data_pnts)
                    error( 'Format error in ASCII Operating Point or Transfer Function LTspice file "%s" ... Incorrect number of data values read', filename );
                end
                raw_data.variable_mat = raw_data.variable_mat(2:end,1:downsamp_N:end).';
                % First variable is normal variable, not time or step parameter
                raw_data.selected_vars = 0:(NumVars - 1);
                raw_data.variable_name_list = variable_name_list;
                raw_data.variable_type_list = variable_type_list;
                raw_data.num_variables = NumVars;

            end

        else
            try fclose( fid );  catch, end
            error( 'Simulation type (%s) not currently supported', simulation_type );
        end

        if downsamp_N == 1
            raw_data.conversion_notes = 'Converted from ASCII format';
        else
            raw_data.conversion_notes = sprintf( 'Converted from ASCII format. Downsampled from %.0f to %.0f points', NumPnts, NumPnts_DS );
        end

    else
        try fclose( fid );  catch, end
        error( 'Format error in LTspice file "%s" ... Data type ID tag not found', filename );
    end

    try fclose( fid );  catch, end


    % Deal with potential compression in Transient Analysis simulations
    if strcmpi( simulation_type, '.tran' ) && ...       % Check to see if the time vector is monotonically increasing.
            (isvector(raw_data.time_vect) && (min(diff(raw_data.time_vect)) < 0.0) || ...
            ~isvector(raw_data.time_vect) && (min(min(diff(raw_data.time_vect, 1, 2))) < 0.0))

        if downsamp_N ~= 1  % If we have already downsampled then we can't uncompress.
            raw_data.time_vect = abs(raw_data.time_vect);

        else
            % The binary file contains 2nd order compression ... use 2nd-order interpolation to add data points in the vicinity of negative time points
            if step_flag
               warning( 'LTspice2Matlab:stepCompressed', 'Stepped simulations in combination with compression might not be read correctly. Consider not using compression.' );
            end

            t_vect = raw_data.time_vect;   % We will add in the offset later.
            neg_pnt_index = find( t_vect < 0.0 & [0,ones(1,length(t_vect)-1)] );
            t_vect = abs(t_vect);

            x1 = t_vect(neg_pnt_index-1);  x2 = t_vect(neg_pnt_index);  x3 = t_vect(neg_pnt_index+1);
            x_new = [(2*x1 + x2)/3;  (x1 + 2*x2)/3;  (2*x2 + x3)/3;  (x2 + 2*x3)/3];   %New sample points

            t_vect_big = NaN*zeros(6,length(t_vect));
            t_vect_big(1,:) = t_vect;
            t_vect_big(4,neg_pnt_index) = t_vect(neg_pnt_index);
            t_vect_big(1,neg_pnt_index) = NaN;
            t_vect_big([2 3 5 6],neg_pnt_index) = x_new;

            full_index = find(~isnan(t_vect_big));
            time_vect_new = t_vect_big(full_index).';   % This is the new time vector with the inserted points.
            t_vect_big([1,4],:) = NaN;
            nan_vect = isnan(t_vect_big(full_index));
            new_index = find( ~nan_vect );  % Index into time_vect_new indicating the new points only.
            old_index = find( nan_vect );
            clear t_vect t_vect_big full_index nan_vect;

            x1sqr = repmat( x1.^2, [4,1] );  x2sqr = repmat( x2.^2, [4,1] );  x3sqr = repmat( x3.^2, [4,1] );
            x1 = repmat( x1, [4,1] );  x2 = repmat( x2, [4,1] );  x3 = repmat( x3, [4,1] );
            denom = (x1sqr-x2sqr).*(x2-x3) - (x2sqr-x3sqr).*(x1-x2);
            r1 = (x_new.^2 - x1sqr)./denom;
            r2 = (x_new - x1)./denom;
            p1 = (x2-x3).*r1 - (x2sqr-x3sqr).*r2;
            p3 = (x1-x2).*r1 - (x1sqr-x2sqr).*r2;
            p2 = -p1 - p3;
            p1 = p1 + 1;
            clear x_new x1sqr x2sqr x3sqr x1 x2 x3 denom r1 r2;

            raw_data.variable_mat(:,end+1:length(time_vect_new)) = 0.0;   % Init the memory
            for k=1:size(raw_data.variable_mat,1)
                y_vect = raw_data.variable_mat(k,1:length(raw_data.time_vect));
                raw_data.variable_mat(k,old_index) = y_vect;
                y_new = repmat(y_vect(neg_pnt_index-1),[4,1]).*p1  +  repmat(y_vect(neg_pnt_index),[4,1]).*p2  +  repmat(y_vect(neg_pnt_index+1),[4,1]).*p3;
                raw_data.variable_mat(k,new_index) = y_new(:).';
            end
            raw_data.time_vect = time_vect_new;
            clear time_vect_new y_vect y_new new_index old_index neg_pnt_index p1 p2 p3;

            raw_data.conversion_notes = sprintf( 'Converted from Binary format with 2nd Order compression.  Upsampled waveforms from %.0f to %.0f points', ...
                raw_data.num_data_pnts, length(raw_data.time_vect) );
            raw_data.num_data_pnts = length(raw_data.time_vect);
        end

    end

    if isfield( raw_data, 'time_vect' )
        raw_data.time_vect = raw_data.time_vect + general_offset;
    elseif isfield( raw_data, 'freq_vect' )
        raw_data.freq_vect = raw_data.freq_vect + general_offset;
    elseif isfield( raw_data, 'source_vect' )
        raw_data.source_vect = raw_data.source_vect + general_offset;
    elseif isfield( raw_data, 'param_vect' )
        raw_data.param_vect = raw_data.param_vect + general_offset;
    end
