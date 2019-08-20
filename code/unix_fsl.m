function y = unix_fsl(version, command, varargin)

% Wrapper for calling fsl functions given a particular version

% 2019-08-08: Made compatible with openmind

if optInputs(varargin, 'nohup')
    command = ['nohup ' command ' &'];
end

switch version
    case '5.0.9'
        fsl_dir = '/cm/shared/openmind/fsl/5.0.9';
        full_command = [...
            '. ' fsl_dir '/etc/fslconf/fsl.sh; ' ...
            fsl_dir '/bin/' command];
        fprintf([full_command '\n']);
        [~,y] = unix(full_command);
    case 'openmind-5.0.6'
        error('Needs to be setup');
        fprintf(['/usr/bin/fsl5.0-' command '\n']);
        [~,y] = unix(['/usr/bin/fsl5.0-' command]);
    case '5.0'
        error('Needs to be setup');
        full_command = [root_directory '/fsl/5.0/' command];
        fprintf([full_command '\n']);
        [~,y] = unix(full_command);
        % fprintf(['/usr/bin/fsl5.0-' command '\n']);
        % [~,y] = unix(['/usr/bin/fsl5.0-' command]);
    case '4.1'
        error('Needs to be setup');
        full_command = [root_directory '/fsl/4.1/' command];
        fprintf([full_command '\n']);
        [~,y] = unix(full_command);
        % fprintf(['/usr/bin/fsl4.1-' command '\n']);
        % [~,y] = unix(['/usr/bin/fsl4.1-' command]);
    otherwise
        [~,y] = unix(command);
end