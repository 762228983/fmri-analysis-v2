function y = unix_fsl(version, command, varargin)

if optInputs(varargin, 'nohup')
    command = ['nohup ' command ' &'];
end

switch version
    case '5.0'
        fprintf(['/usr/bin/fsl5.0-' command '\n']);
        [~,y] = unix(['/usr/bin/fsl5.0-' command]);
    case '4.1'
        fprintf(['/usr/bin/fsl4.1-' command '\n']);
        [~,y] = unix(['/usr/bin/fsl4.1-' command]);
    otherwise
        [~,y] = unix(command);
end