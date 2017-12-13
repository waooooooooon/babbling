function [succ] = HIImageConvert2GIF( regu_path, path_out, delay_ms )

setenv('PATH', [getenv('PATH') ':/usr/local/bin']) ;
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/') ;

if nargin < 3
    delay_ms = 1000 ;
end

[status, cmdout] = system('convert') ;
if status ~= 0
    warning( cmdout ) ;
    succ = false ;
    return 
end

delay = round( delay_ms / 10 ) ;
[status, cmdout] = system(['convert -colorspace sRGB -loop 0 -delay ' num2str(delay) ' ' regu_path ' ' path_out]) ;
if status ~= 0
    warning( cmdout ) ;
    succ = false ;
    return 
end

succ = true ;
