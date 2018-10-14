function [lineCount] = line_counter(path)
% This function finds the number of lines in a file. Note that it is
% designed to only work on linux and Mac machines. Install the GnuWin32 on
% your PC and modify this function to use it on windows.

if (~ispc)
    [status, cmdout]= system(['wc -l ', path]);
    if(status~=1)
        scanCell = textscan(cmdout,'%u %s');
        lineCount = scanCell{1};
    else
        fprintf(1,'Failed to find line count of %s\n',filenameOfInterest.txt);
        lineCount = -1;
    end
else
    fprintf(1,'Sorry, I don''t know what the equivalent is for a windows system\n');
    lineCount = -1;
end
end