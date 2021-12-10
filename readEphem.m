function [t,X,Y,Z] = readEphem(fn)
    raw_ephem = fileread(fn);
    ephem = extractBetween(raw_ephem,'$$SOE','$$EOE');
    ephem = ephem{1}; % turn cell array to character vector

    %pull data from ephemeris readout (csv format)
    expr = '%f %*s %f %f %f %*f %*f %*f %*f %*f %*f';
    C = textscan(ephem,expr, 'Delimiter', ',');
    A = cell2mat(C).';
    t = A(1,:); X = A(2,:); Y = A(3,:); Z = A(4,:);
end

