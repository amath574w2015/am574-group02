% =============================
% Reads output data
% =============================

function [output aux] = readData(outDir,whichRun,frameid)
    frameid = sprintf('%04d',frameid);
    
    % Read in solution data
    filename = [outDir 'run' num2str(whichRun) 'out.q' frameid];
    filename
    fid = fopen(filename);
    delimiterIn = ' ';
    headerlinesIn = 8;

    importedData = importdata(filename,delimiterIn,headerlinesIn);
    fclose(fid);

    data = importedData.data;
    header = importedData.textdata;
   
    x = data(:,1);
    out = data(:,2:end);
    
    output.x = x;
    output.data = out;
    
    foo = strread(header{1},'%s','delimiter',' ');
    maxPolyDegree = foo{1};
    
    foo = strread(header{2},'%s','delimiter',' ');
    nex = foo{1};
    
    foo = strread(header{3},'%s','delimiter',' ');
    xLeft = foo{1};
    foo = strread(header{4},'%s','delimiter',' ');
    xRight = foo{1};
    domain = [xLeft xRight];
    
    foo = strread(header{5},'%s','delimiter',' ');
    time = str2double(foo{1});
    output.t = time;
    
    % Read in auxilliary data
    filename = [outDir 'run' num2str(whichRun) 'out.aux' frameid];
    fid = fopen(filename);

    headerlinesIn = 3;
    importedData = importdata(filename,delimiterIn,headerlinesIn);
    fclose(fid);
    
    data = importedData.data;
    aux.elemCent = data(:,1);
    aux.localPolyDegree = data(:,2);
end
