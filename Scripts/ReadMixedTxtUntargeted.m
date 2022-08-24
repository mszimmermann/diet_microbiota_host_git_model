function TempStruct = ReadMixedTxtUntargeted(datafile, delimiter)
  TempStruct = struct;
  fid = fopen(datafile,'r');   %# Open the file
  lineArray = cell(50000,1);     %# Preallocate a cell array (ideally slightly   %#   larger than is needed)
  lineIndex = 1;               %# Index of cell to place the next line in
  nextLine = fgetl(fid);       %# Read the first line from the file
  while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
    lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
    lineIndex = lineIndex+1;          %# Increment the line index
    nextLine = fgetl(fid);            %# Read the next line from the file
  end
  fclose(fid);                 %# Close the file
  lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
  for iLine = 1:lineIndex-1              %# Loop over lines
    lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                        'Delimiter',delimiter);
    lineData = lineData{1};              %# Remove cell encapsulation
    if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
      lineData{end+1} = '';                     %#   ends with a delimiter
    end
    lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
  end
  %get the lines and columns containin 'Compound'
  % data should be between the columns
  lineCompound = find(cellfun(@(x) contains(x, 'Compound'),lineArray(:,1)));
  columnCompound = find(cellfun(@(x) contains(x, 'Compound'),lineArray(lineCompound,:)));
  columnRT = find(cellfun(@(x) contains(x, 'Retention'),lineArray(lineCompound,:)));
  columnSpectrum = find(cellfun(@(x) contains(x, 'Spectrum'),lineArray(lineCompound,:)));
  columnMass = find(cellfun(@(x) contains(x, 'Mass'),lineArray(lineCompound,:)));
  
  % if there are two "Compound" columns, take the data between them.
  % else, take the data between first column and next column
  if length(columnCompound)>1
      nondatacolumn = min(min(columnRT,columnCompound(2)),min(columnSpectrum,columnMass));
  else
      nondatacolumn = min(columnRT,min(columnSpectrum,columnMass));
  end
  TempStruct.IntensitiesRaw = cellfun(@(x) str2double(x), ...
                    lineArray(lineCompound+1:end, ...
                              columnCompound+1:nondatacolumn-1));
  TempStruct.SampleNames =  lineArray(lineCompound,...
                                  columnCompound+1:nondatacolumn-1);      
  
  TempStruct.IntensitiesRaw(TempStruct.IntensitiesRaw==0) = 1;
  % import RT data
  TempStruct.RT =  cellfun(@(x) str2double(x),...
                         lineArray(lineCompound+1:end, columnRT));
  TempStruct.Mass =  cellfun(@(x) str2double(x),...
                         lineArray(lineCompound+1:end, columnMass));

  % import compounds
  TempStruct.Compounds =  lineArray(lineCompound+1:end, 1);
  TempStruct.CompositeSpectrum =  lineArray(lineCompound+1:end,...
                                            columnSpectrum);

end