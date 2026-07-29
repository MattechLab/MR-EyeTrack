function [twix_obj rawData refScan] = fReadSiemensPerfusionRawData( filepathRawData )
%--------------------------------------------------------------------------
%
%   fReadRawDataSiemens     read raw data from Siemens scanner
%
%     [twix_obj rawData refScan] = fReadSiemensRawData( filepathRawData );
%
%     INPUT:    filepathRawData Raw file path
%
%     OUTPUT:   twix_obj        Twix object containing all the header info
%               rawData         Raw data
%               refScan         Reference scan data
%
%--------------------------------------------------------------------------

  % Add ReadRawDataSiemens directory to matlab path
    tmpdir = which('fReadSiemensRawData'); % find fReadSiemensRawData.m
	[tmpdir ~] = fileparts(tmpdir);
    addpath([tmpdir filesep 'ReadRawDataSiemens']);

    fprintf('Start reading Siemens raw data on %s\n', datestr(now));
    tic;
    
  % Reading RAW DATA header (raw data are actually read only when needed);
    fprintf('... read raw data header ...\n');
    %twix_obj = mapVBVD(filepathRawData);
    twix_obj = mapVBVD_jy(filepathRawData);
    twix_obj.image.flagIgnoreSeg = true; %(essential to read the data correctly)
    twix_obj.refscan.flagIgnoreSeg = true; %(essential to read the data correctly)

  % READ the complete RAW DATA
  %     Raw data format = 2*Np x Nc x Ns
  %     where Np is the number of readout point, Nc is the number of
  %     channels, Ns is the total number of shots, and the factor 2 
  %     accounts for the oversampling in the readout direction.
  %
  % Order of raw data:
  %  1) Columns
  %  2) Channels/Coils
  %  3) Lines
  %  4) Partitions
  %  5) Slices
  %  6) Averages
  %  7) (Cardiac-) Phases
  %  8) Contrasts/Echoes
  %  9) Measurements
  % 10) Sets
  % 11) Segments
  % 12) Ida
  % 13) Idb
  % 14) Idc
  % 15) Idd
  % 16) Ide
    fprintf('... read raw data ...\n');
    rawData = twix_obj.image{''}; 
    refScan = twix_obj.refscan{''};
    
  % Permuting raw data to satisfy the following convention
    [imageInfo rawData]   = permuteDataToSatisfyJYConvention(twix_obj.image, rawData);
    [refscanInfo refScan] = permuteDataToSatisfyJYConvention(twix_obj.refscan, refScan);
  
    twix_obj.image = imageInfo;
    twix_obj.refscan = refscanInfo;
  
    
function [image, rawData] = permuteDataToSatisfyJYConvention(image, rawData)
    
    
  % Permuting raw data to satisfy the following convention
  %     Raw data format = 2*Np x Ns x Nc
  %     where Np is the number of readout point, Nc is the number of
  %     channels, Ns is the total number of shots, and the factor 2 
  %     accounts for the oversampling in the readout direction.
  %
  % Permute data to follow convention [Nx Ny Nz Nt Nc], i.e.,
  %  1) Columns
  %  2) Lines 
  %  3) Slices
  %  4) Measurements
  %  5) Channels/Coils
  %  6) Averages
  %  7) (Cardiac-) Phases 
  %  8) Contrasts/Echoes
  %  9) Partitions
  % 10) Sets
  % 11) Segments
  % 12) Ida
  % 13) Idb
  % 14) Idc
  % 15) Idd
  % 16) Ide
  
  % Size of the non-squeezed raw data
    dataSize = image.dataSize; 
    
  % "Un-squeeze" raw data
    rawData  = reshape( rawData, dataSize ); 
    
  % Permute dimension to follow Jerome's convention
    fprintf('... permute raw data ...\n');
    order_jy = [1 3 5 9 2 6 7 8 4 10 11 12 13 14 15 16];
    rawData = permute( rawData, order_jy );

    tmpCell = image.dataDims;
    for i = 1:length(order_jy)
        tmpCell{i} = image.dataDims{order_jy(i)};
    end
    image.dataDims = tmpCell;
    image.sqzDims  = tmpCell;
    
    fprintf('DONE!!!\n')
    
  % Compute elapsed time
    toc;

    
    
