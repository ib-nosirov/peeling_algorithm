% From 2019, modified 2021 to also save h5create .h5 files
rng(0);
% For HW 10
n   = 1e5;
%n   = 1e4;
m   = n;
% 1e5 x 1e5 is 74.5 GB
% 1e4 x 1e4 is .7 GB

if n < 1e5
    BIG = false;
else
    BIG = true;
end

MATFILE = false;  % both of these can be true, that's OK
H5FILE  = true;  % (so students have an option)

% for the 1e4 test case, 2019, 
%   saving by columns takes 4.7 sec per block (5 blocks)
%   and saving by rows, which I thought might be very slow,
%   takes 5.4 seconds, so not much difference

% for 1e5 case,
%File is 74.5 GB
%Using 38 blocks
%Filename: /rc_scratch/stbe1590/hugeSquareMatrix
%Using 0.1 GB of temporary memory
%Block 1 of 38 (4.4 to make, 58.4 to save)
%  That means, about 38 minutes expected!!
%  If I switch to save by row, took 59 sec to save,
%   so no major difference
% Total, by columns, is 2210.3 seconds = 36.8 min

MB        = 1024^2;
GB        = 1024*MB;
bytes_per_entry    = 8; % a double floating point entry is 8 bytes
%MB_limit  = 2*1024; % How many MB to load into RAM
MB_limit  = 1*1024; % How many MB to load into RAM
% Memory we'll use is BlockSize*m*bytes_per_entry, want < MB_limit

TotalSize = m*n*bytes_per_entry;
fprintf('File is %.1f GB\n', TotalSize/GB );

if TotalSize/GB < 1
    % To test the code, makea  smaller MB_limit
    MB_limit = TotalSize/MB/5;
end

LOCAL = ismac;
if LOCAL
    base = '/tmp';
else
    base = '/rc_scratch/stbe1590';
    %base = '/home/stbe1590';
end

if BIG
    name = 'hugeSquareMatrix';
    if LOCAL
        error('Do not run this on my laptop!');
    end
else
    name = 'largeSquareMatrix';
end

fileName = fullfile( base, name );
%if n < 1e5
%    if ismac % Stephen's laptop
%        fileName    = '/Users/srbecker/largeSquareMatrix';
%    else
%        fileName    = '/rc_scratch/stbe1590/largeSquareMatrix'; % 1e4 x 1e4
%    end
%else
%    fileName    = '/rc_scratch/stbe1590/hugeSquareMatrix'; % 1e5 x 1e5
%end

% How to save? pre-allocate? No, that doesn't really help
% see https://stackoverflow.com/a/27278859

BlockSize = floor( MB_limit*MB/(bytes_per_entry*m) );
NumberOfBlocks  = ceil( n/BlockSize );
fprintf('Using %d blocks\n', NumberOfBlocks );
fprintf('Filename: %s\n', fileName );

% The matrix is A = sigma*Noise + Left*Right
K   = 100;
Left    = randn(m,K)/sqrt(m);
Right   = randn(n,K)/sqrt(n);
sigma   = 1/sqrt(m*n);
fprintf('Using %.1f GB of temporary memory\n', (m+n)*K*bytes_per_entry/GB );

%% Save the data using "matfile"; code from 2019
if MATFILE
    matObj  = matfile( fileName, 'Writable', true );
    timeWrite = 0;
    timeCompute = 0;
    for j = 1:NumberOfBlocks
            rng(j); % so we can easily restart this

            startInd = (j-1)*BlockSize + 1;
            endInd = min( n, j*BlockSize );
            sizeBlock = endInd - startInd+1;

            tic
            A_block     = sigma*randn(m,sizeBlock);
            A_block     = A_block + Left*Right(startInd:endInd,1:K)';
            tC  = toc;
            timeCompute = timeCompute + tC;
            tic
            matObj.('A')(1:m,startInd:endInd ) = A_block;
            % matObj.('A')(startInd:endInd ,1:m) = A_block';% I assumed matfiles are also col major, but maybe not?
            tR  = toc;
            timeWrite = timeWrite + tR;
            fprintf('Block %d of %d (%.1f to make, %.1f to save)\n', j, NumberOfBlocks, tC, tR );
    end
    fprintf('Finished, took %.1f to make randn, %.1f to save to disk\n',timeCompute,timeWrite);
    clear matObj
    % forces it to write to file
    
    [~,result] = system(sprintf('du -h %s',[fileName,'.mat']) );
    disp(result)
    
    %% Read to memory for checking...
    if ~BIG
        mat_object = matfile( fileName );
        tic
        A          = mat_object.A;
        toc
    end
end

%% Save using h5create; added 2021

if H5FILE
    % Different options for Chunking.
    % Documentation at: https://portal.hdfgroup.org/display/HDF5/Chunking+in+HDF5
    % h5write(filename,ds,data,start,count,stride)  Ahh, "count" is how many elements per dimension...
    %   so we read/write from   start    to    start + count
    
    fileName = fullfile( base, [name,'.h5'] );
    %system(sprintf('rm %s',fileName) );
    % h5create( fileName , '/A',[m n], 'Datatype','double' ); % 'Deflate', 0 is default
    h5create( fileName , '/A',[m n], 'Datatype','double', 'Deflate', 0, 'Chunksize', [m, BlockSize] );
    
    timeIO = 0;
    timeCompute = 0;
    for j = 1:NumberOfBlocks
        rng(j); % so we can easily restart this
        
        startInd = (j-1)*BlockSize + 1;
        endInd = min( n, j*BlockSize );
        sizeBlock = endInd - startInd+1;
        
        tic
        A_block     = sigma*randn(m,sizeBlock);
        A_block     = A_block + Left*Right(startInd:endInd,1:K)';
        tC  = toc;
        timeCompute = timeCompute + tC;
        tic
        h5write( fileName, '/A', A_block, [1,startInd], [m,sizeBlock] );
        %h5write( fileName, '/U', U_block, [1,startInd], [m,sizeBlock] );
        tR  = toc;
        timeIO = timeIO + tR;
        fprintf('Block %d of %d (%.1f to make, %.1f to save)\n', j, NumberOfBlocks, tC, tR );
    end
    fprintf('Finished (h5write, no compression), took %.1f to make randn, %.1f to save to disk\n',timeCompute,timeIO);
    % MUCH faster!  0.4 to save to disk (vs 20 to 29 sec using matObj ).
    % Similar size on disk.
    [~,result] = system(sprintf('du -h %s',fileName) );
    disp(result)
    
    h5disp( fileName );
    
    
        %% Read it in, see if we get the right thing
    if ~BIG && MATFILE
        fprintf('\nStarting the checks\n')
        A_full   = zeros(m,n);
        timeIO = 0;
        for j = 1:NumberOfBlocks
            %rng(j); % so we can easily restart this
            
            startInd = (j-1)*BlockSize + 1;
            endInd = min( n, j*BlockSize );
            sizeBlock = endInd - startInd+1;
            
            tic
            A_full(:,startInd:endInd) = h5read( fileName, '/A', [1,startInd], [m,sizeBlock] );
            tR  = toc;
            timeIO = timeIO + tR;
            fprintf('Block %d of %d (%.1f to read)\n', j, NumberOfBlocks, tR );
        end
        fprintf('Finished (h5read, no compression), took %.1f to read from disk\n',timeIO);
        fprintf('Error in the data: %e\n', norm( A - A_full, 'fro' ) );
        
    end
end

%{
on blanca interactive (1 core??)

Small data,

File is 0.7 GB
Using 5 blocks
Filename: /rc_scratch/stbe1590/largeSquareMatrix
Using 0.0 GB of temporary memory
Block 1 of 5 (0.5 to make, 5.5 to save)
Block 2 of 5 (0.4 to make, 5.6 to save)
Block 3 of 5 (0.4 to make, 5.5 to save)
Block 4 of 5 (0.4 to make, 5.4 to save)
Block 5 of 5 (0.4 to make, 5.0 to save)
Finished, took 2.2 to make randn, 26.9 to save to disk
735M       /rc_scratch/stbe1590/largeSquareMatrix.mat

Elapsed time is 4.644349 seconds.
Block 1 of 5 (0.4 to make, 1.4 to save)
Block 2 of 5 (0.5 to make, 1.3 to save)
Block 3 of 5 (0.4 to make, 1.2 to save)
Block 4 of 5 (0.4 to make, 1.1 to save)
Block 5 of 5 (0.4 to make, 1.2 to save)
Finished (h5write, no compression), took 2.2 to make randn, 6.1 to save to disk
764M       /rc_scratch/stbe1590/largeSquareMatrix.h5

HDF5 largeSquareMatrix.h5
Group '/'
Dataset 'A'
        Size:  10000x10000
        MaxSize:  10000x10000
        Datatype:   H5T_IEEE_F64LE (double)
        ChunkSize:  10000x2000
        Filters:  deflate(0)
        FillValue:  0.000000

    Starting the checks
    Block 1 of 5 (0.5 to read)
    Block 2 of 5 (0.5 to read)
    Block 3 of 5 (0.5 to read)
    Block 4 of 5 (0.4 to read)
    Block 5 of 5 (0.4 to read)
    Finished (h5read, no compression), took 2.3 to read from disk
    Error in the data: 0.000000e+00


On HUGE data, n = 1e5
File is 74.5 GB
Using 75 blocks
Filename: /rc_scratch/stbe1590/hugeSquareMatrix
Using 0.1 GB of temporary memory
Block 1 of 75 (2.9 to make, 33.4 to save)
Block 2 of 75 (3.2 to make, 31.6 to save)
...
Block 74 of 75 (3.0 to make, 32.8 to save)
Block 75 of 75 (1.5 to make, 15.9 to save)
Finished, took 226.3 to make randn, 2419.7 to save to disk
73G /rc_scratch/stbe1590/hugeSquareMatrix.mat

[no faster if I use 16 CPUs]


and doing just .h5 file:

File is 74.5 GB
Using 75 blocks
Filename: /rc_scratch/stbe1590/hugeSquareMatrix
Using 0.1 GB of temporary memory
Block 1 of 75 (2.7 to make, 7.9 to save)
Block 2 of 75 (3.0 to make, 7.8 to save)
...
Block 73 of 75 (2.5 to make, 7.8 to save)
Block 74 of 75 (2.5 to make, 7.8 to save)
Block 75 of 75 (1.4 to make, 7.6 to save)
Finished (h5write, no compression), took 193.9 to make randn, 579.2 to save to disk
76G    /rc_scratch/stbe1590/hugeSquareMatrix.h5

HDF5 hugeSquareMatrix.h5
Group '/'
Dataset 'A'
Size:  100000x100000
MaxSize:  100000x100000
Datatype:   H5T_IEEE_F64LE (double)
ChunkSize:  100000x1342
Filters:  deflate(0)
FillValue:  0.000000

%}
