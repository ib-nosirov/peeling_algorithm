% Author: Stephen Becker
function SaveSquareMatrix(A,blockSize)
% Different opts for Chunking.
% Documentation at: https://portal.hdfgroup.org/display/HDF5/Chunking+in+HDF5
% h5write(filename,ds,data,start,count,stride)  Ahh, "count" is how many elements per dimension...
%   so we read/write from   start    to    start + countion
    base = '/home/ibrohim/Research/S22/peeling_algorithm/';
    name = 'BigAMatrix';
    n = length(A);
    numberOfBlocks = ceil(n/blockSize);
    fileName = fullfile(base,[name,num2str(n),'.h5']);
    %system(sprintf('rm %s',fileName) );
    % h5create( fileName , '/A',[m n], 'Datatype','double' ); % 'Deflate', 0 is default
    % 'ChunkSize',[n,blockSize] -> we will need to chunk column-wise.
    h5create(fileName ,'/A',[n,n],'Datatype','double','Deflate',0,...
        'Chunksize',[n,blockSize]);
    
    %timeIO = 0;
    %timeCompute = 0;
    for j = 1:numberOfBlocks        
        startInd = (j-1)*blockSize + 1;
        endInd = min(n,j*blockSize);
        sizeBlock = endInd - startInd+1;
        
%        tic
        A_block = A(:,startInd:endInd);
%        tC  = toc;
%        timeCompute = timeCompute + tC;
%        tic
        h5write(fileName,'/A',A_block,[1,startInd],[n,sizeBlock]);
        %h5write( fileName, '/U', U_block, [1,startInd], [m,sizeBlock] );
%        tR  = toc;
%        timeIO = timeIO + tR;
%        fprintf('Block %d of %d (%.1f to make, %.1f to save)\n',j,...
%            numberOfBlocks,tC,tR);
    end
%    [~,result] = system(sprintf('du -h %s',fileName) );
%    disp(result)
    
%    h5disp(fileName);
end