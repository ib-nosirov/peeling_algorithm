function A_part = ReadBigMatrix(fileName,n,numberOfBlocks,blockSize)
    A_part = [n,blockSize];
    %for j = 1:numberOfBlocks
%        rng(j); % so we can easily restart this
        startInd = (j-1)*blockSize + 1;
        endInd = min(n,j*blockSize);
        sizeBlock = endInd - startInd+1;

    %    tic
        A_part(:,startInd:endInd) = h5read(fileName,'/A',...
            [1,startInd],[m,sizeBlock]);
    %    tR  = toc;
    %    timeIO = timeIO + tR;
    %    fprintf('Block %d of %d (%.1f to read)\n',j,numberOfBlocks,tR);
    %end
    fprintf(['Finished (h5read, no compression), took %.1f to read ',...
        'from disk\n'],timeIO);
end