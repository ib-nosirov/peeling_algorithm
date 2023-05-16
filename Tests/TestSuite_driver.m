%% write a set of 10 tests for different sizes and kinds of matrices.
numTests = 10
fprintf('MakeHODLRatrix_test: \n')
for ii = 1:numTests
    fprintf('Test No. %d\n',ii)
    MakeHODLRMtrx_test()
end

%% write a set of 10 tests for different sizes and kinds of matrices.
for ii = 1:numTests
    fprintf('Test No. %d\n',ii)
    HODLRMatVec_test()
end

%% write a set of 10 tests for different sizes and kinds of matrices.
for ii = 1:numTests
    fprintf('Test No. %d\n',ii)
    Lanczos_test()
end

%% write a set of 10 tests for different sizes and kinds of matrices.
for ii = 1:numTests
    fprintf('Test No. %d\n',ii)
    SLQ_test()
end
%% write a set of 10 tests for different sizes and kinds of matrices.
for ii = 1:numTests
    fprintf('Test No. %d\n',ii)
    TwoDim_test()
end