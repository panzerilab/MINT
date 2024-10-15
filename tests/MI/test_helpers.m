classdef test_helpers < matlab.unittest.TestCase
    methods (Test)
        
        function testEmptyInputFindTrial(testCase)
            testCase.verifyError(@()findtrial(), 'findtrial:notEnoughInput');
        end

        function testValidInputFindTrial(testCase)        
            %Equal Trial Stim
            Nt = [4; 4; 4];
            maxNt = 4;
            Ns = 3;
            expectedResult = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12];
            assert(isequal(findtrial(Nt, maxNt, Ns), expectedResult), 'Expected result not given');
 
            %Diff Trial Stim
            Nt = [3; 5; 2];
            maxNt = 5;
            Ns = 3;
            expectedResult = [1; 2; 3; 6; 7; 8; 9; 10; 11; 12];
            assert(isequal(findtrial(Nt, maxNt, Ns), expectedResult), 'Expected result not given');

            %Zeros
            Nt = zeros(3, 1);  
            maxNt = 0;         
            Ns = 3;            
            assert(isempty(findtrial(Nt, maxNt, Ns)), 'Expected result not given');
        end

        function testNaNInputFindTrial(testCase) 
            Nt = NaN; 
            maxNt = 4;
            Ns = 3; 
            testCase.verifyError(@()findtrial(Nt, maxNt, Ns), 'findtrial:NaNInput');
        end 

        function testInvalidInputFindTrial(testCase)
            Nt = [3; 5; 2];
            maxNt = 4;
            Ns = 3;
            testCase.verifyError(@()findtrial(Nt, maxNt, Ns), 'findtrial:InvalidInput');
        end
    end
end