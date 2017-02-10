% Vincent Yahna
% Matrix Computations
% Instructor: Ronald Taylor
% Storage for a Matrix
% with a diagonal, lower diagonal,
% and two upper diagonals

classdef sMatrix
    properties
        diagonal
        upper1
        upper2
        lower
    end
    methods
         %constructor
         %@param A - a matrix
         %@returns a square sMatrix with four diagonals
         %based on the original matrix
         function B = sMatrix(A)
             B.diagonal = diag(A);
             %make them the same length
             B.upper1 = B.diagonal;
             B.lower = B.diagonal;
             B.upper2 = B.diagonal;
             
             len = length(B.diagonal);
             B.lower(1) = 0;
             B.lower(2:len) = diag(A, -1);
             B.upper1(len) = 0;
             B.upper1(1:len-1) = diag(A, 1);
             B.upper2(len-1:len) = 0;
             B.upper2(1:len-2) = diag(A, 2);
             
         end%method
         
         %overloads the \ operator
         %This will be used to solve
         %the matrix equation Ax=b
         function y = mldivide(A,B)
         if not( isvector(B))
             %can only be used for vectors for now
             msgID = 'sMatrix:mldivide:NotAVector';
             msg = 'The second argument is not a vector';
             except = MException(msgID,msg);
             throw(except);
         end
            
             %Bi = bi - aici-1/Bi-1
             %Ci = ci - ai+1di-1/Bi-1
             beta = A.diagonal;
             gamma = A.upper1;
             y = B; %B must be a vector  
             
             for i =2:length(beta)
                 beta(i) = A.diagonal(i) - A.lower(i) * gamma(i-1)/ beta(i-1);
                 gamma(i) = gamma(i) - A.lower(i) * A.upper2(i-1)/ beta(i-1);
             end%for loop
             

             
             y(1) = y(1) / A.diagonal(1);
             for i=2:length(y)
                 y(i) = (y(i) - A.lower(i) * y(i-1))/beta(i);
             end
             
             
             %back substitution and get solution
             y(length(y) - 1) = y(length(y) - 1) - gamma(length(y) - 1)* y(length(y)) / beta(length(y) - 1);
             for i=length(y)-2:-1:1 
                 y(i) = y(i) - gamma(i)* y(i+1) / beta(i) - (A.upper2(i) * y(i+2))/beta(i);
             end
                 
         end%method
    end%methods
end%class