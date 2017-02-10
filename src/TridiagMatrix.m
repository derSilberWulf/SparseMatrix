%Vincent Yahna
%Matrix Computations
%Instructor: Ronald Taylor
%
%Sparse array for storing three diagonal bands
%in a square matrix
%Call A = TridiagMatrix(B)
%Input: B = an nx3 matrix where
%columns are the three diagonals
%and B(1,1) = 0 and B(n,n)=0
%Output: A = class of type TridiagMatrix
classdef TridiagMatrix


  properties
    % n x 3 matrix with diagonal bands as each column
    diags
  end
  
  methods
      
      function tdiagM = TridiagMatrix(A)
      %constructor

        tdiagM.diags = A;
      end
      
      
      function x = at(this, row, column)
      %indexing operator
      %Call: x = this.at(row, column)
      %Input: this = TridiagMatrix object
      %row = the row index
      %column = the column index
      %Output: x = value from the matrix
      
          if(row == column)
              %get from main diagonal
              x = this.diags(row , 2);
          elseif( row == column -1)
              %get from upper diagonal
              x = this.diags(row, 3);
          elseif( row -1 == column)
              %get from lower diagonal
              x = this.diags( row , 1);
          else
              x = 0;
          end
          
      end
      
      
      function [m,n] = size(X)
      %Call: [m,n] = size(X)
      %input: X - TridiagMatrix object
      %output: m = number of rows in the matrix
      %that is represented by this object
      %n = number of columns in the matrix
      %overloaded size function
      
          m = length(X.diags);
          n = m; %square matrix
      end
      
      
      function n = numel(X)
      %overloaded numel function
      %Call n = numel(X)
      %Input:X = TridiagMatrix
      %Output: n = number of rows/columns
      %in represented matrix
      
          n = length(X.diags);
          
      end
      

      

      function B = subsref(A,S)
      %overloaded subsref for indexing
      %indexing with parenthesis is supported
      %but colons cannot be used to denote
      %selecting an entire row or column
      %NOTE: THIS METHOD MAY BE CAUSING PROBLEMS
      %WITH GETTING CLASS PROPERTIES
      %Do not use this internally, use at instead
      %Call: b = This(m,n)
      %Input: This = TridiagMatrix
      %m = row number
      %n = column number
      %Output: b = number at the indexed position
     
          if(isequal(S.type, '()'))
              %use the at method, previously defined
              B = A.at(S.subs{1}, S.subs{2});
              
          else
              B=builtin('subsref', A,S);
          end
          
      end
      

      function C = mtimes(A,B)
      %overloaded * operator used for matrix multiplication
      %Left multiplications not yet implemented!
      %Call: A * B
      %Input: A = matrix or TridiagMatrix
      %B = matrix or TridiagMatrix
      %Either A or B is a TridiagMatrix (but not both)
      %Since left multilications are not implemented, it should be A
      
          %length check
          
          [Am, An] = size(A);
          [Bm, Bn] = size(B);
          if(An ~= Bm)
            msgID = 'TridiagMatrix:mtimes:MatrixDimensions';
            msg = 'The dimensions of the two matrices do not match';
            except = MException(msgID,msg);
            throw(except);
          end
          
          C = zeros(Am, Bn);
          
          %wasteful approach
          %for i=1:Am
          %    for j=1:Bn
          %        for k=1:An
           %         C(i, j) = C(i,j)+ A.at(i, k) * B(k,j);
           %       end
           %   end
         % end
         
          %determine whether A or B is TridiagMatrix
          if(isa(A, 'TridiagMatrix'))
              %right multiply
              for i=1:Bn
                C(1,i) = A.at(1,1) * B(1,i) + A.at(1,2) * B(2,i);
                C(An,i)= A.at(Am, An-1) *B(An-1,i) + A.at(Am, An) * B(An,i);
              end
              
              pivot = 1; %index where first number is in tridiagonal
              for i=2:An -1
                  for j=1:Bn
                      
                      C(i,j) = A.at(i, pivot) * B(pivot,j);
                      C(i,j) = C(i,j)+ A.at(i, pivot +1) * B(pivot +1,j);
                      C(i,j) = C(i,j)+ A.at(i, pivot +2) * B(pivot +2,j);
                  end
                  pivot = pivot +1;
              end
              
              
          else
              %left multiply
              msgID = 'TridiagMatrix:mtimes:LeftMultiplyNotImplemented';
              msg = 'Left multiplications have not yet been implemented for TridiagMatrix';
              except = MException(msgID,msg);
              throw(except);
          end
          
      
      end
  
  end
end