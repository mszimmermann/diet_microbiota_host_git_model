function [RA,b] = calculateRAmatrix_final(x)
%calculate reverse A matrix with flux formulas for each of the data points
% to solve reverse task RA * M = b where M - metabolite abundance values
% and input x - flux parameter values calculated in the forward problem
RA = [x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
      0    x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
      0 0 x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
      0 0 0 x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
      0 0 0 0 x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
      
      0 0 0 0 0 0 x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
      0 0 0 0 0 0 0    x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
      0 0 0 0 0 0 0 0 x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
      0 0 0 0 0 0 0 0 0 x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 
      0 0 0 0 0 0 0 0 0 0 x(1) -x(1) 0 0 0 0 0 0 0 0 0 0 0 0 
      
      0 0 0 0 0 0 0 0 0 0 0 0 x(1) -x(1) 0 0 0 0 0 0 0 0 0 0  
      0 0 0 0 0 0 0 0 0 0 0 0 0    x(1) -x(1) 0 0 0 0 0 0 0 0 0  
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(1) -x(1) 0 0 0 0 0 0 0 0 
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(1) -x(1) 0 0 0 0 0 0 0 
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(1) -x(1) 0 0 0 0 0 0  
      
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(1) -x(1) 0 0 0 0   
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0    x(1) -x(1) 0 0 0   
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(1) -x(1) 0 0 
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(1) -x(1) 0 
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(1) -x(1) 
];

b = [-x(2)
    -x(2)
    -x(3)-x(8)
    -x(4)
    -x(4)
    -x(2)
    -x(2)
    -x(3)
    -x(4)
    -x(4)
    -x(5)
    -x(5)
    -x(6)-x(9)
    -x(7)
    -x(7)
    -x(5)
    -x(5)
    -x(6)
    -x(7)
    -x(7)];

end