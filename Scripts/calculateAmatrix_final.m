function [A,coefvalues] = calculateAmatrix_final(kmeanMatrix_joint)
% modelled metabolic flux coefficients - vector f
coefvalues = {'f1'	...
     'HSIctr', 'HLIctr', 'HColctr',...
     'HSIhfd', 'HLIhfd', 'HColhfd',...
     'B1LIctr', 'B1LIhfd'};
% coefficient matrix to solve linear optimization task in the form Af=0
 A = [kmeanMatrix_joint(1,1)-kmeanMatrix_joint(1,2)	1	0 0 0	0	0	0 0 
kmeanMatrix_joint(1,2)-kmeanMatrix_joint(1,3)	1	0 0 0	0	0	0	0 
kmeanMatrix_joint(1,3)-kmeanMatrix_joint(1,4)	0	1	0 0 0	0	1	0 
kmeanMatrix_joint(1,4)-kmeanMatrix_joint(1,5)	0	0 1 0	0	0	0 0
kmeanMatrix_joint(1,5)-kmeanMatrix_joint(1,6)	0	0 1 0   0	0	0	0 
kmeanMatrix_joint(2,1)-kmeanMatrix_joint(2,2)	1	0	0	0 0 0	0	0 
kmeanMatrix_joint(2,2)-kmeanMatrix_joint(2,3)	1	0	0	0	0 0 0	0 
kmeanMatrix_joint(2,3)-kmeanMatrix_joint(2,4)	0	1	0 0 0	0	0	0 
kmeanMatrix_joint(2,4)-kmeanMatrix_joint(2,5)	0	0 1	0 0	0	0	0
kmeanMatrix_joint(2,5)-kmeanMatrix_joint(2,6)	0	0	1 0 0	0	0	0 
kmeanMatrix_joint(3,1)-kmeanMatrix_joint(3,2)	0	0	0 1	0 0	0	0
kmeanMatrix_joint(3,2)-kmeanMatrix_joint(3,3)	0	0	0 1	0 0	0	0
kmeanMatrix_joint(3,3)-kmeanMatrix_joint(3,4)	0	0	0	0 1 0 	0 1
kmeanMatrix_joint(3,4)-kmeanMatrix_joint(3,5)	0	0	0	0 0 1	0	0
kmeanMatrix_joint(3,5)-kmeanMatrix_joint(3,6)	0	0	0	0 0 1	0	0 
kmeanMatrix_joint(4,1)-kmeanMatrix_joint(4,2)	0	0	0 1	0 0	0	0
kmeanMatrix_joint(4,2)-kmeanMatrix_joint(4,3)	0	0	0 1	0 0	0	0
kmeanMatrix_joint(4,3)-kmeanMatrix_joint(4,4)	0	0	0	0 1	 0 0	0 
kmeanMatrix_joint(4,4)-kmeanMatrix_joint(4,5)	0	0	0	0 0 1	0	0 
kmeanMatrix_joint(4,5)-kmeanMatrix_joint(4,6)	0	0	0	0 0	1   0   0 
];

end