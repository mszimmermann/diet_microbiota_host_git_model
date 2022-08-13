function [A,coefvalues] = calculateAmatrix_final(kmeanMatrix_joint)
% coefvalues = {'f1'	...
%      'HSI', 'HLI',...
%      'B1LIctr'	'B1LIhfd'};
%      A = [kmeanMatrix_joint(1,1)-kmeanMatrix_joint(1,2)	1	0	0	0
%         kmeanMatrix_joint(1,2)-kmeanMatrix_joint(1,3)	1	0	0	0
%         kmeanMatrix_joint(1,3)-kmeanMatrix_joint(1,4)	0	1	1	0
%         kmeanMatrix_joint(1,4)-kmeanMatrix_joint(1,5)	0	1	1	0
%         kmeanMatrix_joint(1,5)-kmeanMatrix_joint(1,6)	0	0	0	0
%         kmeanMatrix_joint(2,1)-kmeanMatrix_joint(2,2)	1	0	0	0
%         kmeanMatrix_joint(2,2)-kmeanMatrix_joint(2,3)	1	0	0	0
%         kmeanMatrix_joint(2,3)-kmeanMatrix_joint(2,4)	0	1	0	0
%         kmeanMatrix_joint(2,4)-kmeanMatrix_joint(2,5)	0	1	0	0
%         kmeanMatrix_joint(2,5)-kmeanMatrix_joint(2,6)	0	0	0	0
%         kmeanMatrix_joint(3,1)-kmeanMatrix_joint(3,2)	1	0	0	0
%         kmeanMatrix_joint(3,2)-kmeanMatrix_joint(3,3)	1	0	0	0
%         kmeanMatrix_joint(3,3)-kmeanMatrix_joint(3,4)	0	1	0	1
%         kmeanMatrix_joint(3,4)-kmeanMatrix_joint(3,5)	0	1	0	1
%         kmeanMatrix_joint(3,5)-kmeanMatrix_joint(3,6)	0	0	0	0
%         kmeanMatrix_joint(4,1)-kmeanMatrix_joint(4,2)	1	0	0	0
%         kmeanMatrix_joint(4,2)-kmeanMatrix_joint(4,3)	1	0	0	0
%         kmeanMatrix_joint(4,3)-kmeanMatrix_joint(4,4)	0	1	0	0
%         kmeanMatrix_joint(4,4)-kmeanMatrix_joint(4,5)	0	1	0	0
%         kmeanMatrix_joint(4,5)-kmeanMatrix_joint(4,6)	0	0	0	0
%         ];
% end

 
% coefvalues = {'f1'	...
%     'H1SIctr'	'H2LIctr'	'H1hfd'	'H2hfd'	...
%     'BLI1ctr'	'BLI2ctr'	'BLI3ctr'	...
%     'BLI1hfd'	'BLI2hfd'	'BLI3hfd'};
% Aeq = [0 0 0 0 0 1 -1 0 0 0 0;...
%        0 0 0 0 0 0 -1 1 0 0 0;...
%        0 0 0 0 0 0 0 0 1 -1 0;...
%        0 0 0 0 0 0 0 0 0 -1 1];
% beq = [0; 0; 0; 0];   
%      
% A = [kmeanMatrix_joint(1,1)-kmeanMatrix_joint(1,2)	1	0	0	0	0	0	0	0	0	0
% kmeanMatrix_joint(1,2)-kmeanMatrix_joint(1,3)	1	0	0	0	0	0	0	0	0	0
% kmeanMatrix_joint(1,3)-kmeanMatrix_joint(1,4)	0	1	0	0	1	0	0	0	0	0
% kmeanMatrix_joint(1,4)-kmeanMatrix_joint(1,5)	0	1	0	0	0	1	0	0	0	0
% kmeanMatrix_joint(1,5)-kmeanMatrix_joint(1,6)	0	1	0	0	0	0	1	0	0	0
% kmeanMatrix_joint(2,1)-kmeanMatrix_joint(2,2)	1	0	0	0	0	0	0	0	0	0
% kmeanMatrix_joint(2,2)-kmeanMatrix_joint(2,3)	1	0	0	0	0	0	0	0	0	0
% kmeanMatrix_joint(2,3)-kmeanMatrix_joint(2,4)	0	1	0	0	0	0	0	0	0	0
% kmeanMatrix_joint(2,4)-kmeanMatrix_joint(2,5)	0	1	0	0	0	0	0	0	0	0
% kmeanMatrix_joint(2,5)-kmeanMatrix_joint(2,6)	0	1	0	0	0	0	0	0	0	0
% kmeanMatrix_joint(3,1)-kmeanMatrix_joint(3,2)	0	0	1	0	0	0	0	0	0	0
% kmeanMatrix_joint(3,2)-kmeanMatrix_joint(3,3)	0	0	1	0	0	0	0	0	0	0
% kmeanMatrix_joint(3,3)-kmeanMatrix_joint(3,4)	0	0	0	1	0	0	0	0	0	0
% kmeanMatrix_joint(3,4)-kmeanMatrix_joint(3,5)	0	0	0	1	0	0	0	0	0	0
% kmeanMatrix_joint(3,5)-kmeanMatrix_joint(3,6)	0	0	0	1	0	0	0	0	0	0
% kmeanMatrix_joint(4,1)-kmeanMatrix_joint(4,2)	0	0	1	0	0	0	0	0	0	0
% kmeanMatrix_joint(4,2)-kmeanMatrix_joint(4,3)	0	0	1	0	0	0	0	0	0	0
% kmeanMatrix_joint(4,3)-kmeanMatrix_joint(4,4)	0	0	0	1	0	0	0	1	0	0
% kmeanMatrix_joint(4,4)-kmeanMatrix_joint(4,5)	0	0	0	1	0	0	0	0	1	0
% kmeanMatrix_joint(4,5)-kmeanMatrix_joint(4,6)	0	0	0	1	0	0	0	0	0	1];
% end

% coefvalues = {'f1'	...
%      'HSIabs', 'HLIabs',...
%      'H1SIctr'	'H2LIctr'...
%      'H1SIhfd'	'H2LIhfd'...
%      'B1LIctr'	'B1LIhfd'};
% 
% A = [kmeanMatrix_joint(1,1)-kmeanMatrix_joint(1,2)	-kmeanMatrix_joint(1,2)	0	1	0	0	0	0	0
% kmeanMatrix_joint(1,2)-kmeanMatrix_joint(1,3)	-kmeanMatrix_joint(1,3)	0	1	0	0	0	0	0
% kmeanMatrix_joint(1,3)-kmeanMatrix_joint(1,4)	0	-kmeanMatrix_joint(1,4)	0	1	0	0	1	0
% kmeanMatrix_joint(1,4)-kmeanMatrix_joint(1,5)	0	-kmeanMatrix_joint(1,5)	0	1	0	0	1	0
% kmeanMatrix_joint(1,5)-kmeanMatrix_joint(1,6)	0	-kmeanMatrix_joint(1,6)	0	1	0	0	1	0
% kmeanMatrix_joint(2,1)-kmeanMatrix_joint(2,2)	-kmeanMatrix_joint(2,2)	0	1	0	0	0	0	0
% kmeanMatrix_joint(2,2)-kmeanMatrix_joint(2,3)	-kmeanMatrix_joint(2,3)	0	1	0	0	0	0	0
% kmeanMatrix_joint(2,3)-kmeanMatrix_joint(2,4)	0	-kmeanMatrix_joint(2,4)	0	1	0	0	0	0
% kmeanMatrix_joint(2,4)-kmeanMatrix_joint(2,5)	0	-kmeanMatrix_joint(2,5)	0	1	0	0	0	0
% kmeanMatrix_joint(2,5)-kmeanMatrix_joint(2,6)	0	-kmeanMatrix_joint(2,6)	0	1	0	0	0	0
% kmeanMatrix_joint(3,1)-kmeanMatrix_joint(3,2)	-kmeanMatrix_joint(3,2)	0	0	0	1	0	0	0
% kmeanMatrix_joint(3,2)-kmeanMatrix_joint(3,3)	-kmeanMatrix_joint(3,3)	0	0	0	1	0	0	0
% kmeanMatrix_joint(3,3)-kmeanMatrix_joint(3,4)	0	-kmeanMatrix_joint(3,4)	0	0	0	1	0	1
% kmeanMatrix_joint(3,4)-kmeanMatrix_joint(3,5)	0	-kmeanMatrix_joint(3,5)	0	0	0	1	0	1
% kmeanMatrix_joint(3,5)-kmeanMatrix_joint(3,6)	0	-kmeanMatrix_joint(3,6)	0	0	0	1	0	1
% kmeanMatrix_joint(4,1)-kmeanMatrix_joint(4,2)	-kmeanMatrix_joint(4,2)	0	0	0	1	0	0	0
% kmeanMatrix_joint(4,2)-kmeanMatrix_joint(4,3)	-kmeanMatrix_joint(4,3)	0	0	0	1	0	0	0
% kmeanMatrix_joint(4,3)-kmeanMatrix_joint(4,4)	0	-kmeanMatrix_joint(4,4)	0	0	0	1	0	0
% kmeanMatrix_joint(4,4)-kmeanMatrix_joint(4,5)	0	-kmeanMatrix_joint(4,5)	0	0	0	1	0	0
% kmeanMatrix_joint(4,5)-kmeanMatrix_joint(4,6)	0	-kmeanMatrix_joint(4,6)	0	0	0	1	0	0
% ];
% end
% A = [kmeanMatrix_joint(1,1)-kmeanMatrix_joint(1,2)	1	0	0	0	0	0
%     kmeanMatrix_joint(1,2)-kmeanMatrix_joint(1,3)	1	0	0	0	0	0
%     kmeanMatrix_joint(1,3)-kmeanMatrix_joint(1,4)	0	1	0	0	1	0
%     kmeanMatrix_joint(1,4)-kmeanMatrix_joint(1,5)	0	0	0	0	0	0
%     kmeanMatrix_joint(1,5)-kmeanMatrix_joint(1,6)	0	0	0	0	0	0
%     kmeanMatrix_joint(2,1)-kmeanMatrix_joint(2,2)	1	0	0	0	0	0
%     kmeanMatrix_joint(2,2)-kmeanMatrix_joint(2,3)	1	0	0	0	0	0
%     kmeanMatrix_joint(2,3)-kmeanMatrix_joint(2,4)	0	1	0	0	0	0
%     kmeanMatrix_joint(2,4)-kmeanMatrix_joint(2,5)	0	0	0	0	0	0
%     kmeanMatrix_joint(2,5)-kmeanMatrix_joint(2,6)	0	0	0	0	0	0
%     kmeanMatrix_joint(3,1)-kmeanMatrix_joint(3,2)	0	0	1	0	0	0
%     kmeanMatrix_joint(3,2)-kmeanMatrix_joint(3,3)	0	0	1	0	0	0
%     kmeanMatrix_joint(3,3)-kmeanMatrix_joint(3,4)	0	0	0	1	0	1
%     kmeanMatrix_joint(3,4)-kmeanMatrix_joint(3,5)	0	0	0	0	0	0
%     kmeanMatrix_joint(3,5)-kmeanMatrix_joint(3,6)	0	0	0	0	0	0
%     kmeanMatrix_joint(4,1)-kmeanMatrix_joint(4,2)	0	0	1	0	0	0
%     kmeanMatrix_joint(4,2)-kmeanMatrix_joint(4,3)	0	0	1	0	0	0
%     kmeanMatrix_joint(4,3)-kmeanMatrix_joint(4,4)	0	0	0	1	0	0
%     kmeanMatrix_joint(4,4)-kmeanMatrix_joint(4,5)	0	0	0	0	0	0
%     kmeanMatrix_joint(4,5)-kmeanMatrix_joint(4,6)	0	0	0	0	0	0
%     ];
% end
% coefvalues = {'f1'	...
%      'HSIctr', 'HLIctr',...
%      'HSIhfd', 'HLIhfd',...
%      'B1LIctr'	'B1LIhfd'};
%  A = [kmeanMatrix_joint(1,1)-kmeanMatrix_joint(1,2)	1	0	0	0	0	0
% kmeanMatrix_joint(1,2)-kmeanMatrix_joint(1,3)	1	0	0	0	0	0
% kmeanMatrix_joint(1,3)-kmeanMatrix_joint(1,4)	0	1	0	0	1	0
% kmeanMatrix_joint(1,4)-kmeanMatrix_joint(1,5)	0	1	0	0	1	0
% kmeanMatrix_joint(1,5)-kmeanMatrix_joint(1,6)	0	0	0	0	0	0
% kmeanMatrix_joint(2,1)-kmeanMatrix_joint(2,2)	1	0	0	0	0	0
% kmeanMatrix_joint(2,2)-kmeanMatrix_joint(2,3)	1	0	0	0	0	0
% kmeanMatrix_joint(2,3)-kmeanMatrix_joint(2,4)	0	1	0	0	0	0
% kmeanMatrix_joint(2,4)-kmeanMatrix_joint(2,5)	0	1	0	0	0	0
% kmeanMatrix_joint(2,5)-kmeanMatrix_joint(2,6)	0	0	0	0	0	0
% kmeanMatrix_joint(3,1)-kmeanMatrix_joint(3,2)	0	0	1	0	0	0
% kmeanMatrix_joint(3,2)-kmeanMatrix_joint(3,3)	0	0	1	0	0	0
% kmeanMatrix_joint(3,3)-kmeanMatrix_joint(3,4)	0	0	0	1	0	1
% kmeanMatrix_joint(3,4)-kmeanMatrix_joint(3,5)	0	0	0	1	0	1
% kmeanMatrix_joint(3,5)-kmeanMatrix_joint(3,6)	0	0	0	0	0	0
% kmeanMatrix_joint(4,1)-kmeanMatrix_joint(4,2)	0	0	1	0	0	0
% kmeanMatrix_joint(4,2)-kmeanMatrix_joint(4,3)	0	0	1	0	0	0
% kmeanMatrix_joint(4,3)-kmeanMatrix_joint(4,4)	0	0	0	1	0	0
% kmeanMatrix_joint(4,4)-kmeanMatrix_joint(4,5)	0	0	0	1	0	0
% kmeanMatrix_joint(4,5)-kmeanMatrix_joint(4,6)	0	0	0	0	0	0
% ];
% coefvalues = {'f1'	...
%      'HSIctr', 'HLIctr', 'HColctr',...
%      'HSIhfd', 'HLIhfd', 'HColhfd',...
%      'B1LIctr', 'B1Colctr',	'B1LIhfd', 'BColhfd'};
%  A = [kmeanMatrix_joint(1,1)-kmeanMatrix_joint(1,2)	1	0 0 0	0	0	0 0 0	0
% kmeanMatrix_joint(1,2)-kmeanMatrix_joint(1,3)	1	0 0 0	0	0	0	0 0 0
% kmeanMatrix_joint(1,3)-kmeanMatrix_joint(1,4)	0	1	0 0 0	0	1	0 0 0
% kmeanMatrix_joint(1,4)-kmeanMatrix_joint(1,5)	0	0 1 0	0	0	0	1 0 0
% kmeanMatrix_joint(1,5)-kmeanMatrix_joint(1,6)	0	0	0 0 0	0	0	0 0 0
% kmeanMatrix_joint(2,1)-kmeanMatrix_joint(2,2)	1	0	0	0 0 0	0	0 0 0
% kmeanMatrix_joint(2,2)-kmeanMatrix_joint(2,3)	1	0	0	0	0 0 0	0 0 0
% kmeanMatrix_joint(2,3)-kmeanMatrix_joint(2,4)	0	1	0 0 0	0	0	0 0 0
% kmeanMatrix_joint(2,4)-kmeanMatrix_joint(2,5)	0	0 1	0 0	0	0	0 0 0
% kmeanMatrix_joint(2,5)-kmeanMatrix_joint(2,6)	0	0	0 0 0	0	0	0 0 0
% kmeanMatrix_joint(3,1)-kmeanMatrix_joint(3,2)	0	0	0 1	0 0	0	0 0 0
% kmeanMatrix_joint(3,2)-kmeanMatrix_joint(3,3)	0	0	0 1	0 0	0	0 0 0
% kmeanMatrix_joint(3,3)-kmeanMatrix_joint(3,4)	0	0	0	0 1 0 	0	0 1 0
% kmeanMatrix_joint(3,4)-kmeanMatrix_joint(3,5)	0	0	0	0 0 1	0	0 0 1
% kmeanMatrix_joint(3,5)-kmeanMatrix_joint(3,6)	0	0	0	0 0 0	0	0 0 0
% kmeanMatrix_joint(4,1)-kmeanMatrix_joint(4,2)	0	0	0 1	0 0	0	0 0 0
% kmeanMatrix_joint(4,2)-kmeanMatrix_joint(4,3)	0	0	0 1	0 0	0	0 0 0
% kmeanMatrix_joint(4,3)-kmeanMatrix_joint(4,4)	0	0	0	0 1	 0 0	0 0 0
% kmeanMatrix_joint(4,4)-kmeanMatrix_joint(4,5)	0	0	0	0 0 1	0	0 0 0
% kmeanMatrix_joint(4,5)-kmeanMatrix_joint(4,6)	0	0	0	0	0	0 0 0 0 0
% ];

coefvalues = {'f1'	...
     'HSIctr', 'HLIctr', 'HColctr',...
     'HSIhfd', 'HLIhfd', 'HColhfd',...
     'B1LIctr', 'B1LIhfd'};
%  A = [kmeanMatrix_joint(1,1)-kmeanMatrix_joint(1,2)	1	0 0 0	0	0	0 0 
% kmeanMatrix_joint(1,2)-kmeanMatrix_joint(1,3)	1	0 0 0	0	0	0	0 
% kmeanMatrix_joint(1,3)-kmeanMatrix_joint(1,4)	0	1	0 0 0	0	1	0 
% kmeanMatrix_joint(1,4)-kmeanMatrix_joint(1,5)	0	0 1 0	0	0	1 0
% kmeanMatrix_joint(1,5)-kmeanMatrix_joint(1,6)	0	0	0 0 0	0	0	0 
% kmeanMatrix_joint(2,1)-kmeanMatrix_joint(2,2)	1	0	0	0 0 0	0	0 
% kmeanMatrix_joint(2,2)-kmeanMatrix_joint(2,3)	1	0	0	0	0 0 0	0 
% kmeanMatrix_joint(2,3)-kmeanMatrix_joint(2,4)	0	1	0 0 0	0	0	0 
% kmeanMatrix_joint(2,4)-kmeanMatrix_joint(2,5)	0	0 1	0 0	0	0	0
% kmeanMatrix_joint(2,5)-kmeanMatrix_joint(2,6)	0	0	0 0 0	0	0	0 
% kmeanMatrix_joint(3,1)-kmeanMatrix_joint(3,2)	0	0	0 1	0 0	0	0
% kmeanMatrix_joint(3,2)-kmeanMatrix_joint(3,3)	0	0	0 1	0 0	0	0
% kmeanMatrix_joint(3,3)-kmeanMatrix_joint(3,4)	0	0	0	0 1 0 	0 1
% kmeanMatrix_joint(3,4)-kmeanMatrix_joint(3,5)	0	0	0	0 0 1	0	1
% kmeanMatrix_joint(3,5)-kmeanMatrix_joint(3,6)	0	0	0	0 0 0	0	0 
% kmeanMatrix_joint(4,1)-kmeanMatrix_joint(4,2)	0	0	0 1	0 0	0	0
% kmeanMatrix_joint(4,2)-kmeanMatrix_joint(4,3)	0	0	0 1	0 0	0	0
% kmeanMatrix_joint(4,3)-kmeanMatrix_joint(4,4)	0	0	0	0 1	 0 0	0 
% kmeanMatrix_joint(4,4)-kmeanMatrix_joint(4,5)	0	0	0	0 0 1	0	0 
% kmeanMatrix_joint(4,5)-kmeanMatrix_joint(4,6)	0	0	0	0	0	0 0 0 
% ];

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
% coefvalues = {'f1'	...
%      'HSIctr', 'HLIctr', 'HColctr', 'HFecctr',...
%      'HSIhfd', 'HLIhfd', 'HColhfd', 'HFechfd',...
%      'B1LIctr', 'B1LIhfd'};
%  
%  A = [kmeanMatrix_joint(1,1)-kmeanMatrix_joint(1,2)	1	0 0 0 0 0	0	0	0 0 
% kmeanMatrix_joint(1,2)-kmeanMatrix_joint(1,3)	1	0 0 0 0 0	0	0	0	0 
% kmeanMatrix_joint(1,3)-kmeanMatrix_joint(1,4)	0	1	0 0 0 0 0	0	1	0 
% kmeanMatrix_joint(1,4)-kmeanMatrix_joint(1,5)	0	0 1 0 0 0	0	0	0 0
% kmeanMatrix_joint(1,5)-kmeanMatrix_joint(1,6)	0	0 0 1 0 0   0	0	0	0 
% kmeanMatrix_joint(2,1)-kmeanMatrix_joint(2,2)	1	0	0 0 0	0 0 0	0	0 
% kmeanMatrix_joint(2,2)-kmeanMatrix_joint(2,3)	1	0	0 0 0	0	0 0 0	0 
% kmeanMatrix_joint(2,3)-kmeanMatrix_joint(2,4)	0	1	0 0 0 0 0	0	0	0 
% kmeanMatrix_joint(2,4)-kmeanMatrix_joint(2,5)	0	0 1	0 0 0 0	0	0	0
% kmeanMatrix_joint(2,5)-kmeanMatrix_joint(2,6)	0	0	0 1 0 0 0	0	0	0 
% kmeanMatrix_joint(3,1)-kmeanMatrix_joint(3,2)	0	0	0 0 1	0 0 0	0	0
% kmeanMatrix_joint(3,2)-kmeanMatrix_joint(3,3)	0	0	0 0 1	0 0 0	0	0
% kmeanMatrix_joint(3,3)-kmeanMatrix_joint(3,4)	0	0	0	0 0 1 0 0 	0 1
% kmeanMatrix_joint(3,4)-kmeanMatrix_joint(3,5)	0	0	0	0 0 0 1	0 0	0
% kmeanMatrix_joint(3,5)-kmeanMatrix_joint(3,6)	0	0	0	0 0 0 0 1	0	0 
% kmeanMatrix_joint(4,1)-kmeanMatrix_joint(4,2)	0	0	0 0 1 0	0 0	0	0
% kmeanMatrix_joint(4,2)-kmeanMatrix_joint(4,3)	0	0	0 0 1	0 0 0	0	0
% kmeanMatrix_joint(4,3)-kmeanMatrix_joint(4,4)	0	0	0	0 0 1 0 	 0 0	0 
% kmeanMatrix_joint(4,4)-kmeanMatrix_joint(4,5)	0	0	0	0 0 0 1 0 	0	0 
% kmeanMatrix_joint(4,5)-kmeanMatrix_joint(4,6)	0	0	0	0 0	0 0 1   0   0 
% ];

end