function [gitfit] = fitGITmodel(kmeanMatrix_joint, ncond, shuffle_flag)

fitwarning = '';
% minimal value to replace 0 to avoid numerical errors
minval = 0.01;
% number of times to repeat optimization with trust-region-reflective method
nrepeat = 100;

% replace small values with noise
kmeanMatrix_joint_orig = kmeanMatrix_joint;
kmeanMatrix_joint(kmeanMatrix_joint<minval)=...
kmeanMatrix_joint(kmeanMatrix_joint<minval)+rand(nnz(kmeanMatrix_joint<minval),1)*minval;


% calculate fluxes for average profile
Aorig = calculateAmatrix_final(kmeanMatrix_joint_orig);
[A,coefvalues] = calculateAmatrix_final(kmeanMatrix_joint);
        
% check whether there are 2 conditions (GF and Coonized) or 4 conditions
% (GF and Colonized with two different diets)
if ncond==2
    % half the matrices
    Aorig = Aorig(1:size(Aorig,1)/2,:);
    A = A(1:size(A,1)/2,:);
    % remove columns for second diet
    Aorig(:,[5:7 9]) = [];
    A(:,[5:7 9]) = [];
    % end half the matrices
end
        
A(Aorig(:,1)==0,:)=[];
% calculate vector b and upper and lower bounds on x
b = zeros(size(A,1),1);

xlowerlim = zeros(size(A,2),1);
xlowerlim(2:end)=-Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%
xupperlim = Inf*ones(size(A,2),1);
%%%%%%%%%%%%%%%%%%%%%%%%%
% check if A is empty an if yes move to the next metabolite
if isempty(A)
    % there are no values - return epmty fit
    gitfit = [];
    return
end

% run trust-region-reflective algorith multiple times with different initial solutions
testx = zeros(size(A,2),nrepeat);
%testxres = zeros(1,nrepeat);
testCorrRev = zeros(1,nrepeat);
testCorrRevSI = zeros(1,nrepeat);
testCorrRevLI = zeros(1,nrepeat);
testCorrRev_shuffled = zeros(1,nrepeat);
% store restored data for each solution in a matrix
testdataR = zeros(numel(kmeanMatrix_joint), nrepeat);
testx_idx=1;

% test nrepeat solutions from trust-region-reflective
if size(A,2)<=size(A,1)
    for nrepeat_i = 1:nrepeat
        options = optimoptions(@lsqlin,'Display', 'off',...
            'Algorithm','trust-region-reflective','MaxIterations',1500);
        rand_x0=rand(size(A,2),1);
% obsolete: if repeat the diet twice, set initial solution equal to each
% other for two diets
%                 rand_x0=rand(5,1);
%                 rand_x0 = [rand_x0(1); rand_x0(2:4);rand_x0(2:4);rand_x0(5);rand_x0(5)];
%                 if testx_idx==10
%                     rand_x0 = testx(:,1); % interior point solution
%                 end
        [x] = lsqlin(A,b,[],[],[],[],xlowerlim, xupperlim, rand_x0, options);
        testx(:, testx_idx) = x;
        %testxres(testx_idx) = xres;

        % if there are two conditions, they have to be 
        % duplicated to solve the reverse problem for now
        % because it expects four conditions
        if ncond==2 
            x = [x(1); x(2:4); x(2:4); x(5); x(5)];
        end
        % calculate reverse problem
        % use interior point since number of variables is larger than the
        % number of equations
        options = optimoptions(@lsqlin,'Display', 'off',...
            'Algorithm','interior-point','MaxIterations',1500);
        % mulriply solution x to 1000 to avoid numerical issues
        [Ra,rb] = calculateRAmatrix_final(x*1000);
        dataR = lsqlin(Ra,rb,[],[],[],[],zeros(1,size(Ra,2)), [], [], options);
        dataR = reshape(dataR,[],4)';
        % calculate correlations between original data and restored data
        % for the whole GIT profile
        testCorrRev(testx_idx) = corr(kmeanMatrix_joint_orig(:), dataR(:));
        % separately for small intestine
        testCorrRevSI(testx_idx) = corr(reshape(kmeanMatrix_joint_orig(:,1:3),[],1),...
                         reshape(dataR(:,1:3),[],1));
        % separately for large intestine
        testCorrRevLI(testx_idx) = corr(reshape(kmeanMatrix_joint_orig(:,4:end),[],1),...
                         reshape(dataR(:,4:end),[],1));
        % save restored dataR in a matrix as well
        testdataR(:,testx_idx) = dataR(:);

        if shuffle_flag
        % calculate shuffled results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate reverse problem for shuffled x
            x_shuffled = x(randperm(length(x)));
            [Ra_shuffled,rb_shuffled] = calculateRAmatrix_final(x_shuffled*1000);
            dataR_shuffled = lsqlin(Ra_shuffled,rb_shuffled,[],[],[],[],...
                                    zeros(1,size(Ra_shuffled,2)), [], [], options);
            dataR_shuffled = reshape(dataR_shuffled,[],4)';

            testCorrRev_shuffled(testx_idx) = corr(kmeanMatrix_joint_orig(:), ...
                dataR_shuffled(:));
        end

        testx_idx = testx_idx+1;
    end
else
     fitwarning = [fitwarning ...
         'Warning: number of variables > number of equations, '...
         ' only interior-point algorithm can be used\n'];
end

% calculate solution with interior-point algorithm
options = optimoptions(@lsqlin,'Display', 'off',...
            'Algorithm','interior-point','MaxIterations',1500);
rand_x0=rand(size(A,2),1);
[x_ip] = lsqlin(A,b,[],[],[],[],xlowerlim, xupperlim, rand_x0, options);
x_ip_orig = x_ip;
 
% if there are two conditions, they have to be 
% duplicated to solve the reverse problem for now
% because it expects four conditions
if ncond==2 
    x_ip = [x_ip(1); x_ip(2:4); x_ip(2:4); x_ip(5); x_ip(5)];
end
% calculate reverse problem
% use interior point since number of variables is larger than the
% number of equations
options = optimoptions(@lsqlin,'Display', 'off',...
    'Algorithm','interior-point','MaxIterations',1500);
% mulriply solution x to 1000 to avoid numerical issues
[Ra,rb] = calculateRAmatrix_final(x_ip*1000);
dataR_ip = lsqlin(Ra,rb,[],[],[],[],zeros(1,size(Ra,2)), [], [], options);
dataR_ip = reshape(dataR_ip,[],4)';
% calculate correlations between original data and restored data
% for the whole GIT profile
x_ip_CorrRev = corr(kmeanMatrix_joint_orig(:), dataR_ip(:));
% separately for small intestine
x_ip_CorrRevSI = corr(reshape(kmeanMatrix_joint_orig(:,1:3),[],1),...
                 reshape(dataR_ip(:,1:3),[],1));
% separately for large intestine
x_ip_CorrRevLI = corr(reshape(kmeanMatrix_joint_orig(:,4:end),[],1),...
                 reshape(dataR_ip(:,4:end),[],1));
             
% return gitfit structure with solutions
gitfit.kmeanMatrix_joint_orig = kmeanMatrix_joint_orig;
gitfit.A = A;
gitfit.b = A;
gitfit.xlowerlim = xlowerlim;
gitfit.xupperlim = xupperlim;

gitfit.x_ip = x_ip_orig;
gitfit.dataR_ip = dataR_ip;
gitfit.x_ip_CorrRev = x_ip_CorrRev;
gitfit.x_ip_CorrRevSI = x_ip_CorrRevSI;
gitfit.x_ip_CorrRevLI = x_ip_CorrRevLI;

% add vector of trust region solutions
gitfit.testx = testx;
gitfit.testdataR = testdataR;
gitfit.testCorrRev = testCorrRev;
gitfit.testCorrRevSI = testCorrRevSI;
gitfit.testCorrRevLI = testCorrRevLI;
gitfit.testCorrRev_shuffled = testCorrRev_shuffled;

% save coefvalues returned by the calculate A matrix function
gitfit.coefvalues = coefvalues;
gitfit.fitwarning = fitwarning;
