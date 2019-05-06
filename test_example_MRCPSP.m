%%%%% test_example_MRCPSP.m

close all
clear
clc


%== an example CPM graph and its input data format
% % (0)->(1)->(4)->(6)
% %     \(2)->(5)/
% %     \(3)/
% % (0)-dummy start
% % (6)-dummy end

%%%==========prepare input data
% multi mode
th = 40; % time horizon
Nt = th; % number of time step
 
%== task duration
% tmpdd: task duration in two modes
tmpdd{1} = [4,2,3,2,2]'; 
tmpdd{2} = [6,2,5,2,4]';

%== resource demand
% tmprr: renewable resource demand of every task in two modes
tmprr{1} = [2,1,3,1,2]'; 
tmprr{2} = [1,1,1,1,1]';

% tmpnrr: non-resource demand of every task in two modes
tmpnrr{1} = [3,2,4,2,3]';
tmpnrr{2} = [1,1,2,1,2]';    
        
I = length(tmpdd{1});   % number of task
m = length(tmpdd);      % number of task mode

%== p: binary precedence matrix
% if the icol is the precedensor task of the irow, p(irow,icol) = 1
% otherwise p(irow,icol) = 0
p = zeros(I,I);
p(4,1) = 1;
p(5,2) = 1;
p(5,3) = 1;        
        
%== define resource constraints
ar = 4;   % renewable resource constraint
anr = 8;  % non-renewable resource constraint
  
% constaint non-renewable resource constraint
resourceNR = anr;
    
% constant renewable resource constraint over time
resourceRConstant = repmat(ar,Nt,1); 

% renewable resource varies over time?
resourceRVary = resourceRConstant;
tvary = linspace(2,4,3);
rvary = 2*ar;
        
for it = 1:length(tvary)
    idx = tvary(it);
    resourceRVary(idx,:) = rvary; 
end
                
%== construct the cell of "task"
for ii = 1:m
    tmptask{ii} = [tmpdd{ii},tmprr{ii},tmpnrr{ii}]; 
end


%% assign input to use teh MRCPSP optimization 
task = tmptask;
time_horizon = th;
precedence = p;
resourceRenewable = resourceRConstant;
% resourceRenewable = resourceRVary; 
resourceNonRenewable = resourceNR;

[CT,SFT,Ur, Unr] = Library.optimizationMRCPSP(task, precedence, resourceRenewable, resourceNonRenewable, time_horizon);


%% post-process and plot figures
st = SFT(1:I,1);
ft = SFT(1:I,2);
trackM = SFT(1:I,3); 
dm = SFT(1:I,4); 

figure('Renderer', 'painters', 'Position', [20 50 1000 600])

% starting time and ending time
subplot(2,3,1),
barh(ft);
hold on, barh(st,'FaceColor','w',...
       'EdgeColor','w');
   grid on;
xlabel('t');
ylabel('task ID');
title('Starting Time and Finishing Time');

% task mode selected
subplot(2,3,2),
barh(trackM);
xlim([0,m+1]);
ylabel('task ID');
xlabel('mode');
title('Task Mode Adopted');

% task duration in the corresponding mode
subplot(2,3,3),
barh(dm);
xlabel('duration');
ylabel('task ID');
title('Duration Adopted');


% renewable resource usage
tx = linspace(1,Nt,Nt)-0.5;

ir = 1;
subplot(2,3,4), bar(Ur(ir,:));
hold on, stairs(tx,resourceRenewable(:,1));
ylim([0,10]);
xlabel('t');
ylabel('Resource usage');
str = strcat('Renewable Resource Type ', num2str(ir));
title(str);


% non-renewable resource usage
subplot(2,3,5), bar(Unr);
hold on, line([0 th],[resourceNR resourceNR]);
xlim([0,3]);
ylim([0,15]);
ylabel('Resource usage');
title('Non-renewable Resource')




        
        