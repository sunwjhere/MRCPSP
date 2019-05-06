classdef Library
    methods(Static)
        

        %% ================================================================
        % function optimizationMRCPSP
        % updated repair schedule using the formulation of multi-mode 
        % resource-constrained project scheduling problem (MRCPSP)
        %Robert Klein (2000), Scheduling of resource-constrained projects, Kluwer
        %Academic Publishers, page 97-98.  
        % ================================================================
        function [CT,SFT,Ur, Unr] = optimizationMRCPSP(task, precedence, resourceRenewableConstraint, resourceNonRenewableConstraint, time_horizon)
        %=================================================================
        % Input variables:
        % 1.task = cells of restoration task data, 
        %   -1.1 the number of task mode = the mode of cells in "task"
        %   -1.2 the number of actual tasks = the number of rows in every cell of "task"
        %   -1.3 the duration of every task at a mode = the 1st column in every cell of "task"
        %   -1.4 the resource demand(s) at a mode = the 2nd column - the end column in every cell of "task"                         
        % 2.precedence(P)=binary precedence matrix
        % 3.resourceRenewableConstraint (R) = a matrix of renewable resource
        %     constraint(s) over time 
        %     [either constaint value(s) or varying values over time]
        % 4.resourceNonRenewableConstraint (NR): a matrix of non-renewable
        % resource constraint(s) [constaint value(s)]
        % 5.time_horizon (H) = time horizon of interest
        %
        % Output variables:
        % 1.CT = completion time
        % 2.SFT = [task starting time, task finishing time, selcted task
        % mode, task duration for the selcted mode]
        % 3.Ur = a matrix of used amount of renewable resource(s) over time
        % 4.Unr =  a matrix of used amount of non-renewable resource(s)
        %
        % Copyright 2019 Wenjuan Sun
        %=================================================================
            %======================= input parameters===========================
            H = time_horizon;                       % time horizon
            R = resourceRenewableConstraint;        % renewable resource constraints
            NR = resourceNonRenewableConstraint;    % nonrenewable resource constraints
            Kr = size(R,2);                         % the number of renewable resource type
            Knr = size(NR,2);                       % the number of nonrenewable resource type
            Nt = H;                                 % the number of time steps
            m = size(task,2);                       % the number of task modes
            t = linspace(1,Nt,Nt);
           
            % add a dummy end task
            I = size(task{1},1)+1;                  % the number of task (add a dummy end task)
            P = [precedence; ones(1,I-1)]; P =[P, zeros(I,1)];  % precedence

            % assign duration d and resource requirement r for every task 
            % during schedule, task duration d = duration mode
            
            for im = 1:m
                for ii = 1:I-1
                    d{im}(ii,1) = max(round(task{im}(ii, 1)), 1); % make sure the duration greater than 0
                    for ir = 1:Kr
                        r{im}(ii,ir) = task{im}(ii, 1+ir);
                    end
                    for inr = 1:Knr
                        nr{im}(ii,inr) = task{im}(ii, (1+Kr)+inr);
                    end
                end
                ii = I; % dummy end task
                d{im}(ii,1) = 0; r{im}(ii,:) = zeros(1,Kr); nr{im}(ii,:) = zeros(1,Knr);  
            end
             
            dd = []; rr = []; rnr = []; 
            for im = 1:m
                dd = vertcat(dd,d{im});
                rr = vertcat(rr,r{im});
                rnr = vertcat(rnr,nr{im});
            end

            z = linspace(1,I,I); 
            
            % compute the time bounds as ELFT(EFT&LFT) for every task considing
            % both fast mode(s) and slow mode(s)
            [ELSTm, ELFTm, CFTm] = Library.CPMmm(d,P,H);
            
            % index of low bound and upper bound for the time variable  
            % from time to time index (0->1, 1->2, ...)
            for im = 1:m
                lbm{im} = ELFTm{im}(:,1);
                ubm{im} = ELFTm{im}(:,2);
            end
            
            %% using "yalmip" to construct the optimization problem
            % ====================== decision variables===========================
            x = binvar(m*I,Nt); % which time step task i(=1,2,...,I,I+1) is finished within [0, H]


            % =========================constraints===============================
            constraints = [];
            
            % c1: every task only executes once. [equation(3.70)] 
            % considering t \in [0, H], t(time_index+1)=time_index: t(1)=0, t(H+1) = H;
    
            for ii = 1:I   
                for im = 1:m
                    lb1 = lbm{im}(ii);
                    ub1 = ubm{im}(ii);
                    iinew = ii + (im-1)*I;
                    tmp1(ii,im) = sum(x(iinew, lb1:ub1)); 
                end
                
                constraints = [constraints, sum(tmp1,2) - 1 == 0]; 
            end

            % c2: precedence
            % [equation(3.71)] 
            tmp21 = x.*repmat(t,m*I,1);
            tmp22 = x.*(repmat(t,m*I,1)-repmat(dd,1,Nt)); 
            for ii = 1:I
                for im = 1:m
                    lb = lbm{im}(ii);
                    ub = ubm{im}(ii);
                    iinew = ii + (im-1)*I;
                    tmppre(ii,im) = sum(tmp21(iinew,lb:ub)); 
                    tmpsuc(ii,im) = sum(tmp22(iinew,lb:ub)); 
                end
            end

            tpre = sum(tmppre,2);               % finishing time of every task i in Prej
            tsuc = sum(tmpsuc,2);               % starting time of every task j
            
            for ii = 1:I
                J = P(ii,:).*z;                 %predecesor task of task ii
                J1 = nonzeros(J);
                if ~isempty(J1)
                    for ipre = 1:length(J1)     % for every precedence task idx = J1(ipre) for task ii
                        idx = J1(ipre);
                        constraints=[constraints, tsuc(ii)-tpre(idx) >= 0];
                    end
                end
            end
                        
            % c3: renewable resource usage should be less than the resource constraint for every resource type and for every time step
            % [equation(3.72)] 
            if Kr>0 % if there is renewable resource
                ub3 = []; lb3 = [];  tmpub =[]; tmplb =[];
                for im = 1:m
                    for it = 1:Nt
                        tmpub{im}(:,it) = min([it-1+d{im}, ubm{im}],[],2); 
                        tmplb{im}(:,it) = max([it+zeros(I,1), lbm{im}],[],2);   
                    end
                    ub3 = vertcat(ub3,tmpub{im});
                    lb3 = vertcat(lb3,tmplb{im});
                end


                % the following implementation works when the decision variable 
                % x=1 representing the finish of the task!

                for it = 1:Nt
                    for ii = 1:I*m 
                        lb = lb3(ii,it);
                        ub = ub3(ii,it);
                        if ge(ub,lb)==1
                            tmpx(ii,it) = sum(x(ii,lb:ub));
                        end
                    end
                end

                tmpx = [tmpx; zeros(1,Nt)]; 

                for it = 1:Nt
                    for ir = 1:Kr
                        tmpr = [];
                        tmpr = rr(:,ir).*tmpx(:,it);
                        ru = sum(tmpr);
                        constraints=[constraints, R(it,ir)-ru >= 0];
                    end
                end
            
            end
            
           % c4: nonrenewable resource usage
           if Knr > 0 % if there is nonrenewable resource consumed
               ub4 = []; lb4 = []; 
               
               for im = 1:m
                   ub4 = vertcat(ub4,ubm{im});
                   lb4 = vertcat(lb4,lbm{im});
               end

               for ii = 1:m*I
                   lb = lb4(ii);
                   ub = ub4(ii);
                   if ge(ub,lb)==1
                       tmpxn(ii,1) = sum(x(ii,lb:ub));
                   end
               end
               
               
               for ir = 1:Knr
                   tmprn = rnr(:,ir).*tmpxn;
                   rnu = sum(tmprn);
                   constraints = [constraints, NR(ir)-rnu >= 0];
               end
           end
             
           % =========================objective=============================
            % minimize the finishing time of the dummy end task
            objective = 0;
            
            for ii = 1:I
                for im = 1:m
                    lb1 = lbm{im}(ii);
                    ub1 = ubm{im}(ii);
                    iinew = ii + (im-1)*I;
                    tmpo = t.*x(iinew,:);
                    ft(iinew) = sum(tmpo(lb1:ub1));
                end
            end
            objective = objective + max(ft); 


           % ==========================solve=================================
            %ops = sdpsettings('solver','gurobi','gurobi.timelimit',200,'verbose',0,'showprogress',1);
            ops = sdpsettings('solver','gurobi','verbose',0,'showprogress',1);
            optimize(constraints,objective,ops)
            opx = value(x); % value of variable x 
            obj = value(objective); % value of objective function

            % ==========================post-process=========================
            X = opx;     X(isnan(X))=0;
            
            % the finishing time of every task
            tX = zeros(I,1);
            for im = 1:m
                ii = (im-1)*I+1;
                jj = im*I;
                tX = tX + X(ii:jj,:);
            end
            ft = sum(tX.*repmat(t,I,1),2);
            
            % track which mode a task has used
            M = [];
            for im = 1:m
                tmpm = im+zeros(I,1);
                M = [M; tmpm]; 
            end
            
            tmpM = reshape(sum(X.*repmat(M,1,Nt),2),I,m); 
            trackM = sum(tmpM,2);                   % the mode number of eevry task
            %ratioM = nnz(trackM-1)/length(trackM);  % ratio of the amount of tasks using mode 1 among all tasks
            
            if ~all(trackM)
                disp('The optimization solution is infeasbible or wrong. Terminate the function.')       
                CT=[]; SFT=[]; Ur=[]; Unr=[];
                return
            end
            

            % st: starting time of every task
            % rm: required resource amount in the optimal result (a mixture of resource requirements from different modes)
            for ii = 1:I
                im = trackM(ii);
                dm(ii,1) = d{im}(ii);
                rrm(ii,:) = r{im}(ii,:); 
                rnrm(ii,:) = nr{im}(ii,:); 
                st(ii,1) = ft(ii)-dm(ii);     
            end
            %ttmp = [st ft cell2mat(d) trackM dm]
            
            % ru: actual used resoruce amount over every time step for every resource type 
            XR = zeros(I,Nt);
            for ii = 1:I
                ilb = st(ii)+1;
                iub = ft(ii);
                blength = dm(ii);
                XR(ii, ilb:iub) = ones(1,blength);
            end
            Ur = rrm'*XR; 
            Unr = sum(rnrm);
            %Unr = rnrm'*XR; 

            SFT = [st,ft, trackM, dm];    % starting time (1st column) and finishing time (2nd column)
            CT = max(ft);       % completion time of all tasks

        end
        
        %% ================================================================
        % function CPMsm
        % get EFT&LFT for every task in scheduling using critial path method
        % (CPM) for scheduling single-mode tasks 
        % =================================================================
        function [ELST, ELFT, CFT] = CPMsm(d,p,th)
        %=================================================================
        % Input variables:
        % 1.d = task duration vector
        % 2.p = precedence matrix
        % 3.th = time horizon
        %
        % Output variables:
        % 1.ELST = [task earliest starting time, task latest starting time]
        % 2.ELFT = [task earliest finishting time, task latest finishing time]    
        % 3.CFT = completion time of all tasks   
        %
        % Copyright 2019 Wenjuan Sun
        %=================================================================
            Nt = size(d,1);
            
            ES = zeros(Nt,1);
            EC = ES+d;
            LC = th+zeros(Nt,1);
            LS = LC-d;
            
            % Given a task #it, find its all precence task(s)# 
            % the starting time of this task >= the finishing time of 
            % its predecessor task(s)
            % the finishing time of this task = starting time + its duration 
            for it = 1:Nt
                ipt = find(p(it,:));
                if ~isempty(ipt)
                    for ii = 1:length(ipt)
                        idx = ipt(ii); 
                        ES(it,1) = max([EC(idx),ES(it,1)]);
                    end
                    EC(it,1) = ES(it,1)+d(it);
                    
                end
            end
            
            % is LC>=EC?
            % if not (LC<EC), display error
            tmp = lt(LC-EC,0);
            check = find(tmp);
            if ~isempty(check) 
                msg = strcat('Function CPM: completion time error (ect > lct) for task', num2str(check));
                disp(msg); 
            else % LC>=EC
                % output
                ELST = [ES,LS];
                ELFT = [EC,LC];
                CFT = max(ELFT);
            end
            ELST = ELST; ELFT = ELFT; CFT = CFT;
        end
   
        %% ================================================================
        % function CPMmm
        % get EFT&LFT for every task in scheduling using critial path method
        % (CPM) for scheduling multi-mode tasks 
        % =================================================================
        function [ELSTm, ELFTm, CFTm] = CPMmm(d,p,th) 
        %=================================================================
        % Input variables:
        % 1.d = task duration vector
        % 2.p = precedence matrix
        % 3.th = time horizon
        %
        % Output variables:
        % 1.ELSTm = [task earliest starting time, task latest starting time]
        % 2.ELFTm = [task earliest finishting time, task latest finishing time]    
        % 3.CFTm = completion time of all tasks   
        %
        % Copyright 2019 Wenjuan Sun
        %=================================================================

            m = size(d,2);
            for im = 1:m               
                [ELSTm{im}, ELFTm{im}, CFTm{im}] = Library.CPMsm(d{im},p,th);
            end
            ELSTm = ELSTm;
            ELFTm = ELFTm;
            CFTm = CFTm;
        end

        
    end
end
