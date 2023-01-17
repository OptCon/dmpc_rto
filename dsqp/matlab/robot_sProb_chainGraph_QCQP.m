% MIT License
% 
% Copyright (c) 2023 Goesta Stomberg, Henrik Ebel, Timm Faulwasser, Peter Eberhard
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function sProb = robot_sProb_chainGraph_QCQP(Nrobot,N,T,x0,u0,xd,xinit)

import casadi.*

nx = 2; %states per robot
nu = 2; %inputs per robot

xmin = -inf(nx,1);
xmax = inf(nx,1);
umin = -0.2*ones(nu,1);
umax = 0.2*ones(nu,1);

Ad = eye(2);
Bd = T*eye(2);

% declare decision variables
for i=1:Nrobot
    % states
    XX{i}  = SX.sym(['x' num2str(i)], [nx N+1]);
    ZZZ{i} = cell(1,Nrobot);
    % state copies
    s = 1:Nrobot;
    for k=s(abs(s-i)==1)
       ZZZ{i}{k} =  SX.sym(['z' num2str(i) num2str(k)], [nx N+1]);
       llbz{i}{k} = -inf(nx,N+1);
       uubz{i}{k} = inf(nx,N+1);
    end
    % inputs
    UU{i}  = SX.sym(['u' num2str(i)], [nu N]);
    % inital condition
    XX0{i} = SX.sym(['xx0' num2str(i)], [nx 1]);
    UU0{i} = SX.sym(['uu0' num2str(i)], [nu 1]);
    % setpoint
    if i < Nrobot
        XXd{i} = SX.sym(['xxd' num2str(i)], [nx 2]);
    else
        XXd{i} = SX.sym(['xxd' num2str(i)], [nx 1]);
    end
    if i > 1
        slacks{i} = SX.sym(['s' num2str(i)], [1 1]); %slack for minimum-distance constraint
    end
    
end


% setup individual OCPs        
for i=1:Nrobot
    JJ{i}   = 0;
    gg{i}   = [];
    hh{i}   = [];
    llbx{i} = [];
    uubx{i} = [];
    llbu{i} = [];
    uubu{i} = [];

    % over horizon ...
    
    % initial condition
    gg{i} = [gg{i}; XX{i}(:,1) - XX0{i}];
    gg{i} = [gg{i}; UU{i}(:,1) - UU0{i}];
    % dynamics
    for j=1:N
        gg{i}      = [ gg{i}; XX{i}(:,j+1) - ( Ad*XX{i}(:,j) + Bd*UU{i}(:,j) )];
    end
    
    %terminal constraint
%     if i == 1
%         gg{i} = [gg{i}; XX{i}(:,N) - xd(:,i)];
%     else
%         gg{i} = [gg{i}; XX{i}(:,N)- (ZZZ{i}{i-1}(:,N) - [delta_x;delta_y] )];
%     end
    
    %box constraints
    for j= 1:N+1                            
    % state constraints
        llbx{i} = [ llbx{i}, xmin];
        uubx{i} = [ uubx{i}, xmax];
    end
    % input constraints
    for j = 1:N
        llbu{i} = [ llbu{i}, umin];
        uubu{i} = [ uubu{i}, umax];
    end

    %minimum distance constraint
    dmin = 0.4; %meter
    for j = 1:N+1
        if i > 1
            hh{i} = [hh{i}; (-(XX{i}(:,j) - ZZZ{i}{i-1}(:,j)).' * (XX{i}(:,j) - ZZZ{i}{i-1}(:,j)) + dmin^2) - slacks{i}];            
        end
    end
    if i > 1
        JJ{i} = JJ{i} + 10^6* slacks{i}.'*slacks{i};
    end
end

%cost
Q = 10*eye(2);
R = 1*eye(2);
beta = 10;%30;
P = beta*eye(2);

%stage cost
for i = 1:Nrobot
   ZZZ{i}{i} = XX{i};
   for k =  1:N
       if i == 1
           JJ{i} = JJ{i} + 0.5 * ( XX{i}(:,k) - XXd{i}(:,1) ).' * Q * ( XX{i}(:,k) - XXd{i}(:,1) );
       end

       if i > 1
           JJ{i} = JJ{i} + 0.25*( XX{i}(:,k) - ZZZ{i}{i-1}(:,k) - XXd{i}(:,1) ).'* Q * ( XX{i}(:,k) - ZZZ{i}{i-1}(:,k) - XXd{i}(:,1) );     
       end

       if i < Nrobot
           JJ{i} = JJ{i} + 0.25*( ZZZ{i}{i+1}(:,k) - XX{i}(:,k) - XXd{i}(:,2) ).'* Q * ( ZZZ{i}{i+1}(:,k) - XX{i}(:,k) - XXd{i}(:,2) );     
       end
       JJ{i} = JJ{i} + 0.5*UU{i}(:,k).'*R*UU{i}(:,k);
   end
end

%terminal cost
for i = 1:Nrobot
    k = N+1;
    if i == 1
       JJ{i} = JJ{i} + 0.5 * ( XX{i}(:,k) - XXd{i}(:,1) ).' * P * ( XX{i}(:,k) - XXd{i}(:,1) );
    end

    if i > 1
       JJ{i} = JJ{i} + 0.25*( XX{i}(:,k) - ZZZ{i}{i-1}(:,k) - XXd{i}(:,1) ).'* P * ( XX{i}(:,k) - ZZZ{i}{i-1}(:,k) - XXd{i}(:,1) );     
    end

    if i < Nrobot
       JJ{i} = JJ{i} + 0.25*( ZZZ{i}{i+1}(:,k) - XX{i}(:,k) - XXd{i}(:,2) ).'* P * ( ZZZ{i}{i+1}(:,k) - XX{i}(:,k) - XXd{i}(:,2) );     
    end
end



    
        
for i=1:Nrobot
   ZZZ{i}{i} = XX{i};
   llbz{i}{i} = llbx{i};
   uubz{i}{i} = uubx{i};
   llbzi = vertcat(vertcat(llbz{i}{:}));
   uubzi = vertcat(vertcat(uubz{i}{:}));
   
   ZZZi = [];
   for j = 1:Nrobot
       if ~isempty(ZZZ{i}{j})
           for k = 1:size(ZZZ{i}{j},2)
              ZZZi = [ZZZi; ZZZ{i}{j}(:,k)]; 
           end
       end       
   end
   UUi = [];
   for k = 1:size(UU{i},2)
       UUi = [UUi; UU{i}(:,k)];
   end
   XXU{i}    = [ ZZZi; UUi];
   llbxu{i} = [llbzi(:); llbu{i}(:)];
   uubxu{i} = [uubzi(:); uubu{i}(:)];
end

% set up consensus constraints
AA = cell(1,Nrobot);
for i = 1:Nrobot
    for j = setdiff(1:Nrobot,i)
        tmp = ZZZ{i}{j};
        if ~isempty(tmp) %i has copies of j
            ncopy = numel(tmp);
            nbefore  = numel(vertcat(vertcat(ZZZ{i}{1:j-1})));            
            nafter = numel(XXU{i})-ncopy-nbefore;
            AA{i} = [AA{i}; zeros(ncopy,nbefore),-eye(ncopy),zeros(ncopy,nafter)];
            
            nbefore = numel(vertcat(vertcat(ZZZ{j}{1:j-1})));
            nafter = numel(XXU{j})-ncopy-nbefore;
            AA{j} = [AA{j}; zeros(ncopy,nbefore),eye(ncopy),zeros(ncopy,nafter)];
            
            
            for k = setdiff(1:Nrobot,[i,j])
                AA{k} = [AA{k}; zeros(ncopy,numel(XXU{k}))];
            end
                
            clear ncopy ZZZi nbefore nafter tmp
        end        
    end

    pp{i} = [XX0{i};UU0{i};vertcat(XXd{i}(:))];
end

%account for slack variables in AA
for i = 2:Nrobot
    XXU{i} = [XXU{i};vertcat(slacks{i}(:))];
    AA{i} = [AA{i},zeros(size(AA{i},1),1)];
    llbxu{i} = [llbxu{i}; 0];
    uubxu{i} = [uubxu{i}; inf];
end



% set up sProb
% convert expressions to MATLAB functions
for i=1:Nrobot    
    sProb.locFuns.ffi{i} = Function(['f' num2str(i)],{XXU{i},pp{i}},{JJ{i}});
    sProb.locFuns.ggi{i} = Function(['g' num2str(i)],{XXU{i},pp{i}},{gg{i}});
    sProb.locFuns.hhi{i} = Function(['h' num2str(i)],{XXU{i},pp{i}},{hh{i}});
    
    % set up dIP parameters
    sProb.llbx{i}  = llbxu{i};
    sProb.uubx{i}  = uubxu{i};
    sProb.AA{i}    = AA{i};
    
    if isempty(xinit)
        sProb.zz0{i}   = zeros(length(XXU{i}),1);
    else
        sProb.zz0{i} = xinit{i};
    end
    
    sProb.llam0{i} = zeros(size(AA{i},1),1);
    sProb.pp{i} = [x0{i};u0{i};vertcat(xd{i}(:))];

end

end