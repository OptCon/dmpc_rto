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

function gen_c_sProb( sProb )

import casadi.*

NsubSys = size(sProb.AA,2);

opts = struct('main', false,...
              'mex', false);

C = CodeGenerator('locFuns.c',opts);

%% sensitivities, HessL, Phi
Ncons   = size(sProb.AA{1},1); 
for i=1:NsubSys
    nx{i} = size(sProb.AA{i},2);
    xCas{i} = SX.sym(['x',num2str(i)],nx{i},1);
    pCas{i} = SX.sym(['p',num2str(i)],size(sProb.pp{i}));
    

% separate llbx, uubx from Hcas and Gas
    Hcas{i} = sProb.locFuns.hhi{i}(xCas{i},pCas{i});
    Gcas{i} = sProb.locFuns.ggi{i}(xCas{i},pCas{i});
           
    Fcas{i}  = sProb.locFuns.ffi{i}(xCas{i},pCas{i});

    ng{i} = size(Gcas{i},1);
    nh{i} = size(Hcas{i},1);
    Cond.ng{i} = ng{i};
    Cond.nh{i} = nh{i};
    nuCas{i} = SX.sym(['nu',num2str(i)],ng{i},1);
    muCas{i}  = SX.sym(['mu',num2str(i)],nh{i},1);
    lamCas{i} = SX.sym('lam',Ncons,1);


    gradFCas{i} = jacobian(Fcas{i}, xCas{i});
    JGCas{i}    = jacobian(Gcas{i}, xCas{i});
    JHCas{i}    = jacobian(Hcas{i}, xCas{i});

    HessFCas{i} = jacobian(gradFCas{i}, xCas{i});
    
    gradLCas{i} = gradFCas{i} + nuCas{i}'*JGCas{i} + muCas{i}'*JHCas{i};
    HessLCas{i} = jacobian(gradLCas{i}, xCas{i});


    %d-SQP stopping criterion
    PhiCas{i} = [   gradLCas{i}.';
                    Gcas{i};];
    gradPhiCas{i} = [jacobian(PhiCas{i}, xCas{i}), jacobian(PhiCas{i}, nuCas{i})]; 

    gradF{i} = Function(['gradFun',num2str(i)],{xCas{i},pCas{i}},{gradFCas{i}});
    JG{i}    = Function(['JGfun',num2str(i)],{xCas{i},pCas{i}},{JGCas{i}});
    JH{i}    = Function(['JHfun',num2str(i)],{xCas{i},pCas{i}},{JHCas{i}});
    HessF{i} = Function(['HessFfun',num2str(i)],{xCas{i},pCas{i}},{HessFCas{i}});
    eq{i}    = Function(['eqfun',num2str(i)],{xCas{i},pCas{i}},{Gcas{i}});
    ineq{i}  = Function(['ineqfun',num2str(i)],{xCas{i},pCas{i}},{Hcas{i}});
    
    HessL{i} = Function(['HessLfun',num2str(i)],{xCas{i},nuCas{i},muCas{i},pCas{i}},{HessLCas{i}});

    Phi{i}   = Function(['Phi',num2str(i)],{xCas{i},nuCas{i},muCas{i},lamCas{i},pCas{i}},{PhiCas{i}});
    JPhi{i}  = Function(['JPhi',num2str(i)],{xCas{i},nuCas{i},muCas{i},lamCas{i},pCas{i}},{gradPhiCas{i}});



    C.add(gradF{i});
    C.add(JG{i});
    C.add(JH{i});
    C.add(HessF{i});
    C.add(eq{i});
    C.add(ineq{i});
    C.add(HessL{i});
    C.add(Phi{i});
    C.add(JPhi{i});
end

C.generate();

for i = 1:NsubSys
    str = sprintf('A%i.csv',i);
    writematrix(sProb.AA{i},str);

    try
        str = sprintf('ub%i.csv',i);
        writematrix(sProb.uubx{i},str);
        str = sprintf('lb%i.csv',i);
        writematrix(sProb.llbx{i},str);
    catch
    end

end

end