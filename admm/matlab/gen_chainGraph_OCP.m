% MIT License
% Copyright (c) 2023 Goesta Stomberg, Henrik Ebel, Timm Faulwasser, Peter Eberhard
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

close all
clear all
clc

import casadi.*

% scenario
Nrobot = 4;
N = 7;          %horizon
Ts = 0.2;       %sampling interval seconds
xx0 = {[0.5;0],[0.15;0.3],[-0.15;-0.3],[-0.5;0]};
uu0 = {[0;0],[0;0],[0;0],[0;0]};

%setpoints
x1d = [0.7; -0.7]; %desired setpoint for x1
x2d = [-0.4; 0.0];  %desired setpoint for (x2-x1)
x3d = [-0.4; 0.0];  %desired setpoint for (x3-x2)
x4d = [-0.4; 0.0];  %desired setpoint for (x4-x3)

xxd{1} = [x1d, x2d];
xxd{2} = [x2d, x3d];
xxd{3} = [x3d, x4d];
xxd{4} = [x4d];

xinit = []; %solver initialization

% setup OCP
sProb = robot_sProb_chainGraph(Nrobot,N,Ts,xx0,uu0,xxd,xinit);

%export C code
gen_c_sProb(sProb);

