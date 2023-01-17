/*MIT License
Copyright (c) 2023 Goesta Stomberg, Henrik Ebel, Timm Faulwasser, Peter Eberhard
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <thread>
#include <lcm/lcm-cpp.hpp>
#include <chrono>

#include "sProb.hpp"
#include "admmAgent.hpp"
#include "dmpcAgent.hpp"

void lcm_handler_threadfunction(lcm::LCM* lcm, bool* finished) {
    while (!(*finished)) {
        lcm->handleTimeout(1000);
    }
}


int main(int argc, char *argv[]) {

    int32_t my_id = atoi(argv[1]);

    lcm::LCM lcm("udpm://239.255.76.67:7667?ttl=1");
    //lcm::LCM lcm;
    if (!lcm.good()) {
        return 1;
    }
    // start LCM listener / handler thread
    bool finished = false;
    std::thread lcm_handling_thread(lcm_handler_threadfunction, &lcm, &(finished));


    
    std::cout << "My ID is " << my_id << std::endl; 

    unsigned int Nagents = 4;
    std::string folderName = "sProb_chain_QP_N7";

    sProb mysProb = sProb(Nagents);
    mysProb.read_sProb(folderName,Nagents);
    
    //Setup meta data specifically for chain-coupling
    mpcMetaData meta_;

    meta_.nx = 2;
    meta_.nu = 2;
    meta_.idx_eqx0 = 0;
    meta_.idx_equ0 = 2;
    int N = 7; //MPC horizon

    if (my_id == 0 || my_id == Nagents - 1){
        meta_.idx_u0 = 2*(N+1)*meta_.nx;
        meta_.idx_u1 = 2*(N+1)*meta_.nx + meta_.nu;
        meta_.idx_u2 = 2*(N+1)*meta_.nx + 2*meta_.nu;

        //x1
        meta_.shift_seq_start.push_back(0);
        meta_.shift_seq_length.push_back((N+1)*meta_.nx);
        meta_.shift.push_back(meta_.nx);

        //x2
        meta_.shift_seq_start.push_back((N+1)*meta_.nx);
        meta_.shift_seq_length.push_back((N+1)*meta_.nx);
        meta_.shift.push_back(meta_.nx);

        //u1
        meta_.shift_seq_start.push_back(2*(N+1)*meta_.nx);
        meta_.shift_seq_length.push_back(N*meta_.nu);
        meta_.shift.push_back(meta_.nu);

    }else{
        meta_.idx_u0 = 3*(N+1)*meta_.nx;
        meta_.idx_u1 = 3*(N+1)*meta_.nx + meta_.nu;
        meta_.idx_u2 = 3*(N+1)*meta_.nx + 2*meta_.nu;

        //x1
        meta_.shift_seq_start.push_back(0);
        meta_.shift_seq_length.push_back((N+1)*meta_.nx);
        meta_.shift.push_back(meta_.nx);

        //x2
        meta_.shift_seq_start.push_back((N+1)*meta_.nx);
        meta_.shift_seq_length.push_back((N+1)*meta_.nx);
        meta_.shift.push_back(meta_.nx);

        //x3
        meta_.shift_seq_start.push_back(2*(N+1)*meta_.nx);
        meta_.shift_seq_length.push_back((N+1)*meta_.nx);
        meta_.shift.push_back(meta_.nx);

        //u1
        meta_.shift_seq_start.push_back(3*(N+1)*meta_.nx);
        meta_.shift_seq_length.push_back(N*meta_.nu);
        meta_.shift.push_back(meta_.nu);
    }

    VectorXd p;
    if (my_id < Nagents-1){
        p = Eigen::VectorXd::Zero(8); //p = [x0; u0; xd]; xd = [x_my_id_d; x_my_id+1_d];
        meta_.nxd = 4;
    } else {
        p = Eigen::VectorXd::Zero(6); //p = [x0; u0; xd]; xd = [x_my_id_d];
        meta_.nxd = 2;
    }


//initial positions
    Vector2d x10; x10 << 0.5, 0.0;
    Vector2d x20; x20 << 0.15, 0.3;
    Vector2d x30; x30 << -0.15, -0.3;
    Vector2d x40; x40 << -0.5, 0.0;

//initial velocities
    Vector2d u10; u10 << 0.0, 0.0;
    Vector2d u20; u10 << 0.0, 0.0;
    Vector2d u30; u10 << 0.0, 0.0;
    Vector2d u40; u10 << 0.0, 0.0;

//desired positions
    Vector2d x1d; x1d << 0.7, -0.7;
    Vector2d x2d; x2d << -0.4, 0.0;
    Vector2d x3d; x3d << -0.4, 0.0;
    Vector2d x4d; x4d << -0.4, 0.0;

    if (my_id == 0){
        p[0] = x10[0]; p[1] = x10[1];
        p[2] = u10[0]; p[3] = u10[1];
        p[4] = x1d[0]; p[5] = x1d[1];
        p[6] = x2d[0]; p[7] = x2d[1]; 
    } else if (my_id == 1){
        p[0] = x20[0]; p[1] = x20[1];
        p[2] = u20[0]; p[3] = u20[1];
        p[4] = x2d[0]; p[5] = x2d[1];
        p[6] = x3d[0]; p[7] = x3d[1];     
    } else if (my_id == 2){
        p[0] = x30[0]; p[1] = x30[1];
        p[2] = u30[0]; p[3] = u30[1];
        p[4] = x3d[0]; p[5] = x3d[1];
        p[6] = x4d[0]; p[7] = x4d[1]; 
    } else if (my_id == 3){
        p[0] = x40[0]; p[1] = x40[1];
        p[2] = u40[0]; p[3] = u40[1];
        p[4] = x4d[0]; p[5] = x4d[1];
    }

    std::chrono::milliseconds Ts = std::chrono::milliseconds(200);
    dmpcAgent my_dmpc_agent(folderName,my_id,1.0,&lcm,meta_,5,Ts,4,p);

    std::cout << "starting initQP" << std::endl;
    my_dmpc_agent.build_QP();
    std::cout << "build_QP finished" << std::endl;

    my_dmpc_agent.init_DMPC();
    my_dmpc_agent.run_DMPC();

    finished = true;
    lcm_handling_thread.join();    

    return 0;
}
