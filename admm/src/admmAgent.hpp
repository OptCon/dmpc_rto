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

#ifndef ADMMAGENT_H
#define ADMMAGENT_H


#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <numeric> //accumulate
#include <vector>
#include <math.h> //pow, max
#include <iostream>
#include <fstream>
#include <thread>

#include "sProb.hpp"
#include <lcm/lcm-cpp.hpp>
#include <qpOASES.hpp>
#include <cfloat> //DBL_MAX
#include "vcpy.hpp"
#include "handler_vector_idx_t.hpp"
#include "doptTimer.hpp"

using namespace Eigen;

class admmAgent{
    public:
        //problem data
        int Nagents;
        int my_id;

        int nx;
        int ng;
        int nh;
        int Ncons;
        sProb m_sprob;

        Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> H_bar;
        Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A;
        VectorXd g;
        VectorXd g_bar;
        VectorXd lbA;
        VectorXd ubA;
        VectorXd ub_vec;
        VectorXd lb_vec;

        Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> z_log;
        Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> zbar_log;
        Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> gam_log;

        qpOASES::Options myOptions;
        qpOASES::QProblem loc_prob;
        int nWSR;
        double rho;
        VectorXd lb;
        VectorXd ub; 

        //problem metadata
        Vector<bool, Eigen::Dynamic> isOriginal;    //nx x 1
        Vector<bool, Eigen::Dynamic> isCopy;        //nx x 1
        Vector<int, Eigen::Dynamic> numCopies;      //nx x 1 (how many copies there are of each original variable)    
        int N_og;                                   //Number of original variables this agent owns
        std::map<int,int> og_idx_to_idx;            //original variable index -> local variable index
        std::map<int,int> idx_to_og_idx;            //local variable index -> original variable index

        //ADMM iterates
        VectorXd z;                                 //nx x 1
        VectorXd z_bar;                             //nx x 1
        VectorXd gam;                               //nx x 1

        vCpy* v_in;     //copy from in-neighbors       //N_in_neighbors x 1
        vCpy* v_out;    //copy for out-neighbors       //N_out_neighbors x 1
        std::vector<double>* XV;                       //N_og x n_out_neighbors(i)

        std::vector<int> in_neighbors;              //N_in_neighbors x 1
        std::vector<int> out_neighbors;             //N_out_neighbors x 1
        int N_in_neighbors;
        int N_out_neighbors;

        admmAgent(const sProb& sprob_, const int& my_id_, const double& rho_, lcm::LCM* lcm);  
        ~admmAgent();                                           
        
        void init_coupling();   //extract coupling metadata from coupling matrices        
        
        int solve(unsigned int maxiter_, Eigen::VectorXd& z, Eigen::VectorXd& zbar_, Eigen::VectorXd& nu_, Eigen::VectorXd& mu_, Eigen::VectorXd& gam_, bool createLog);

        //before averaging
        void update_v_in();     //place z into v_in
        void send_vin_receive_vout(int64_t seq_number_);

        //after averaging
        void update_v_out();
        void send_vout_receive_vin(int64_t seq_number_);

        void init_comms();

        //communication
        lcm::LCM* comms; 
        itmessage::vector_idx_t*  v_outMsg;
        itmessage::vector_idx_t*  v_outMsg_old;
        VectorXi*  v_outMsg_idx_firstReceived;
        std::string* v_outTopic;                     //send all originals to all out-neighbors. could be replaced by something more efficient (goesta)
        std::string* v_inTopic; 

        itmessage::vector_idx_t *v_inMsg;
        itmessage::vector_idx_t *v_inMsg_old;
        VectorXi* v_inMsg_idx_firstReceived;

        
        lcm::Subscription** v_outSubscriptions;
        handler_vector_idx_t** v_outHandlers;

        lcm::Subscription** v_inSubscriptions;
        handler_vector_idx_t** v_inHandlers;


        //timing
        doptTimer loc_timer;        //time for solving the local QP
        doptTimer iter_timer;       //time for running one iteration
        doptTimer z_comm_timer;     //time for communicating z in one iteration
        doptTimer zbar_comm_timer;  //time for communicating zbar in one iteration
        doptTimer send_vin_timer;
        doptTimer receive_vout_timer;

        int update_g_beq();

        void print_time_measurements();
        void clear_time_measurements();
        void reserve_time_measurements(unsigned int new_cap);

};


#endif