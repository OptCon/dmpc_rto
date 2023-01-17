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

#ifndef DSQP_AGENT_H
#define DSQP_AGENT_H

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
#include "admmAgent.hpp"
#include <casadi/casadi.hpp>



class dsqpAgent{
    public:
        //problem data
        unsigned int Nagents;
        int my_id;

        int nx;
        int ng;
        int nh;
        int Ncons;
        sProb m_qp;
        admmAgent* m_admm;
        std::string m_folderName;
        double m_rho;

        Eigen::VectorXd m_p;
        int m_sqp_maxiter;
        int m_admm_maxiter;

        //SQP iterates
        VectorXd z;                                 //nx x 1; "zbar" in ADMM
        VectorXd nu;                                //ng x 1
        VectorXd mu;                                //nh x 1
        VectorXd gam;                               //nx x 1        

        dsqpAgent(const std::string& folderName, int my_id_, unsigned int Nagents,  double rho_, unsigned int sqp_maxiter_, unsigned int admm_maxiter_,lcm::LCM* lcm);  
        ~dsqpAgent();                                           
        
        
        int solve(bool createLog);       //solve NLP
        int solve(int sqp_maxiter_, int admm_maxiter_, bool createLog); //solve NLP

        int buildQP(Eigen::VectorXd* z, Eigen::VectorXd* nu, Eigen::VectorXd* mu, bool GN, bool eval_HessF);

        int init();
        int init(Eigen::VectorXd* z, Eigen::VectorXd* nu, Eigen::VectorXd* mu);

        lcm::LCM* m_lcm;

        Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> z_log;

        doptTimer buildQP_timer;    //time for building the local QP
        doptTimer reg_timer;        //time for regularizing the QP
        doptTimer updateQP_timer;   //time for updating the QP in the admmAgent
        doptTimer iter_timer;       //time for running one sqp iteration
        doptTimer admm_timer;       //time for running ADMM in one sqp iteration

        void clear_time_measurements();
        void reserve_time_measurements(unsigned int new_cap);
    
    private:
        Eigen::MatrixXd casadi2Eigen ( const casadi::DM& A );
        Eigen::VectorXd casadi2EigenVector ( const casadi::DM& A );
        casadi::DM Eigen2casadi( const Eigen::VectorXd& in);

        VectorXd z_;                                //nx x 1; "z" in ADMM
        VectorXd delz;                              //nx x 1; 

        casadi::DM z_cas;
        casadi::DM nu_cas;
        casadi::DM mu_cas;
        casadi::DM p_cas;
    
        // Use CasADi's "external" to load the compiled function
        casadi::Function f; // = casadi::external("gradFun1","sProb_chain/locFuns.so");
        std::vector<casadi::DM> arg; // = {z};
        std::vector<casadi::DM> res; // = f(arg);
        std::string functionLibrary;
        std::string str;

        VectorXd ub_file; //upper bounds as specified in csv file
        VectorXd lb_file; //lower bounds as specified in csv file

};

#endif 