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

#ifndef DMPCAGENT_H
#define DMPCAGENT_H
#include "dsqpAgent.hpp"
#include "Eigen/Dense"
#include "doptTimer.hpp"
#include <vector>
#include <thread>
#include <chrono>
#include <lcm/lcm-cpp.hpp>
#include "handler_stateVector.hpp"
#include <casadi/casadi.hpp>
struct mpcMetaData{
    int nx; //state dimension
    int nu; //input dimension
    int nxd; //dimension of setpoint
    int idx_eqx0; //equality constraint for initial condition
    int idx_equ0; //equality constraint for u0
    int idx_u0;
    int idx_u1;
    int idx_u2;
    std::vector<int> shift_seq_start;
    std::vector<int> shift_seq_length;
    std::vector<int> shift;
}; 


class dmpcAgent{

    public:
        VectorXd m_x_meas;
        VectorXd m_u0;      //current control input
        VectorXd m_u1;      //next control input
        VectorXd m_u2;      //backup
        VectorXd m_sol;
        VectorXd m_z;
        VectorXd m_zbar;
        VectorXd m_gam;
        VectorXd m_nu;
        VectorXd m_mu;

        mpcMetaData m_meta;
        VectorXd m_p; //parameters for OCP ([x0; u0; xd])
        VectorXd m_xd;
        int m_numSteps;

        int m_sqp_maxiter;
        int m_admm_maxiter;
        double m_rho;

        int init_DMPC();
        int run_DMPC();

        std::chrono::milliseconds m_Ts; //ms
        int m_my_id;
        std::string m_folderName; //folder with sProb
        dmpcAgent(const std::string& folderName, int my_id_, double rho_, lcm::LCM* lcm, const mpcMetaData& meta_, int sqp_maxiter_, int admm_maxiter_,std::chrono::milliseconds Ts_, unsigned int Nagents, const Eigen::VectorXd& p);
        ~dmpcAgent();

        // communication
        //! comms_init
        /* Initializes communication for the DMPC controller. Do not forget to call comms_destroy() at the end of the programme if comms_init has been called. 
        * \param[in] lcm pointer to the LCM instance to be used for communication
        * \param[in] TsState sampling time of the state/pose measurements (e.g., frame time of the tracking system employed);  TsC/TsState should be an even number, otherwise problems may occur
        * \param[out] state_step integer - how many pose updates received per control sampling interval
        */
        int comms_init(lcm::LCM* lcm, double TsState);
        int update_setpoint_in_ocp();
        void comms_destroy();
        int idx_agent = -1;
        int Nagents = -1;
        lcm::LCM* m_lcm;
        double TsC; // control sampling time in seconds
        lcm::LCM* comms;
        itmessage::vector_t controlInputMsg;
        std::string controlInputTopic;
        lcm::Subscription** state_vector_subscriptions;

        void print_time_measurements();
        void print_ocp_solution_log(const Eigen::MatrixXd& z_log);

    private: 
        dsqpAgent* m_dsqp;
        doptTimer m_mpc_timer;
        doptTimer m_receivex_timer;
        doptTimer m_sendu_timer;
        doptTimer m_dsqp_timer;
        doptTimer m_sleep_timer;
        

        handler_stateVector* stateVectorHandler;
        void set_xmeas_in_ocp();
        void set_u0_in_ocp();
        void get_u_from_sol();

        Eigen::MatrixXd m_x_log;
        Eigen::MatrixXd m_u_log;
        Eigen::MatrixXd m_xd_log;

        // communication
        //! comms_init
        /* Initializes communication for the DMPC controller. Do not forget to call comms_destroy() at the end of the programme if comms_init has been called. 
        * \param[in] t current time step
        * \param[in] state_step how many pose updates received per control sampling interval (TsState/TsC - if the latter results in an integer)
        * \param[in] experimental_mode if true, only wait waiting_timeout seconds for new pose measurements until continuing with (partly) old values
        * \param[in] waiting_timeout see experimental_mode; value in seconds
        */
        void receive_xmeas(int t, int state_step, bool experimental_mode, double waiting_timeout); //receive from sensor (hebel)
        
        // communication
        //! send_u_to_robot
        /* Sends the robot's control input (channel /robot[idx_agent]/u).
        * \param[in] u three-dimensional control input vector (first two dimensions translational velocity, third dimension angular velocity)
        */
        void send_u_to_robot(const Eigen::Ref<const Vector3d>& u); // (hebel)
        void warm_start_solver(); 


        Eigen::MatrixXd casadi2Eigen ( const casadi::DM& A );
        Eigen::VectorXd casadi2EigenVector ( const casadi::DM& A );
        casadi::DM Eigen2casadi( const Eigen::VectorXd& in);


};

#endif