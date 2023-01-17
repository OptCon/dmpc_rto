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


#include "dmpcAgent.hpp"


dmpcAgent::dmpcAgent(const std::string& folderName_, int my_id_, double rho_, lcm::LCM* lcm, const mpcMetaData& meta_, int sqp_maxiter_, int admm_maxiter_,std::chrono::milliseconds Ts_, unsigned int Nagents_, const Eigen::VectorXd& p){

    m_folderName = folderName_;
    m_meta = meta_;
    m_sqp_maxiter = sqp_maxiter_;
    m_admm_maxiter = admm_maxiter_;
    m_Ts = Ts_;
    m_my_id = my_id_;

    idx_agent = my_id_;
    Nagents = Nagents_;
    TsC = ((double) Ts_.count())*0.001; // ms to s
    m_u0 = Eigen::VectorXd::Zero(meta_.nu);
    m_u1 = Eigen::VectorXd::Zero(meta_.nu);
    m_u2 = Eigen::VectorXd::Zero(meta_.nu);
    m_x_meas = Eigen::VectorXd::Zero(meta_.nx);

    m_lcm = lcm;
    m_p = p;
    m_rho = rho_;

    m_xd = m_p.segment(m_meta.nx+m_meta.nu,m_meta.nxd);

    m_dsqp = new dsqpAgent(m_folderName, m_my_id, Nagents, m_rho, m_sqp_maxiter, m_admm_maxiter, m_lcm);
    m_dsqp->m_p = m_p;
}


dmpcAgent::~dmpcAgent(){
    if (m_dsqp != nullptr){
        delete m_dsqp;
    }
}

int dmpcAgent::comms_init(lcm::LCM* lcm, double TsState) {
    comms = lcm;

    controlInputMsg.seq_number = 0;
    controlInputMsg.length = 3;
    controlInputMsg.value = {0.0, 0.0, 0.0};
    controlInputMsg.id_sender = idx_agent + 1;
    controlInputTopic = "/robot";
    controlInputTopic.append(std::to_string(idx_agent + 1)); // robots are numbered beginning at 1
    controlInputTopic.append("/u");
    // usage: comms->publish(controlInputTopic, &controlInputMsg);
    const int state_step = std::round(TsC/TsState); // TsC/TsState should be an even number, otherwise problems may occur
    std::cout << "state_step = " << state_step << std::endl;
    int num_state_aggregators = Nagents;
	stateVectorHandler = new handler_stateVector(-state_step+1, idx_agent+1, Nagents, state_step); // agentNumber+1 is Matlab numbering of agents
	state_vector_subscriptions = new lcm::Subscription*[num_state_aggregators];
	for (int idx_state_aggregator=0; idx_state_aggregator<num_state_aggregators; idx_state_aggregator++) {
        // deal with the robots' states
        std::string cur_topic = "/robot";
        cur_topic.append(std::to_string(idx_state_aggregator + 1)); // robots are numbered beginning at 1
        cur_topic.append("/euler");
        state_vector_subscriptions[idx_state_aggregator] = lcm->subscribe(cur_topic, &handler_stateVector::handleMessage, stateVectorHandler);

        state_vector_subscriptions[idx_state_aggregator]->setQueueCapacity(1);
        std::cout << "subscribed to " << cur_topic << std::endl;
    }

    return state_step;
}

void dmpcAgent::comms_destroy() {
    delete[] state_vector_subscriptions;
    delete stateVectorHandler;
}

void lcm_handler_threadfunction(lcm::LCM* lcm, int64_t* t, int64_t* numSteps, bool* finished) {
    while (((*t) <= (*numSteps)) && (!(*finished))) {
        lcm->handleTimeout(1);
    }
}

int dmpcAgent::run_DMPC(){
    bool experiment_mode = false; // do not wait forever for state measurement

    double TsState = 0.01; // sampling time of the state/pose measurements (e.g., frame time of the tracking system employed);  TsC/TsState should be an even number, otherwise problems may occur
    int state_step = comms_init(m_lcm, TsState) ;
    Vector3d u = Vector3d::Zero(3);


    int64_t t = 1;
    m_numSteps = 50; // rectangle: 500
    int64_t numSteps = m_numSteps;
    bool finished = false;
    std::thread lcm_handling_thread(lcm_handler_threadfunction,m_lcm,&t,&numSteps,&finished);

    m_dsqp->clear_time_measurements();

    m_mpc_timer.reserve(m_numSteps);
    m_receivex_timer.reserve(m_numSteps);
    m_dsqp_timer.reserve(m_numSteps);
    m_sleep_timer.reserve(m_numSteps);
    m_sendu_timer.reserve(m_numSteps);
    m_dsqp->reserve_time_measurements(m_numSteps*m_sqp_maxiter);
    m_dsqp->m_admm->reserve_time_measurements(m_numSteps*m_sqp_maxiter*m_admm_maxiter);

    m_x_meas[0] = m_p[0]; //in simulation. uncomment in experiments
    m_x_meas[1] = m_p[1]; //in simulation. uncomment in experiments

    MatrixXd z_log = -MatrixXd::Ones(m_dsqp->nx+1,numSteps);
    m_x_log = -MatrixXd::Ones(m_numSteps,m_meta.nx);
    m_u_log = -MatrixXd::Ones(m_numSteps,m_meta.nu);
    m_xd_log = -MatrixXd::Ones(m_numSteps,m_meta.nxd);

    Eigen::VectorXd x_next = Eigen::VectorXd::Zero(2);

    int dsqp_flag = 0;

    std::chrono::milliseconds send_u_time_ = std::chrono::milliseconds(20); //ms
    int iter0 = t;
    auto t0 = steady_clock::now();

    m_dsqp->z = VectorXd::Zero(m_dsqp->nx);
    m_dsqp->gam = VectorXd::Zero(m_dsqp->nx);
    m_z = VectorXd::Zero(m_dsqp->nx);

    std::chrono::milliseconds awake_at_ms_since_t0;
    auto awake_time = t0;

    while (t <= m_numSteps) {

        m_mpc_timer.tic();
       
        m_receivex_timer.tic();
        // receive_xmeas(t, state_step, experiment_mode, 0.05); //run in experiments
        m_receivex_timer.toc();

        std::cout << "t = " << t << " , x = " << m_x_meas.transpose() << std::endl;
        
        set_xmeas_in_ocp();
        set_u0_in_ocp();
        

        if (m_my_id == 0) {
            m_xd(0) = 1.0;
            m_xd(1) = 0.9;
            m_xd(2) = -1.7; 
            m_xd(3) = -0.6;
        }
        if (m_my_id == 1) {
            m_xd(0) = -1.7;
            m_xd(1) = -0.6;
            m_xd(2) = 0.75;
            m_xd(3) = 0.0;
        }
        if (m_my_id == 2) {
            m_xd(0) = 0.75;
            m_xd(1) = 0.0;
            m_xd(2) = -1.625;
            m_xd(3) = -1.2;
        }
        if (m_my_id == 3) {
            m_xd(0) = -1.625;
            m_xd(1) = -1.2;
        }

        double step1 = 20.0/TsC;
        if (t>=step1) {
            if (m_my_id == 0) {
                m_xd(0) = 0.875;
                m_xd(1) = -0.4;
                m_xd(2) = -0.5;
                m_xd(3) = 0.0;
            }
            if (m_my_id == 1) {
                m_xd(0) = -0.5;
                m_xd(1) = 0.0;
                m_xd(2) = -0.5;
                m_xd(3) = 0.0;
            }
            if (m_my_id == 2) {
                m_xd(0) = -0.5;
                m_xd(1) = 0.0;
                m_xd(2) = -0.5;
                m_xd(3) = 0.0;
            }
            if (m_my_id == 3) {
                m_xd(0) = -0.5;
                m_xd(1) = 0.0;
            }
        }

        
        // rectangle
        // if(m_my_id == 0){
        //     m_xd(0) = -0.6;
        //     m_xd(1) = 0.0;
        //     double step1 = 15.0/TsC;
        //     double divisor = 4.0/3.0;
        //     double v_des = 0.2/divisor;
        //     if (t>=step1) {
        //         m_xd(0) = -0.6;
        //         m_xd(1) = 0.0+v_des*(t-step1)*TsC;
        //     }
        //     double step2 = 5.0*divisor/TsC+step1;
        //     if (t>=step2) {
        //         m_xd(0) = -0.6+v_des*(t-step2)*TsC;
        //         m_xd(1) = 1.0;
        //     }
        //     double step3 = 11.0*divisor/TsC+step2;
        //     if (t>=step3) {
        //         m_xd(0) = 1.6;
        //         m_xd(1) = 1.0-v_des*(t-step3)*TsC;
        //     }
        //     double step4 = 10.0*divisor/TsC+step3;
        //     if (t>=step4) {
        //         m_xd(0) = 1.6-v_des*(t-step4)*TsC;
        //         m_xd(1) = -1.0;
        //     }
        //     double step5 = 11.0*divisor/TsC+step4;
        //     if (t>=step5) {
        //         m_xd(0) = -0.6;
        //         m_xd(1) = -1.0+v_des*(t-step5)*TsC;
        //     }
        //     double step6 = 5.0*divisor/TsC+step5;
        //     if (t>=step6) {
        //         m_xd(0) = -0.6;
        //         m_xd(1) = 0.0;
        //     }
        //     update_setpoint_in_ocp();
        // }

        update_setpoint_in_ocp();

        m_x_log.block(t-1,0,1,m_meta.nx) = m_x_meas.transpose();
        m_u_log.block(t-1,0,1,m_meta.nu) = m_u0.transpose();
        m_xd_log.block(t-1,0,1,m_meta.nxd) = m_xd.transpose();

        m_dsqp_timer.tic();       
        dsqp_flag = m_dsqp->solve(false);        
        m_dsqp_timer.toc();
        //compute new t0 if last mpc step was too slow (e.g. because receivex was too slow)
        if (t>1){
            if (m_mpc_timer.m_d_us > m_Ts - send_u_time_){ //last mpc step was too slow
                t0 = steady_clock::now();
                iter0 = t;
            }
        }

        m_zbar = m_dsqp->z;
        m_nu = m_dsqp->nu;
        m_mu = m_dsqp->mu;
        m_gam = m_dsqp ->gam;
        m_sol = m_zbar; //plain ADMM: z; dSQP: zbar
        z_log(0,t-1) = t;
        z_log.block(1,t-1,m_dsqp->nx,1) = m_sol;
                
        get_u_from_sol();
       
        u.head(m_meta.nu) = m_u1.head(m_meta.nu); 
        warm_start_solver();
             
        m_mpc_timer.toc();
       
        m_sleep_timer.tic();
        awake_at_ms_since_t0 = (t-iter0)*m_Ts - send_u_time_;
        awake_time = t0 + awake_at_ms_since_t0;
        std::this_thread::sleep_until(awake_time);
        m_sleep_timer.toc();

        m_sendu_timer.tic();
        // send_u_to_robot(u); //run in experiments
        x_next = m_x_meas + m_u0* (double) m_Ts.count()/1000; //only in cpp simulation
        m_u0 = m_u1;
        m_x_meas = x_next; //only in cpp simulation
        t++;
        m_sendu_timer.toc();

        awake_at_ms_since_t0 = (t-iter0-1)*m_Ts; //-1 because the increment is in line 339
        awake_time = t0 + awake_at_ms_since_t0;
        std::this_thread::sleep_until(awake_time);

    }

    u.setZero();
    for (int i=0; i<20; i++) {
        if (!experiment_mode) {
            controlInputMsg.value[0] = u(0);
            controlInputMsg.value[1] = u(1);
            controlInputMsg.value[2] = u(2);
            comms->publish(controlInputTopic, &controlInputMsg);
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        } else {
            send_u_to_robot(u);
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }

    print_time_measurements();
    print_ocp_solution_log(z_log);


    std::cout << "finished" << std::endl;
    lcm_handling_thread.join();
    comms_destroy();
    return 0;
}

void dmpcAgent::set_xmeas_in_ocp(){
    m_p.segment(0,m_meta.nx) = m_x_meas;
    m_dsqp->m_p = m_p;
    return;
}

void dmpcAgent::set_u0_in_ocp(){
    m_p.segment(m_meta.nx,m_meta.nu) = m_u0;
    m_dsqp->m_p = m_p;
    return;
}

void dmpcAgent::get_u_from_sol(){
    m_u0 = m_sol.block(m_meta.idx_u0,0,m_meta.nu,1);
    m_u1 = m_sol.block(m_meta.idx_u1,0,m_meta.nu,1);
    m_u2 = m_sol.block(m_meta.idx_u2,0,m_meta.nu,1);
}

void dmpcAgent::warm_start_solver(){

    int start = 0;
    int length = 0;
    int shift = 0;

    unsigned int N_seq = m_meta.shift.size();
    for (unsigned int i = 0; i < N_seq; i++){
        start = m_meta.shift_seq_start[i];
        length = m_meta.shift_seq_length[i];
        shift = m_meta.shift[i];
        m_zbar.block(start,0,length-shift,1) = m_zbar.block(start+shift,0,length-shift,1);
        m_gam.block(start,0,length-shift,1) = m_gam.block(start+shift,0,length-shift,1);
    }

    m_dsqp->z = m_zbar;
    m_dsqp->gam = m_gam;
}

void dmpcAgent::receive_xmeas(int t, int state_step, bool experimental_mode, double waiting_timeout) {
    //std::cout << "receive_xmeas, t = " << t << ", stateVectorHandler->StateSeq = " << stateVectorHandler->StateSeq << ", state_step * (t - 1) = " << state_step * (t - 1)<< std::endl;
    bool waitingTimeOver= false;
    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
    double elapsedTime = ((double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()) / 1000000.0;
    bool receivedAllStates = true;
    while ((stateVectorHandler->StateSeq < state_step * (t - 1)) && ((!experimental_mode) || (!waitingTimeOver))) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        if (!experimental_mode) {
            comms->publish(controlInputTopic, &controlInputMsg);
        }
        if (t > 1) {
            // later: change to t > 1 to, in the first time step, wait "indefinitely"
            endTime = std::chrono::steady_clock::now();
            elapsedTime = ((double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()) / 1000000.0;
            if (elapsedTime > waiting_timeout) {
                receivedAllStates = false;
                waitingTimeOver = true;
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        // std::cout << "stateVectorHandler->StateSeq = " << stateVectorHandler->StateSeq << ">=!" << state_step * (t - 1) << std::endl;
    }

    if (experimental_mode) {
        stateVectorHandler->processMessages(false);
    } else {
        stateVectorHandler->processMessages(true);
    }

    // stateVectorHandler->robot_poses now contains a robot pose in each column (planar position plus orientation in radian), ordered by the robot id
    m_x_meas[0] = stateVectorHandler->robot_poses(0,m_my_id); //x-position
    m_x_meas[1] = stateVectorHandler->robot_poses(1,m_my_id); //y-position
    std::cout << "m_x_meas = " << m_x_meas[0] << ", " << m_x_meas[1] << std::endl; 
    return;
}

void dmpcAgent::send_u_to_robot(const Eigen::Ref<const Vector3d>& u){
    controlInputMsg.value[0] = u(0);
    controlInputMsg.value[1] = u(1);
    controlInputMsg.value[2] = u(2);
    controlInputMsg.seq_number = controlInputMsg.seq_number + 1;
    comms->publish(controlInputTopic, &controlInputMsg);
    return;
}

int dmpcAgent::init_DMPC(){
    
    int nx = m_dsqp->m_qp.A[m_my_id].cols();
    if (m_z.size() != nx){
        m_z = VectorXd::Zero(nx);
    }

    int flag = m_dsqp->init(&m_z,nullptr,nullptr);
    // int sqp_maxiter_ = 10;
    // int admm_maxiter_ = 10;
    // flag = m_dsqp->solve(sqp_maxiter_,admm_maxiter_);
    m_z = VectorXd::Zero(nx);

    return 0;
}

int dmpcAgent::update_setpoint_in_ocp(){

    m_p.segment(0,m_meta.nx) = m_x_meas;
    m_p.segment(m_meta.nx,m_meta.nu) = m_u0;
    m_p.segment(m_meta.nx+m_meta.nu,m_meta.nxd) = m_xd;

    return 0;


}

void dmpcAgent::print_time_measurements(){
    std::time_t t = std::time(0);   // get time now
    std::tm* now = std::localtime(&t);
 
    std::ostringstream fileName;
    fileName << "dmpc_time_measurement" << "_agent" << m_my_id << "_" << (now->tm_year + 1900) << '_' << (now->tm_mon + 1) << '_' <<  now->tm_mday << "_" << now->tm_hour << "_" << now->tm_min << "_" << now->tm_sec <<".csv";
    std::ofstream file(fileName.str());
    if (file.is_open())
    {
    file << "mpc step, mpc step time (us), x0(1), x0(2), u0(1), u0(2), xd(1), xd(2), receivex time (us), sqp time (us), sleep time (us), sendu time (us), sqp iter, sqp iter time (us), buildQP time(us), reg time (us), updateQP time(us), admm time (us), admm iter, admm iter time (us), loc_qp time(us), zcomm time (us), zbarcomm time(us), sendvin time (us), receivevout time (us)\n";
    file.close();
    }

    int N_rows = m_dsqp->iter_timer.m_log.size();
    int k = 1; //MPC step
    unsigned int row = 0;

    int dsqp_iter = 0;
    int admm_iter = 0;
    for (int k = 1; k <= m_numSteps; k++){
        for (unsigned int i = 0; i < m_sqp_maxiter; i++){
            for (unsigned int j = 0; j < m_admm_maxiter; j++){
                file.open(fileName.str(),std::ios_base::app);
                if (file.is_open())
                {
                    file << k << "," << m_mpc_timer.m_log[k-1] << "," << m_x_log(k-1,0) << "," << m_x_log(k-1,1) << "," << m_u_log(k-1,0) << "," << m_u_log(k-1,1) << "," << m_xd_log(k-1,0) << "," << m_xd_log(k-1,1) << "," << m_receivex_timer.m_log[k-1] << "," << m_dsqp_timer.m_log[k-1] << "," << m_sleep_timer.m_log[k-1] << "," << m_sendu_timer.m_log[k-1] << "," << i+1 << "," << m_dsqp->iter_timer.m_log[dsqp_iter] << "," << m_dsqp->buildQP_timer.m_log[dsqp_iter] << "," << m_dsqp->reg_timer.m_log[dsqp_iter] << "," << m_dsqp->updateQP_timer.m_log[dsqp_iter] << "," << m_dsqp->admm_timer.m_log[dsqp_iter] << "," << j+1 << "," << m_dsqp->m_admm->iter_timer.m_log[admm_iter] << "," << m_dsqp->m_admm->loc_timer.m_log[admm_iter] << "," << m_dsqp->m_admm->z_comm_timer.m_log[admm_iter] << "," << m_dsqp->m_admm->zbar_comm_timer.m_log[admm_iter] << "," << m_dsqp->m_admm->send_vin_timer.m_log[admm_iter] << "," << m_dsqp->m_admm->receive_vout_timer.m_log[admm_iter] <<"\n";        
                }
                file.close();
                admm_iter+=1;
            }
        dsqp_iter+=1;
        }
    }

    
    return;
}

void dmpcAgent::print_ocp_solution_log(const MatrixXd& z_log_){
    std::time_t t = std::time(0);
    std::tm* now = std::localtime(&t);
 
    std::ostringstream fileName;
    fileName << "ocp_solution_log" << "_agent" << m_my_id << "_" << (now->tm_year + 1900) << '_' << (now->tm_mon + 1) << '_' <<  now->tm_mday << "_" << now->tm_hour << "_" << now->tm_min << "_" << now->tm_sec <<".csv";
    std::ofstream file(fileName.str());
    if (file.is_open())
    {
        file << z_log_;
        file.close();
    }
}

// Dense to Dense
Eigen::MatrixXd dmpcAgent::casadi2Eigen ( const casadi::DM& A ){
    // This method is based on code by Petr Listov, see https://groups.google.com/g/casadi-users/c/npPcKItdLN8
    casadi::Sparsity SpA = A.get_sparsity();
    std::vector<long long int> output_row, output_col;
    SpA.get_triplet(output_row, output_col);
    std::vector<double> values = A.get_nonzeros();
    using T = Eigen::Triplet<double>;
    std::vector<T> TripletList;
    TripletList.resize(values.size());
    for(int k = 0; k < values.size(); ++k){
        if (values[k] == casadi::inf){
            values[k] = DBL_MAX;
        }
        else if (values[k] == -casadi::inf){
            values[k] = -DBL_MAX;
        }
        TripletList[k] = T(output_row[k], output_col[k], values[k]);
    }
    Eigen::SparseMatrix<double> SpMatrx(A.size1(), A.size2());
    SpMatrx.setFromTriplets(TripletList.begin(), TripletList.end());

    Eigen::MatrixXd DeMatrx = SpMatrx.toDense();
    return DeMatrx;
}

// Dense to Dense
Eigen::VectorXd dmpcAgent::casadi2EigenVector ( const casadi::DM& A ){
    // This method is based on code by Petr Listov, see https://groups.google.com/g/casadi-users/c/npPcKItdLN8
    casadi::Sparsity SpA = A.get_sparsity();
    std::vector<long long int> output_row, output_col;
    SpA.get_triplet(output_row, output_col);
    std::vector<double> values = A.get_nonzeros();
    using T = Eigen::Triplet<double>;
    std::vector<T> TripletList;
    TripletList.resize(values.size());
    for(int k = 0; k < values.size(); ++k){
        if (values[k] == casadi::inf){
            values[k] = DBL_MAX;
        }
        else if (values[k] == -casadi::inf){
            values[k] = -DBL_MAX;
        }
        TripletList[k] = T(output_row[k], output_col[k], values[k]);
    }
    
    Eigen::VectorXd DeVec;
    if (A.size2() == 1){ // column-vector
        DeVec = Eigen::VectorXd::Zero(A.size1());
        
        for (int k = 0; k < values.size(); k++){
            DeVec[output_row[k]] = values[k];
        }
    } else { //row-vector (will be transposed)
        DeVec = Eigen::VectorXd::Zero(A.size2());
        for (int k = 0; k < values.size(); k++){
            DeVec[output_col[k]] = values[k];
        }
    }

    for (int k = 0; k < DeVec.size(); k++){
         if (DeVec[k] == casadi::inf){
            DeVec[k] = DBL_MAX;
        }
    }

    return DeVec;
}

casadi::DM dmpcAgent::Eigen2casadi( const Eigen::VectorXd& in){
    int length = in.size();
    casadi::DM out = casadi::DM::zeros(length,1);
    std::memcpy(out.ptr(), in.data(), sizeof(double)*length*1);
    return out;
}
