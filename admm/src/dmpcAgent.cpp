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

dmpcAgent::dmpcAgent(const std::string& folderName_, int my_id_, double rho_, lcm::LCM* lcm, const mpcMetaData& meta_, int maxiter_, std::chrono::milliseconds Ts_, unsigned int Nagents_, const Eigen::VectorXd& p){

    m_folderName = folderName_;
    m_admm = nullptr;
    m_meta = meta_;
    m_maxiter = maxiter_;
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
}


dmpcAgent::~dmpcAgent(){
    if (m_admm != nullptr){
        delete m_admm;
    }
}

int dmpcAgent::build_QP(){

    sProb qp = sProb(Nagents);
    qp.read_AA(m_folderName,Nagents);
    qp.read_ublb(m_folderName,Nagents);

    std::string functionLibrary = m_folderName + "/locFuns.so";

    casadi::DM z_cas;
    casadi::DM p_cas;

    // Use CasADi's "external" to load the compiled function
    casadi::Function f; // = casadi::external("gradFun1","sProb_chain/locFuns.so");
    std::string str;

    std::vector<casadi::DM> arg; // = {z};
    std::vector<casadi::DM> res; // = f(arg);

    int nx = qp.A[m_my_id].cols();


    VectorXd z = Eigen::VectorXd::Zero(nx);
  
    z_cas = Eigen2casadi(z);
    p_cas = Eigen2casadi(m_p);

    //H
    str = "HessFfun" + std::to_string(m_my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    arg = {z_cas,p_cas};
    res = f(arg);
    casadi::DM H = res[0];
    MatrixXd Heig = casadi2Eigen(H);
    qp.H[m_my_id] = Heig;

    //g
    str = "gradFun" + std::to_string(m_my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    res = f(arg);
    qp.g[m_my_id] = Eigen::VectorXd::Zero(nx);
    qp.g[m_my_id] = casadi2EigenVector(res[0]);

    //Aeq
    str = "JGfun" + std::to_string(m_my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    res = f(arg);
    casadi::DM Aeq = res[0];
    qp.Aeq[m_my_id] = casadi2Eigen(Aeq);

    //Aineq
    str = "JHfun" + std::to_string(m_my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    res = f(arg);
    casadi::DM Aineq = res[0];
    qp.Aineq[m_my_id] = casadi2Eigen(Aineq);

    //beq
    str = "eqfun" + std::to_string(m_my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    res = f(arg);
    int ng = res[0].size1();
    int cols = res[0].size2();
    qp.beq[m_my_id] = Eigen::VectorXd::Zero(ng);        
    std::memcpy(qp.beq[m_my_id].data(), res.at(0).ptr(), sizeof(double)*ng*cols);
    qp.beq[m_my_id] = -qp.beq[m_my_id];

    //bineq
    str = "ineqfun" + std::to_string(m_my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    res = f(arg);
    int nh = res[0].size1();
    cols = res[0].size2();
    qp.bineq[m_my_id] = Eigen::VectorXd::Zero(nh);        
    std::memcpy(qp.bineq[m_my_id].data(), res.at(0).ptr(), sizeof(double)*nh*cols);
    qp.bineq[m_my_id] = -qp.bineq[m_my_id];

    VectorXd mlb = qp.lb[m_my_id];

    for (int k = 0; k < qp.lb[m_my_id].size(); k++){
        if (qp.lb[m_my_id][k] == -casadi::inf){
            qp.lb[m_my_id][k] = -pow(10,20); 
        } else if (qp.lb[m_my_id][k] == casadi::inf){
            qp.lb[m_my_id][k] = pow(10,20);
        }
    }

    for (int k = 0; k < qp.ub[m_my_id].size(); k++){
        if (qp.ub[m_my_id][k] == -casadi::inf){
            qp.ub[m_my_id][k] = -pow(10,20);
        } else if (qp.ub[m_my_id][k] == casadi::inf){
            qp.ub[m_my_id][k] = pow(10,20);
        }
    }

    m_admm = new admmAgent(qp,m_my_id,m_rho,m_lcm);

    return 0;
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

    m_admm->clear_time_measurements();
    clear_time_measurements();

    int64_t t = 1;
    int64_t numSteps = 50;
    bool finished = false;
    std::thread lcm_handling_thread(lcm_handler_threadfunction,m_lcm,&t,&numSteps,&finished);

    m_mpc_timer.reserve(numSteps);
    m_receivex_timer.reserve(numSteps);
    m_sleep_timer.reserve(numSteps);
    m_sendu_timer.reserve(numSteps);
    m_admm_timer.reserve(numSteps);
    m_admm->reserve_time_measurements(numSteps*m_maxiter);

    m_x_meas[0] = m_p[0];
    m_x_meas[1] = m_p[1];


    MatrixXd z_log = -MatrixXd::Ones(m_admm->nx+1,numSteps);
    m_x_log = -MatrixXd::Ones(numSteps,m_meta.nx);
    m_u_log = -MatrixXd::Ones(numSteps,m_meta.nu);
    m_xd_log = -MatrixXd::Ones(numSteps,m_meta.nxd);

    Eigen::VectorXd x_next = Eigen::VectorXd::Zero(2);

    int flag = 0;
    auto t0 = steady_clock::now();
    int iter0 = t;
    std::chrono::milliseconds send_u_time_ = std::chrono::milliseconds(20); //ms

    while (t <= numSteps) {

        m_mpc_timer.tic();
       
        m_receivex_timer.tic();
        // receive_xmeas(t, state_step, experiment_mode, 0.05); //uncomment in cpp simulations; include this line in experiments
        m_receivex_timer.toc();



        std::cout << "t = " << t << " , x = " << m_x_meas.transpose() << std::endl;
        
        set_xmeas_in_ocp();
        set_u0_in_ocp();

        if(m_my_id == 0){
            m_xd(0) = -0.6;
            m_xd(1) = 0.0;
            double step1 = 15.0/TsC;
            double divisor = 4.0/3.0;
            double v_des = 0.2/divisor;
            if (t>=step1) {
                m_xd(0) = -0.6;
                m_xd(1) = 0.0+v_des*(t-step1)*TsC;
            }
            double step2 = 5.0*divisor/TsC+step1;
            if (t>=step2) {
                m_xd(0) = -0.6+v_des*(t-step2)*TsC;
                m_xd(1) = 1.0;
            }
            double step3 = 11.0*divisor/TsC+step2;
            if (t>=step3) {
                m_xd(0) = 1.6;
                m_xd(1) = 1.0-v_des*(t-step3)*TsC;
            }
            double step4 = 10.0*divisor/TsC+step3;
            if (t>=step4) {
                m_xd(0) = 1.6-v_des*(t-step4)*TsC;
                m_xd(1) = -1.0;
            }
            double step5 = 11.0*divisor/TsC+step4;
            if (t>=step5) {
                m_xd(0) = -0.6;
                m_xd(1) = -1.0+v_des*(t-step5)*TsC;
            }
            double step6 = 5.0*divisor/TsC+step5;
            if (t>=step6) {
                m_xd(0) = -0.6;
                m_xd(1) = 0.0;
            }
            update_setpoint_in_ocp();
        }

        m_x_log.block(t-1,0,1,m_meta.nx) = m_x_meas.transpose();
        m_u_log.block(t-1,0,1,m_meta.nu) = m_u0.transpose();
        m_xd_log.block(t-1,0,1,m_meta.nxd) = m_xd.transpose();

        m_admm_timer.tic();
        flag = m_admm->solve(m_maxiter,m_z,m_zbar,m_nu,m_mu,m_gam,false);
        m_admm_timer.toc();
        //compute new t0 if last mpc step was too slow (e.g. because receivex was too slow)
        if (t>1){
            if (m_mpc_timer.m_d_us > m_Ts - send_u_time_){ //last mpc step was too slow
                t0 = steady_clock::now();
                iter0 = t;
            }
        }
        m_sol = m_z; //plain ADMM: z; dSQP: zbar
        z_log(0,t-1) = t;
        z_log.block(1,t-1,m_admm->nx,1) = m_sol;

        get_u_from_sol();

        
        u.head(m_meta.nu) = m_u1.head(m_meta.nu);

        warm_start_solver();        
        
        m_mpc_timer.toc();
       
        m_sleep_timer.tic();
        std::chrono::milliseconds awake_at_ms_since_t0 = (t-iter0)*m_Ts - send_u_time_;
        auto awake_time = t0 + awake_at_ms_since_t0;
        std::this_thread::sleep_until(awake_time);
        m_sleep_timer.toc();

        x_next = m_x_meas + m_u0* (double) m_Ts.count()/1000; //cpp simulation (uncomment this line in experiments)
        m_sendu_timer.tic();
        // send_u_to_robot(u);
        m_u0 = m_u1;
        t++;
        m_sendu_timer.toc();

        awake_at_ms_since_t0 = (t-iter0-1)*m_Ts; //t-1 because the increment appears above
        awake_time = t0 + awake_at_ms_since_t0;
        std::this_thread::sleep_until(awake_time);
        
        m_x_meas = x_next; //cpp simulation (uncomment this line in experiments)
        
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
    m_admm->lbA.block(m_meta.idx_eqx0,0,m_meta.nx,1) = m_x_meas;
    m_admm->ubA.block(m_meta.idx_eqx0,0,m_meta.nx,1) = m_x_meas;
    return;
}

void dmpcAgent::set_u0_in_ocp(){
    m_p.segment(m_meta.nx,m_meta.nu) = m_u0;
    m_admm->lbA.block(m_meta.idx_equ0,0,m_meta.nu,1) = m_u0;
    m_admm->ubA.block(m_meta.idx_equ0,0,m_meta.nu,1) = m_u0;
    return;
}

void dmpcAgent::get_u_from_sol(){
    m_u0 = m_sol.block(m_meta.idx_u0,0,m_meta.nu,1);
    m_u1 = m_sol.block(m_meta.idx_u1,0,m_meta.nu,1);
    m_u2 = m_sol.block(m_meta.idx_u2,0,m_meta.nu,1);
}

void dmpcAgent::warm_start_solver(){

    // m_zbar = m_admm->z_bar;
    // m_gam = m_admm->gam;

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
}

void dmpcAgent::receive_xmeas(int t, int state_step, bool experimental_mode, double waiting_timeout) {
    bool waitingTimeOver= false;
    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
    double elapsedTime = ((double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()) / 1000000.0;
    bool receivedAllStates = true;
    while ((stateVectorHandler->StateSeq < (state_step * (t - 1)+1)) && ((!experimental_mode) || (!waitingTimeOver))) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        if (!experimental_mode) {
            comms->publish(controlInputTopic, &controlInputMsg);
        }
        if (t > 1) {
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
    }

    if (experimental_mode) {
        stateVectorHandler->processMessages(false);
    } else {
        stateVectorHandler->processMessages(true);
    }
    if (!receivedAllStates) {
        std::cout << "Did not receive all states within waiting_timeout s. stateVectorHandler->StateSeq = " << stateVectorHandler->StateSeq << ", state_step * (t - 1) = " << (state_step * (t - 1)+1) << " (should be equal)" << std::endl;
    } else {
        std::cout << "Received all states." << std::endl;
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
    m_admm->init_comms();

    m_z     = VectorXd::Zero(m_admm->nx);
    m_zbar  = VectorXd::Zero(m_admm->nx);
    m_gam   = VectorXd::Zero(m_admm->nx);
    m_nu    = VectorXd::Zero(m_admm->ng);
    m_mu    = VectorXd::Zero(m_admm->nh);

    int flag = m_admm->solve(10,m_z,m_zbar,m_nu,m_mu,m_gam,false); //initializes data structures in admmAgent

    //reset all ADMM variables to zero. 
    m_z     = VectorXd::Zero(m_admm->nx);
    m_zbar  = VectorXd::Zero(m_admm->nx);
    m_gam   = VectorXd::Zero(m_admm->nx);
    m_nu    = VectorXd::Zero(m_admm->ng);
    m_mu    = VectorXd::Zero(m_admm->nh);
    return 1;
}

int dmpcAgent::update_setpoint_in_ocp(){

    m_p.segment(0,m_meta.nx) = m_x_meas;
    m_p.segment(m_meta.nx,m_meta.nu) = m_u0;
    m_p.segment(m_meta.nx+m_meta.nu,m_meta.nxd) = m_xd;

    std::string functionLibrary = m_folderName + "/locFuns.so";

    casadi::DM z_cas;
    casadi::DM p_cas;

    // Use CasADi's "external" to load the compiled function
    casadi::Function f; // = casadi::external("gradFun1","sProb_chain/locFuns.so");
    std::string str;

    std::vector<casadi::DM> arg; // = {z};
    std::vector<casadi::DM> res; // = f(arg);

    int nx = m_admm->g.size();


    VectorXd z = Eigen::VectorXd::Zero(nx);
  
    z_cas = Eigen2casadi(z);
    p_cas = Eigen2casadi(m_p);

    //g
    str = "gradFun" + std::to_string(m_my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    arg = {z_cas,p_cas};
    res = f(arg);
    m_admm->m_sprob.g[m_my_id]= casadi2EigenVector(res[0]);
    
    //beq
    str = "eqfun" + std::to_string(m_my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    // arg = {z_cas};
    res = f(arg);
    int ng = res[0].size1();
    int cols = res[0].size2();
    m_admm->m_sprob.beq[m_my_id] = Eigen::VectorXd::Zero(ng);        
    std::memcpy(m_admm->m_sprob.beq[m_my_id].data(), res.at(0).ptr(), sizeof(double)*ng*cols);
    m_admm->m_sprob.beq[m_my_id] = -m_admm->m_sprob.beq[m_my_id];

    m_admm->update_g_beq();
   

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
    file << "mpc step, mpc step time (us), x0(1), x0(2), u0(1), u0(2), xd(0), xd(1), receivex time (us), admm time (us), sleep time (us), sendu time (us), admm iter, admm iter time (us), loc_qp time(us), zcomm time (us), zbarcomm time(us), sendvin time (us), receivevout time (us)\n";
    file.close();
    }
    

    int N_rows = m_admm->iter_timer.m_log.size();
    int k = 1; //MPC step
    unsigned int row = 0;
    while (row < N_rows){
        for (unsigned int i = 0; i < m_maxiter; i++){       
            file.open(fileName.str(),std::ios_base::app);
            if (file.is_open())
            {
                file << k << "," << m_mpc_timer.m_log[k-1] << "," << m_x_log(k-1,0) << "," << m_x_log(k-1,1) << "," << m_u_log(k-1,0) << "," << m_u_log(k-1,1) << "," << m_xd_log(k-1,0) << "," << m_xd_log(k-1,1) << "," << m_receivex_timer.m_log[k-1] << "," << m_admm_timer.m_log[k-1] << "," << m_sleep_timer.m_log[k-1] << "," << m_sendu_timer.m_log[k-1] << "," << i << "," << m_admm->iter_timer.m_log[row] << "," << m_admm->loc_timer.m_log[row] << "," << m_admm->z_comm_timer.m_log[row] << "," << m_admm->zbar_comm_timer.m_log[row] << "," << m_admm->send_vin_timer.m_log[row] << "," << m_admm->receive_vout_timer.m_log[row] << "\n";        
            }
            file.close();

            row++;
        }
        k++;
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

void dmpcAgent::clear_time_measurements(){
    m_mpc_timer.clear();
    m_receivex_timer.clear();
    m_sendu_timer.clear();
    m_admm_timer.clear();
    m_sleep_timer.clear();
    return;
}