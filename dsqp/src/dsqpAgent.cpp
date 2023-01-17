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



#include "dsqpAgent.hpp"

dsqpAgent::dsqpAgent ( const std::string& folderName_, int my_id_, unsigned int Nagents_, double rho_, unsigned int sqp_maxiter_, unsigned int admm_maxiter_,lcm::LCM* lcm) : m_qp(Nagents_){

    my_id = my_id_;
    Nagents = Nagents_;
    m_folderName = folderName_;

    m_qp.read_AA(m_folderName,Nagents);

    nx = m_qp.A[my_id].cols();

    std::cout << "dsqpAgent::nx = " << nx << std::endl;

    z = VectorXd::Zero(nx);
    gam = VectorXd::Zero(nx);
    z_ = VectorXd::Zero(nx);
    delz = VectorXd::Zero(nx);

    m_admm = nullptr;

    m_lcm = lcm;
    m_sqp_maxiter = sqp_maxiter_;
    m_admm_maxiter = admm_maxiter_;
    m_rho = rho_;

    functionLibrary = m_folderName + "/locFuns.so";  

}

int dsqpAgent::init(){

    nx = m_qp.A[my_id].cols();
    Eigen::VectorXd ztmp = Eigen::VectorXd::Zero(nx);
    Eigen::VectorXd* nutmp = nullptr;
    Eigen::VectorXd* mutmp = nullptr;
    init(&ztmp, nutmp, mutmp);

    return 0;
}

int dsqpAgent::init(Eigen::VectorXd* z_, Eigen::VectorXd* nu_, Eigen::VectorXd* mu_){

    m_qp.read_AA(m_folderName,Nagents);
    m_qp.read_ublb(m_folderName,Nagents);
    
    for (int k = 0; k < m_qp.lb[my_id].size(); k++){
        if (m_qp.lb[my_id][k] == -casadi::inf){
            m_qp.lb[my_id][k] = -pow(10,20);
        } else if (m_qp.lb[my_id][k] == casadi::inf){
            m_qp.lb[my_id][k] = pow(10,20);
        }
    }

    for (int k = 0; k < m_qp.ub[my_id].size(); k++){
        if (m_qp.ub[my_id][k] == -casadi::inf){
            m_qp.ub[my_id][k] = -pow(10,20);
        } else if (m_qp.ub[my_id][k] == casadi::inf){
            m_qp.ub[my_id][k] = pow(10,20);
        }
    }

    ub_file = m_qp.ub[my_id];
    lb_file = m_qp.lb[my_id];

    bool GN = true;

    buildQP(z_,nu_,mu_,GN,true);
   
    z = *z_;
    
    m_admm = new admmAgent(m_qp,my_id,m_rho,m_lcm);

    std::cout << "starting init_comms" << std::endl;
    m_admm->init_comms();
    std::cout << "init_comms is finished" << std::endl;

    return 0;

}

dsqpAgent::~dsqpAgent(){
    if (m_admm != nullptr){
        delete m_admm;
        m_admm = nullptr;
    }
}

int dsqpAgent::buildQP(Eigen::VectorXd* z_, Eigen::VectorXd* nu, Eigen::VectorXd* mu, bool GN, bool eval_HessF){

    if(z_ == nullptr){
        z = Eigen::VectorXd::Zero(nx);
    } else {
        z = *z_;
    }

    z_cas = Eigen2casadi(z);
    p_cas = Eigen2casadi(m_p);

    //g
    str = "gradFun" + std::to_string(my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    arg = {z_cas,p_cas};
    res = f(arg);
    m_qp.g[my_id] = Eigen::VectorXd::Zero(nx);
    m_qp.g[my_id] = casadi2EigenVector(res[0]);

    //Aeq
    str = "JGfun" + std::to_string(my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    res = f(arg);
    casadi::DM Aeq = res[0];
    m_qp.Aeq[my_id] = casadi2Eigen(Aeq);

    //Aineq
    str = "JHfun" + std::to_string(my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    res = f(arg);
    casadi::DM Aineq = res[0];
    m_qp.Aineq[my_id] = casadi2Eigen(Aineq);

    //beq
    str = "eqfun" + std::to_string(my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    res = f(arg);
    ng = res[0].size1();
    int cols = res[0].size2();
    m_qp.beq[my_id] = Eigen::VectorXd::Zero(ng);        
    std::memcpy(m_qp.beq[my_id].data(), res.at(0).ptr(), sizeof(double)*ng*cols); //works
    m_qp.beq[my_id] = -m_qp.beq[my_id];

    //bineq
    str = "ineqfun" + std::to_string(my_id+1); //convert to matlab index
    f = casadi::external(str,functionLibrary);
    res = f(arg);
    nh = res[0].size1();
    cols = res[0].size2();
    m_qp.bineq[my_id] = Eigen::VectorXd::Zero(nh);        
    std::memcpy(m_qp.bineq[my_id].data(), res.at(0).ptr(), sizeof(double)*nh*cols); //works
    m_qp.bineq[my_id] = -m_qp.bineq[my_id];

    if(GN == true){
        if (eval_HessF == true){
            //HessF
            if(nu!=nullptr){
                nu_cas = Eigen2casadi(*nu);
            } else {
                nu_cas = casadi::DM::zeros(ng,1);
            }
            if(mu!=nullptr){
                mu_cas = Eigen2casadi(*mu);
            } else {
                mu_cas = casadi::DM::zeros(nh,1);
            }
            str = "HessFfun" + std::to_string(my_id+1); //convert to matlab index
            f = casadi::external(str,functionLibrary);
            arg = {z_cas,p_cas};

            res = f(arg);
            
            casadi::DM HL = res[0];
            MatrixXd HLeig = casadi2Eigen(HL);
            m_qp.H[my_id] = HLeig;
        }

    } else {
        //HessL
        if(nu!=nullptr){
            nu_cas = Eigen2casadi(*nu);
        } else {
            nu_cas = casadi::DM::zeros(ng,1);
        }
        if(mu!=nullptr){
            mu_cas = Eigen2casadi(*mu);
        } else {
            mu_cas = casadi::DM::zeros(nh,1);
        }
        str = "HessLfun" + std::to_string(my_id+1); //convert to matlab index
        f = casadi::external(str,functionLibrary);
        arg = {z_cas,nu_cas,mu_cas,p_cas};

        res = f(arg);
        
        casadi::DM HL = res[0];
        MatrixXd HLeig = casadi2Eigen(HL);
        m_qp.H[my_id] = HLeig;
        m_qp.H[my_id] = 0.5*(m_qp.H[my_id] + m_qp.H[my_id].transpose());
    }

    //regularization
    reg_timer.tic();
    if(GN == false){
        SelfAdjointEigenSolver<MatrixXd> es(m_qp.H[my_id]);
        std::cout << "Eigenvalues before regularization:" << std::endl;
        MatrixXd D_tmp = es.eigenvalues().asDiagonal();
        for (int j = 0; j < D_tmp.cols(); j ++){
            std::cout << D_tmp(j,j) << std::endl;
        }
        
        if (es.eigenvalues()[0] <= -0.0){//pow(10,-8)){

            double reg = pow(10,-4);
            MatrixXd D = es.eigenvalues().asDiagonal();
            for (int j = 0; j < D.cols(); j++){
                if (D(j,j) <= 0.0){
                    D(j,j) = -D(j,j);
                } 
                if (D(j,j) <= reg){
                    D(j,j) = reg;
                }
            }

            MatrixXd V = es.eigenvectors();

            m_qp.H[my_id] = (V * D * V.transpose()).real();
            m_qp.H[my_id] = 0.5*(m_qp.H[my_id] + m_qp.H[my_id].transpose());
        }
        SelfAdjointEigenSolver<MatrixXd> es_tmp(m_qp.H[my_id]);
        std::cout << "Eigenvalues after regularization:" << std::endl;
        D_tmp = es_tmp.eigenvalues().asDiagonal();
        for (int j = 0; j < D_tmp.cols(); j ++){
            std::cout << D_tmp(j,j) << std::endl;
        }
    }
    reg_timer.toc();

    //lb,ub
    for (int k = 0; k < nx; k++){
        if (lb_file[k] == -pow(10,20)){
            m_qp.lb[my_id][k] = -pow(10,20); //DBL_MAX;
        }else{
            m_qp.lb[my_id][k] = lb_file[k] - z[k];
        }
        if (ub_file[k] == pow(10,20)){
            m_qp.ub[my_id][k] = pow(10,20);
        }else{
            m_qp.ub[my_id][k] = ub_file[k] - z[k];
        }        
    }
    return 0;
}

int dsqpAgent::solve(bool createLog_){
    int flag = solve(m_sqp_maxiter,m_admm_maxiter, createLog_);
    return flag;
}

int dsqpAgent::solve(int sqp_maxiter_, int admm_maxiter_, bool createLog_){

    if(m_admm == nullptr){
        return 1;
    }

    int flag = 0;
    bool GN = true;             //use Gauss-Newton for H
    bool eval_HessF = false;    //if Gauss-Newton is used, recalculate HessF
    bool createLog = false;
    if (createLog){
        z_log = MatrixXd::Zero(nx,sqp_maxiter_);
    }

    for (unsigned int iter = 0; iter < sqp_maxiter_; iter++){

        iter_timer.tic();
        delz = VectorXd::Zero(nx);
        z_ = VectorXd::Zero(nx);

        buildQP_timer.tic();
        flag = buildQP(&z,&nu,&mu,GN,eval_HessF);
        buildQP_timer.toc();

        updateQP_timer.tic();
        flag = m_admm->updateQP(m_qp,!GN);
        updateQP_timer.toc();

        admm_timer.tic();
        flag = m_admm->solve(admm_maxiter_, z_, delz, nu, mu, gam, createLog_);
        admm_timer.toc();

        z = z + delz;

        if(createLog){
            z_log.block(0,iter,nx,1) = z;
        }

        iter_timer.toc();     

    }

    if (createLog){
        std::time_t t = std::time(0);   // get time now
        std::tm* now = std::localtime(&t);
    
        std::ostringstream fileName;
        fileName << "dsqp_zlog" << "_agent" << my_id << "_" << (now->tm_year + 1900) << '_' << (now->tm_mon + 1) << '_' <<  now->tm_mday << "_" << now->tm_hour << "_" << now->tm_min << "_" << now->tm_sec <<".csv";
        std::ofstream file(fileName.str());
        if (file.is_open())
        {
        file << z_log << std::endl; //"mpc step, mpc step time (us), x0(1), x0(2), u0(1), u0(2), xd(0), xd(1), receivex time (us), admm time (us), sleep time (us), sendu time (us), admm iter, admm iter time (us), loc_qp time(us), zcomm time (us), zbarcomm time(us), sendvin time (us), receivevout time (us)\n";
        file.close();
        }
    }
    return 0;

}

// Dense to Dense
Eigen::MatrixXd dsqpAgent::casadi2Eigen ( const casadi::DM& A ){
    // This method is based on code by Petr Listov, see https://groups.google.com/g/casadi-users/c/npPcKItdLN8    
    casadi::Sparsity SpA = A.get_sparsity();
    std::vector<long long int> output_row, output_col;
    SpA.get_triplet(output_row, output_col);
    std::vector<double> values = A.get_nonzeros();
    using T = Eigen::Triplet<double>;
    std::vector<T> TripletList;
    TripletList.resize(values.size());
    for(int k = 0; k < values.size(); ++k){
        TripletList[k] = T(output_row[k], output_col[k], values[k]);
    }
    Eigen::SparseMatrix<double> SpMatrx(A.size1(), A.size2());
    SpMatrx.setFromTriplets(TripletList.begin(), TripletList.end());

    Eigen::MatrixXd DeMatrx = SpMatrx.toDense();
    return DeMatrx;
}

// Dense to Dense
Eigen::VectorXd dsqpAgent::casadi2EigenVector ( const casadi::DM& A ){
    // This method is based on code by Petr Listov, see https://groups.google.com/g/casadi-users/c/npPcKItdLN8    
    casadi::Sparsity SpA = A.get_sparsity();
    std::vector<long long int> output_row, output_col;
    SpA.get_triplet(output_row, output_col);
    std::vector<double> values = A.get_nonzeros();
    using T = Eigen::Triplet<double>;
    std::vector<T> TripletList;
    TripletList.resize(values.size());
    for(int k = 0; k < values.size(); ++k){
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

    return DeVec;
}

casadi::DM dsqpAgent::Eigen2casadi( const Eigen::VectorXd& in){
    int length = in.size();
    casadi::DM out = casadi::DM::zeros(length,1);
    std::memcpy(out.ptr(), in.data(), sizeof(double)*length*1);
    return out;
}

void dsqpAgent::clear_time_measurements(){
    iter_timer.clear();
    buildQP_timer.clear();
    reg_timer.clear();
    updateQP_timer.clear();
    admm_timer.clear();
    m_admm->clear_time_measurements();
    return;
}

void dsqpAgent::reserve_time_measurements(unsigned int new_cap){
    buildQP_timer.reserve(new_cap);
    iter_timer.reserve(new_cap);
    admm_timer.reserve(new_cap);
    reg_timer.reserve(new_cap);
    updateQP_timer.reserve(new_cap);
}