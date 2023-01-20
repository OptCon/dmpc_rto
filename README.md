# dmpc_rto

## Distributed Model Predictive Control via Decentralized Real Time Optimization
Efficient implementations of ADMM for quadratic programming and of d-SQP for nonlinear programming.
Please cite https://arxiv.org/abs/2301.07960 if you use this code.

## Build Instructions

0. Install dependencies 
	- lcm (https://lcm-proj.github.io/, https://github.com/lcm-proj)
	- CasADi (www.casadi.org, https://github.com/casadi/casadi)
	- qpOASES (https://github.com/coin-or/qpOASES)
	- Eigen (https://eigen.tuxfamily.org/index.php?title=Main_Page)
	- Hint: ensure that lcm is on your Java classpath, e.g. by adding `export CLASSPATH=${CLASSPATH}:/usr/local/share/java/lcm.jar` to `~/.bashrc`
	
1. Go either to folder `admm/` or to `dsqp/`.

2. `mkdir build`

3. `cd build/`

4. `cmake ../src`
	OR, if you have non-standard installation directories for
	dependencies, use the two commands
	```
		export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
		cmake -DLCM_INCLUDE_DIR=/usr/local/include -DLCM_LIBRARY_DIR=/usr/local/lib -DQPOASES_LIBRARY_DIR=/usr/local/lib -DQPOASES_INCLUDE_DIR=/usr/local/include ../src
	```
	and modify the paths. The paths in the example are the default 
	locations, otherwise
	 - headers of lcm should be in the directory `${LCM_INCLUDE_DIR}/lcm`,
	 - `liblcm.so` should be located in the directory `${LCM_LIBRARY_DIR}`,
	 - `libqpOASES.so` should be located in the directory `${QPOASES_LIBRARY_DIR}`,
	 - pkgconfig files for casadi and related packages need to be in the directory `${PKG_CONFIG_PATH}`, i.e., folder needs to contain `casadi.pc`, `coinmetis.pc`, `conmumps.pc`, `ipopt.pc`, `opencv.pc`).

5. Execute `bash gen_locFuns.sh` in `bin/sProb_chain_QP_N7` or `bin/sProb_chain_QCQP_N7`.

6. Execute `make` or `make -j #jobs` in `build`.

7. Enjoy - the `test_simulation.sh` script in the `bin/` directory gives startup advice. You may inspect network traffic with LCM's `lcm-spy`.
	- Hint: to run `lcm-spy`, ensure that the `itmessage.jar` (generated during the build process in `src/lcm` is on your Java classpath, e.g. by adding `export CLASSPATH=${CLASSPATH}:/usr/local/share/java/lcm.jar:/path/to/itmessage.jar` to `~/.bashrc`
