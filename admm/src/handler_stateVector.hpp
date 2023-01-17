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

#ifndef _USE_MATH_DEFINES 
#define _USE_MATH_DEFINES 
#endif
#ifndef HANDLER_STATEVECTOR_H
#define HANDLER_STATEVECTOR_H  

#include <Eigen/Dense>
#include <vector>
#include <mutex>
#include <iostream>
#include <lcm/lcm-cpp.hpp>
#include "lcm/itmessage/vector_t.hpp"

using namespace Eigen;

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

class handler_stateVector {
	public:
		//itmessage::vector_t msgContent;
		int64_t StateSeq;
        int64_t StateSeq_processed;
		Matrix<int64_t, Dynamic, 1> seq_state_aggregators;
        Matrix<int64_t, Dynamic, 1> seq_state_aggregators_processed;
		int agentNumber; // matlab numbering! (starting with 1, not 0)
		Matrix<double, 3, Dynamic> robot_poses; // first two lines: position; third line: orientation

        int64_t StateSeq_init;
		
        int numRots; // signed
        double lastRawAngle; // between -pi and pi

		int N; // total number of robots
        
        VectorXi numRotsRob;
        VectorXd lastRawAngleRob;

        int state_step; // only every state_step-th state vector is considered (for 'supersampling' measurements)

        itmessage::vector_t* msg_cache;
        bool* cachedMsgNew;

        std::mutex* mutexes_stateVector;

		
		handler_stateVector(int64_t StateSeq_init_, int agentNumber_, int N_, int state_step_) {
			StateSeq_init = StateSeq_init_;
			StateSeq = StateSeq_init;
            StateSeq_processed = StateSeq_init;
			agentNumber = agentNumber_;
            numRots = 0;
            lastRawAngle = 0.0;
			N = N_;
            state_step = state_step_;

            lastRawAngleRob.resize(N);
            lastRawAngleRob.setZero();
            numRotsRob.resize(N);
            numRotsRob.setZero();
            robot_poses.resize(3,N);
            robot_poses.fill(std::numeric_limits<double>::quiet_NaN());

			seq_state_aggregators.resize(N,1);
			seq_state_aggregators.fill(StateSeq_init);
            seq_state_aggregators_processed.resize(N,1);
			seq_state_aggregators_processed.fill(StateSeq_init);

            msg_cache = new itmessage::vector_t[N];
            cachedMsgNew = new bool[N+1];
            mutexes_stateVector = new std::mutex[N];
            for (int i=0; i<N; i++) {
                cachedMsgNew[i] = false;
            }
		};
		
		void handleMessage(const lcm::ReceiveBuffer* rbuf,
				const std::string& chan,
				const itmessage::vector_t* msg);

        void processMessages(bool strict); 
		
		~handler_stateVector() {
            delete[] msg_cache;
            delete[] cachedMsgNew;
            delete[] mutexes_stateVector;
        };
};

#endif
