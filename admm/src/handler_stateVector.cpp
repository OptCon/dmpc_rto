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

#include "handler_stateVector.hpp"

void handler_stateVector::handleMessage(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const itmessage::vector_t* msg) {
    // handle the robots
    // check for newness of message
    // in seq_state_aggregators, the first entry is for the object, the rest for the robots
    // careful: msg->id_sender uses Matlab numbering, counting from 1, and not 0
	//std::cout << "msg->seq_number = " << msg->seq_number << std::endl;
    // if (msg->seq_number > seq_state_aggregators(msg->id_sender-1)) {
    //     std::cout << "Received new message with seq_number = = " << msg->seq_number  << std::endl;
    // }
    if (msg->seq_number >= (seq_state_aggregators(msg->id_sender-1) + state_step)) {
        mutexes_stateVector[msg->id_sender-1].lock();
        msg_cache[msg->id_sender-1].id_sender = msg->id_sender;
        msg_cache[msg->id_sender-1].length = msg->length;
        msg_cache[msg->id_sender-1].seq_number = msg->seq_number;
        msg_cache[msg->id_sender-1].value = msg->value;
        cachedMsgNew[msg->id_sender-1] = true;
        seq_state_aggregators(msg->id_sender-1) = msg->seq_number;
        StateSeq = seq_state_aggregators.minCoeff();
        mutexes_stateVector[msg->id_sender-1].unlock();
	}
}

void handler_stateVector::processMessages(bool strict) {
    for (int idx_cache=0; idx_cache<(N); idx_cache++) {
        mutexes_stateVector[idx_cache].lock();
        if (likely(cachedMsgNew[idx_cache])) {
            // handle the robots
            // in seq_state_aggregators, the first entry is for the object, the rest for the robots
            // careful: msg->id_sender / msg_cache[idx_cache].id_sender uses Matlab numbering, counting from 1, and not 0
	        Eigen::VectorXd value(msg_cache[idx_cache].length);
	        value = VectorXd::Map(msg_cache[idx_cache].value.data(), msg_cache[idx_cache].value.size());
            value(3) = value(3)*M_PI/180.0; // (value(3) received in degree, convert to radian)
            // consider all agents:
            robot_poses.col(msg_cache[idx_cache].id_sender-1).setZero();
            robot_poses.col(msg_cache[idx_cache].id_sender-1).head<2>() = value.head<2>();
            if ((value(3) < -M_PI / 2.0) && (lastRawAngleRob(msg_cache[idx_cache].id_sender - 1) > M_PI / 2.0)) {
                numRotsRob(msg_cache[idx_cache].id_sender - 1) = numRotsRob(msg_cache[idx_cache].id_sender - 1) + 1;
            }
            if ((value(3) > M_PI / 2.0) && (lastRawAngleRob(msg_cache[idx_cache].id_sender - 1) < -M_PI / 2.0)) {
                numRotsRob(msg_cache[idx_cache].id_sender - 1) = numRotsRob(msg_cache[idx_cache].id_sender - 1) - 1;
            }
            lastRawAngleRob(msg_cache[idx_cache].id_sender - 1) = value(3);
            robot_poses(2,msg_cache[idx_cache].id_sender-1) = value(3) + 2.0 * M_PI * ((double)numRotsRob(msg_cache[idx_cache].id_sender - 1));

            cachedMsgNew[idx_cache] = false;
            seq_state_aggregators_processed(msg_cache[idx_cache].id_sender-1) = msg_cache[idx_cache].seq_number;
            StateSeq_processed = seq_state_aggregators_processed.minCoeff();
        } else {
            std::string errorMsg = "";
            errorMsg.append(std::string(__FILE__));
            errorMsg.append(" around line ");
            errorMsg.append(std::to_string(__LINE__));
            errorMsg.append(": ");
            if (strict) {
                errorMsg = "Error in " + errorMsg;
                errorMsg.append(" Programming error, expected message to be new.");
                throw std::runtime_error(errorMsg);
            }
        }
        mutexes_stateVector[idx_cache].unlock();
    }
}