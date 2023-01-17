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

#ifndef HANDLER_VECTOR_IDX_T_H
#define HANDLER_VECTOR_IDX_T_H

#include <Eigen/Dense>
#include <lcm/lcm-cpp.hpp>
#include <limits>
#include <mutex>
#include <vector>
#include <iostream>
#include <atomic>

#include "lcm/itmessage/vector_idx_t.hpp"

class handler_vector_idx_t {
   public:
    std::atomic<std::int64_t> seq_received;
    std::atomic<std::int64_t> seq_processed;
    int id_number;  // to compare with id_sender, do not compare if equal to -1
   private:
    itmessage::vector_idx_t msg_cache; // most current data, updated by handleMessage
   public:
    Eigen::VectorXd value; // most current processed data - written to from the msg_cache by processMessage
    Eigen::VectorXi idx;
    int val_length = -1;
    int idx_length = -1;
    std::atomic<bool> cached_msg_new;
    std::mutex mutex_vector_t;

    void handleMessage(const lcm::ReceiveBuffer* rbuf,
                       const std::string& chan,
                       const itmessage::vector_idx_t* msg);

    bool processMessage();

    handler_vector_idx_t(int64_t seq_init, int agentNumber_) {
        cached_msg_new = false;
        seq_received = seq_init;
        seq_processed = seq_init;
        id_number = agentNumber_;
    }

    ~handler_vector_idx_t() {}
};

#endif
