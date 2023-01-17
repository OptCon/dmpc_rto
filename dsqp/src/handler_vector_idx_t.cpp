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

#include "handler_vector_idx_t.hpp"

void handler_vector_idx_t::handleMessage(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const itmessage::vector_idx_t* msg) {
	if ((msg->seq_number == (seq_processed+1) && (msg->seq_number == (seq_received+1))) && ((msg->id_sender == id_number) || (msg->id_sender == -1))) {
		mutex_vector_t.lock();
        msg_cache.val_length = msg->val_length;
        msg_cache.idx_length = msg->idx_length;
        msg_cache.value = msg->value;
        msg_cache.idx = msg->idx;
        msg_cache.id_sender = msg->id_sender;
        msg_cache.seq_number = msg->seq_number;

        cached_msg_new = true;
        seq_received = msg->seq_number;
        val_length = msg->val_length;
        idx_length = msg->idx_length;
        mutex_vector_t.unlock();
	} 
    else if (chan.find("sync") != std::string::npos) { //s.find('a') == std::string::npos
        if ( (msg->id_sender == id_number) || (msg->id_sender == -1) ) {
		mutex_vector_t.lock();
        msg_cache.val_length = msg->val_length;
        msg_cache.idx_length = msg->idx_length;
        msg_cache.value = msg->value;
        msg_cache.idx = msg->idx;
        msg_cache.id_sender = msg->id_sender;
        msg_cache.seq_number = msg->seq_number;

        cached_msg_new = true;
        seq_received = msg->seq_number;
        val_length = msg->val_length;
        idx_length = msg->idx_length;
        mutex_vector_t.unlock();
        }
    }
}

bool handler_vector_idx_t::processMessage() {
    mutex_vector_t.lock();

    if (cached_msg_new) {
        value.resize(msg_cache.val_length);
		value = Eigen::VectorXd::Map(msg_cache.value.data(), msg_cache.value.size());
        idx.resize(msg_cache.idx_length);
        idx = Eigen::VectorXi::Map(msg_cache.idx.data(), msg_cache.idx.size());
        id_number = msg_cache.id_sender;
        cached_msg_new = false;
        seq_processed = msg_cache.seq_number;
        mutex_vector_t.unlock();
        return true;
    } else {
        mutex_vector_t.unlock();
        return false;
    }
}
