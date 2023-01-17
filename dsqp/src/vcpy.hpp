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

#ifndef VCPY_H
#define VCPY_H

#include <Eigen/Dense>
#include <map>


class vCpy{
    public:
        Eigen::Vector<double, Eigen::Dynamic> val;
        Eigen::Vector<int, Eigen::Dynamic> og_idx;
        Eigen::Vector<int, Eigen::Dynamic> cpy_idx;
        std::map<int,int> og_idx_to_cpy_idx;
        int original_agent;
        int copying_agent;
        int nv;

        vCpy(){
            nv = 0;
            original_agent = 0;
            copying_agent = 0;
            val = Eigen::VectorXd(0);
            og_idx = Eigen::VectorXi(0);
            cpy_idx = Eigen::VectorXi(0);
        }
};


#endif