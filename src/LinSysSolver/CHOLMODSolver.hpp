//
//  CHOLMODSolver.hpp
//  IPC
//
//  Created by Minchen Li on 6/22/18.
//
#pragma once

#ifdef IPC_WITH_CHOLMOD

#include "LinSysSolver.hpp"

#include "cholmod.h"

#include <Eigen/Eigen>

#include <vector>
#include <set>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
class CHOLMODSolver : public LinSysSolver<vectorTypeI, vectorTypeS> {
    typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:
    cholmod_common cm;
    cholmod_sparse* A;
    cholmod_factor* L;
    cholmod_dense *b, *solution;
    cholmod_dense *x_cd, *y_cd; // for multiply

    void *Ai, *Ap, *Ax, *bx, *solutionx, *x_cdx, *y_cdx;

public:
    CHOLMODSolver(void);
    ~CHOLMODSolver(void);

    LinSysSolverType type() const override { return LinSysSolverType::CHOLMOD; }

    void set_pattern(const std::vector<std::set<int>>& vNeighbor, const std::set<int>& fixedVert) override;
    void set_pattern(const Eigen::SparseMatrix<double>& mtr) override; //NOTE: mtr must be SPD
    void load(const char* filePath, Eigen::VectorXd& rhs) override;

    void analyze_pattern(void) override;

    bool factorize(void) override;

    void solve(Eigen::VectorXd& rhs,
        Eigen::VectorXd& result) override;
    void solve(Eigen::VectorXd& rhs,
        Eigen::VectorXd& result, const std::string& outputPath) override;

    virtual void multiply(const Eigen::VectorXd& x,
        Eigen::VectorXd& Ax) override;

    virtual void outputFactorization(const std::string& filePath) override;
    void outputA(const std::string& filePath) override;
    void outputSolve(const std::string& filePath) override;
};

} // namespace IPC

#endif
