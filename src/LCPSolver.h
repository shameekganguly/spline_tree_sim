// LCPSolver.h: Solves the problem of form a = Af + b, a > 0 f > 0, a_i*f_i = 0
// NOTE: For now, we assume that A is positive semi-definite
// NOTE: The implementation is based on Baraff's algorithm for frictionless
// contacts 

#ifndef LCP_SOLVER_H
#define LCP_SOLVER_H

#include <Eigen/Core>
#include <unordered_set>
#include <vector>

// NOTE: class uses fixed size of maximum contacts
#define MAX_LCP_SIZE 64

// NOTE: We currently do not warm start. TODO: consider warm starting as an optimization
class LCPSolver {
	// typedefs
	typedef std::unordered_set<uint> IndexSet;
	typedef std::vector<uint> IndexVec;
	typedef Eigen::Matrix<double, MAX_LCP_SIZE, 1>::Index LCPVectorIndex;
	typedef Eigen::Matrix<double, MAX_LCP_SIZE, MAX_LCP_SIZE>::Index LCPMatIndex;

	// static member functions
public:
	// max problem size
	static uint maxProblemSize() {
		return MAX_LCP_SIZE;
	}

	// regular member functions
public:
	// ctor
	LCPSolver() {
		/* Nothing to do */
	}

	// dtor
	virtual ~LCPSolver() {
		/* Nothing to do */
	}

	// solve
	void solve(Eigen::VectorXd& ret_a, Eigen::VectorXd& ret_f, const Eigen::MatrixXd& A, const Eigen::VectorXd& b);

	// internal functions
public:
	// initialize
	void initialize();

	// solveInternal. called after initialization is complete
	void solveInternal();

	// drive-to-zero
	void driveToZero(uint ind_d, const double tol=1e-3);

	// fdirection
	void fdirection(uint ind_d);

	// move an index from C to NC
	// also reorder the _A, _b, _f and _a matrices to maintain the index order C, NC, d, others
	void moveCtoNC(uint ind);

	// move an index from NC to C
	// also reorder the _A, _b, _f and _a matrices to maintain the index order C, NC, d, others
	void moveNCtoC(uint ind);

	// add from unknown to C
	// also reorder the _A, _b, _f and _a matrices to maintain the index order C, NC, d, others
	void addToC(uint ind);

	// add from unknown to NC
	// also reorder the _A, _b, _f and _a matrices to maintain the index order C, NC, d, others
	void addToNC(uint ind);

	// swap indices in _A, _b, _f and _a. Also update the index order map for bookkeeping
	void swap(uint ind1, uint ind2);

	// data members
public:
	// A mat
	Eigen::Matrix<double, MAX_LCP_SIZE, MAX_LCP_SIZE> _A;

	// b vec
	Eigen::Matrix<double, MAX_LCP_SIZE, 1> _b;

	// f, delta_f vec
	Eigen::Matrix<double, MAX_LCP_SIZE, 1> _f;
	Eigen::Matrix<double, MAX_LCP_SIZE, 1> _delta_f;

	// a, delta_a vec
	Eigen::Matrix<double, MAX_LCP_SIZE, 1> _a;
	Eigen::Matrix<double, MAX_LCP_SIZE, 1> _delta_a;

	// current problem size
	uint _size;

	// index management data structures
	IndexVec _order_map;
	uint _count_C; // number of indices in known contact list
	uint _count_NC; // number of indices in known non-contact list

};

#endif //LCP_SOLVER_H
