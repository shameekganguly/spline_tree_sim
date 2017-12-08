
#include "LCPSolver.h"
#include <Eigen/Cholesky>
#include <stdexcept>

void LCPSolver::solve(Eigen::VectorXd& ret_a, Eigen::VectorXd& ret_f, const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
	// set size and check if less than our limit
	_size = b.rows();
	if(_size > MAX_LCP_SIZE) {
		throw(std::runtime_error("LCP size exceeds max internal limit."));
	}

	ret_a.setZero(_size);
	ret_f.setZero(_size);

	// copy _A, _b
	_A.block(0,0,_size,_size) = A;
	_b.head(_size) = b;

	// initialize
	initialize();

	// solve
	solveInternal();

	// copy return values
	ret_a = _a.head(_size);
	ret_f = _f.head(_size);
}

// ----- internal functions -----

// initialize
void LCPSolver::initialize() {
	// set _a to _b
	_a.head(_size) = _b.head(_size);

	// set _f to 0
	_f.head(_size).setZero();
	_delta_f.head(_size).setZero();
	_delta_a.head(_size).setZero();

	// reset counts
	_count_C = 0;
	_count_NC = 0;

	// initialize index vector for bookkeeping
	_order_map.clear();
	for(uint i=0; i<_size; ++i) {
		_order_map.push_back(i);
	}
}

// internal solve function
void LCPSolver::solveInternal() {
	double max_pen_acc;
	while (true) {
		if (_size-_count_C-_count_NC == 0) {
			// all contacts have been processed
			break;
		}
		LCPVectorIndex index;
		max_pen_acc = _a.segment(_count_C+_count_NC, _size-_count_C-_count_NC).minCoeff(&index);
		if(max_pen_acc < 0) {
			driveToZero(index);
		} else {
			// already the LCP conditions are satisfied
			break;
		}
	}
}

// drive-to-zero
void LCPSolver::driveToZero(uint ind_d, const double tol) {
	while (_a[ind_d] < -tol) {
		fdirection(ind_d); // this sets _delta_f
		// compute _delta_a
		_delta_a.head(_size) = _A.block(0,0,_size,_size)*_delta_f.head(_size);

		// compute max step
		double max_step = -1.0;
		int index_critical = -1;

		if (_delta_a[ind_d] > 0) {
			max_step = -_a[ind_d]/_delta_a[ind_d];
			index_critical = ind_d;
		}
		// loop over known contacts
		double s;
		for(uint i=0; i < _count_C; ++i) {
			s = -_f[i]/_delta_f[i];
			if (s > 0 && (max_step < 0 || s < max_step)) {
				max_step = s;
				index_critical = i;
			}
		}
		// loop over known non-contacts
		for(uint i=_count_C; i < _count_C+_count_NC; ++i) {
			s = -_a[i]/_delta_a[i];
			if (s > 0 && (max_step < 0 || s < max_step)) {
				max_step = s;
				index_critical = i;
			}
		}
		// TODO: is there really no way to vectorize the above?

		// compute a and f
		_a.head(_size) += max_step * _delta_a.head(_size);
		_f.head(_size) += max_step * _delta_f.head(_size);

		// update indices
		if (index_critical == ind_d) {
			addToC(ind_d);
		} else if (index_critical < _count_C) {
			// move index from C to NC
			moveCtoNC(index_critical);
		} else /* index must be in NC */ {
			// move to C
			moveNCtoC(index_critical);
		}
	}
}

// fdirection
void LCPSolver::fdirection(uint ind_d) {
	_delta_f.head(_size).setZero();
	_delta_f[ind_d] = 1;
	_delta_f.head(_count_C) = _A.block(0,0,_count_C,_count_C).llt().solve(-_A.block(0,ind_d,_count_C,1));
}

// move an index from C to NC
// also reorder the _A, _b, _f and _a matrices to maintain the index order C, NC, d, others
void LCPSolver::moveCtoNC(uint ind) {
	//NOTE: we do not check if index is actually C or if it is already NC since
	//this is an internal fucntion

	// get last index of set C
	uint new_ind = _count_C-1;

	// swap
	swap(ind, new_ind);

	// update _count_C, _count_NC
	--_count_C;
	++_count_NC;
}

// move an index from NC to C
// also reorder the _A, _b, _f and _a matrices to maintain the index order C, NC, d, others
void LCPSolver::moveNCtoC(uint ind) {
	//NOTE: we do not check if index is actually C or if it is already NC since
	//this is an internal fucntion

	// get first index of set NC
	uint new_ind = _count_C;

	// swap
	swap(ind, new_ind);

	// update _count_C, _count_NC
	++_count_C;
	--_count_NC;
}

// add from unknown to C
// also reorder the _A, _b, _f and _a matrices to maintain the index order C, NC, d, others
void LCPSolver::addToC(uint ind) {
	// this operation is performed using two swaps
	// first is to swap the first element of NC to the first element of the unknown set if it exists
	// second is to swap the unknown element sitting at the top of the NC list with the target index
	swap(_count_C, _count_C+_count_NC);
	swap(ind, _count_C);

	// increment _count_C
	++_count_C;
}

// add from unknown to NC
// also reorder the _A, _b, _f and _a matrices to maintain the index order C, NC, d, others
void LCPSolver::addToNC(uint ind) {
	// add to the end of the NC list
	swap(ind, _count_C+_count_NC);

	// increment _count_NC
	++_count_NC;
}

// swap indices in _A, _b, _f and _a. Also update the index order map for bookkeeping
void LCPSolver::swap(uint ind, uint new_ind) {
		if (new_ind != ind) {
		// swap rows and cols of _A
		_A.row(ind).swap(_A.row(new_ind));
		_A.col(ind).swap(_A.col(new_ind));

		// swap rows of _a, _f, _b
		_a.row(ind).swap(_a.row(new_ind));
		_f.row(ind).swap(_f.row(new_ind));
		_b.row(ind).swap(_b.row(new_ind));

		// update index ordering for book keeping
		uint swap_temp = _order_map[ind];
		_order_map[ind] = _order_map[new_ind];
		_order_map[new_ind] = swap_temp;
	}
}
