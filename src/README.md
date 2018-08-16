# README

Test Repo to keep track of the Bayesian SEM package.

SUR subdirectory only as of now, while we solve the "smaller" problem of one regression node -- i.e. no graph structure for the moment, "just" classic and general VS on a Normal Mv Regression with correlated outcomes in an efficient conditional way

### What is this repository for?

* Bayesian Mv Regression with Variable Selection in the most general way possible.
  See model description for a list of acronyms used, which implies different VS procedure and different correlation structures in the outcomes.

* Novel proposal mechanism / transition kernel inspired by Bayesian Bandit to sample from binary vectors in high dimensions.


### How do I get set up?

* make clean && make -j X
* then see usage.txt for the various algorithm and/or refer to the test_[\*].cpp for a complete list of arguments
