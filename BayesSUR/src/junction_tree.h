#ifndef JUNCTION_TREE_H
#define JUNCTION_TREE_H

#ifdef CCODE
	#include <iostream>
    #include <armadillo>
#else
    #include <RcppArmadillo.h>
#endif


#include <string>
#include <vector>
#include <deque>
#include <set>
#include <memory>
#include <algorithm> // std::set_difference, std::sort, ...

#include "distr.h"
#include "utils.h"

class JTComponent { 

    public:

        JTComponent( );
        JTComponent( const std::vector<unsigned int>& );
        JTComponent( const std::vector<unsigned int>& , const std::vector<unsigned int>& );
        JTComponent( const std::vector<unsigned int>& , const std::vector<unsigned int>& ,
                          const std::vector<std::shared_ptr<JTComponent>>& , const std::shared_ptr<JTComponent>& );
        JTComponent( const JTComponent& );

        JTComponent copyNode( ); // returns a COPY of this node

        std::vector<unsigned int> getNodes() const;
        std::vector<unsigned int> getSeparator() const;
        std::shared_ptr<JTComponent> getParent() const;
        std::vector<std::shared_ptr<JTComponent>> getChildrens() const;

        void add1Node( const unsigned int );
        void addNodes( const std::vector<unsigned int>& );

        void clearSeparator();
        void add1Separator( const unsigned int );
        void addSeparators( const std::vector<unsigned int>& );
        
        void add1Children( const std::shared_ptr<JTComponent>& );
        void addChildrens( const std::vector<std::shared_ptr<JTComponent>>& );
        
        void setNodes( const std::vector<unsigned int>& );
        void setSeparator( const std::vector<unsigned int>& );
        void setChildrens( const std::vector<std::shared_ptr<JTComponent>>& );
        void setParent( const std::shared_ptr<JTComponent>& );

        void print() const;

    private:

        std::vector<unsigned int> nodes;
        std::vector<unsigned int> separator;
        std::weak_ptr<JTComponent> parent;
        std::vector<std::shared_ptr<JTComponent>> childrens;

};

class JunctionTree {

    public:

        JunctionTree() = default;
        JunctionTree( const unsigned int , const std::string type="" );
        JunctionTree( const unsigned int , const std::deque<std::shared_ptr<JTComponent>>& );

        std::deque<std::shared_ptr<JTComponent>> getPCS() const;
        std::vector<unsigned int> getPEO() const;
        arma::sp_umat getAdjMat() const;
        unsigned int getDimension() const;

        void cloneRoot( std::shared_ptr<JTComponent>& , 
                std::shared_ptr<JTComponent>& );
        void cloneComponent( std::shared_ptr<JTComponent>& , 
                std::shared_ptr<JTComponent>& , 
                std::shared_ptr<JTComponent>& );

        void buildNewPCS( std::deque<std::shared_ptr<JTComponent>>& , unsigned int& );
        void updatePEO();
        void updateAdjMat();

        void copyJT( JunctionTree& ); 
        
        std::pair<bool,double> propose_single_edge_update( );
        std::pair<bool,double> propose_single_edge_update( arma::uvec& );

        std::pair<bool,double> propose_multiple_edge_update( );

        void swapParentChild( std::shared_ptr<JTComponent>& parent , std::shared_ptr<JTComponent>& child );
        void reRoot();
        bool isChild( std::shared_ptr<JTComponent>& , std::shared_ptr<JTComponent>& );
        void randomJTPermutation();

        void print() const;
        // for usability reasons I will have these as public
        // (i.e. I want to use directly the list::insert and other methods)
        std::deque<std::shared_ptr<JTComponent>> perfectCliqueSequence;
        std::vector<unsigned int> perfectEliminationOrder;
        arma::sp_umat adjacencyMatrix;
        unsigned int n;

};

#endif


// note that I treat an empty parent (nullptr, for the root)
// differently than an empty list of childern (which is a vector of length 0, without ANY children ptr initialised)
