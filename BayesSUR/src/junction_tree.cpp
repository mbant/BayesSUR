#include "junction_tree.h"

#ifndef CCODE
  using Rcpp::Rcout;
#else
  #define Rcout std::cout
#endif

/*
I here decide that the Clique Sequence follows the JT in a Depth First manner, 
so we follow each branch till the end, then come back andd follow another branch and so on...
*/

/*
adjacencyMatrix is symmetric without diagonal
*/


/*
update PEO and update AdjMat both are global functions now, so wildly inefficient I believe
as we migt be updating them only locally. Still...
*/

/*
Nodes and Separator sets are SORTED, mainly to help STL algorithms to serach for inclusions and stuff..
*/


JTComponent::JTComponent( )
{
    nodes = std::vector<unsigned int>(0);
    separator = std::vector<unsigned int>(0);
    // parent = nullptr;
    childrens = std::vector<std::shared_ptr<JTComponent>>(0);
}

JTComponent::JTComponent( const std::vector<unsigned int>& nodes_)
{
    this->setNodes(nodes_);
    separator = std::vector<unsigned int>(0);
    // parent = nullptr;
    childrens = std::vector<std::shared_ptr<JTComponent>>(0);
}

JTComponent::JTComponent( const std::vector<unsigned int>& nodes_ , const std::vector<unsigned int>& separator_)
{
    this->setNodes(nodes_);
    this->setSeparator(separator_);
    // parent = nullptr;
    childrens = std::vector<std::shared_ptr<JTComponent>>(0);
}

JTComponent::JTComponent( const std::vector<unsigned int>& nodes_ , const std::vector<unsigned int>& separator_ ,
                          const std::vector<std::shared_ptr<JTComponent>>& childrens_ , const std::shared_ptr<JTComponent>& parent_)
{
    this->setNodes(nodes_);
    this->setSeparator(separator_);
    this->setParent(parent_);
    this->setChildrens(childrens_);
}

JTComponent::JTComponent( const JTComponent& otherJTComponent )
{
    nodes = otherJTComponent.getNodes();
    separator = otherJTComponent.getSeparator();
    parent = otherJTComponent.getParent();
    childrens = otherJTComponent.getChildrens();
}

std::vector<unsigned int> JTComponent::getNodes() const
{
    return nodes;
}
std::vector<unsigned int> JTComponent::getSeparator() const
{
    return separator;
}

std::shared_ptr<JTComponent> JTComponent::getParent() const
{
    return parent.lock();
}

std::vector<std::shared_ptr<JTComponent>> JTComponent::getChildrens() const
{
    return childrens;
}

void JTComponent::add1Node( const unsigned int node_)
{
    if(std::find(nodes.begin(), nodes.end(), node_) == nodes.end())
    {
        nodes.push_back(node_);
        std::sort(nodes.begin(), nodes.end());
    }
}

void JTComponent::addNodes( const std::vector<unsigned int>& nodes_)
{
    for( auto i : nodes_ )
    {
        if(std::find(nodes.begin(), nodes.end(), i) == nodes.end()) // if it's not there
            nodes.push_back(i);
    }
    std::sort(nodes.begin(), nodes.end());
}


void JTComponent::add1Separator( const unsigned int sep_)
{
    if(std::find(separator.begin(), separator.end(), sep_) == separator.end())
    {
        separator.push_back(sep_);
        std::sort(separator.begin(), separator.end());
    }
}

void JTComponent::addSeparators( const std::vector<unsigned int>& seps_)
{
    for( auto i : seps_ )
    {
        if(std::find(separator.begin(), separator.end(), i) == separator.end()) // if it's not there
            separator.push_back(i);
    }
    std::sort(separator.begin(), separator.end());
}


void JTComponent::add1Children( const std::shared_ptr<JTComponent>& otherJTComponentPTR)
{
    if(std::find(childrens.begin(), childrens.end(), otherJTComponentPTR) == childrens.end())
        childrens.push_back( otherJTComponentPTR );
    // else do nothing as it is a duplicate
}

void JTComponent::addChildrens( const std::vector<std::shared_ptr<JTComponent>>& otherJTComponentPTRs )
{
    for( auto i : otherJTComponentPTRs )
    {
        if(std::find(childrens.begin(), childrens.end(), i) == childrens.end())
            childrens.push_back( i );
        // else do nothing as it is a duplicate
    }
}


void JTComponent::setNodes( const std::vector<unsigned int>& nodes_)
{
    nodes = nodes_ ;
    // remove duplicates from nodes
    nodes.erase( std::unique( nodes.begin(), nodes.end() ), nodes.end() );
    std::sort(nodes.begin(), nodes.end());
}

void JTComponent::setSeparator( const std::vector<unsigned int>& sep_)
{
    separator = sep_;
    // remove duplicates from seo
    separator.erase( std::unique( separator.begin(), separator.end() ), separator.end() );
    std::sort(separator.begin(), separator.end());
}
void JTComponent::setChildrens( const  std::vector<std::shared_ptr<JTComponent>>& c)
{
    // remove duplicates from c
    childrens = c;
    childrens.erase( std::unique( childrens.begin(), childrens.end() ), childrens.end() );
}

void JTComponent::clearSeparator()
{
    separator.clear();
}

void JTComponent::setParent( const std::shared_ptr<JTComponent>& otherJTComponentPTR )
{
    parent = otherJTComponentPTR;
}

void JTComponent::print() const
{
    Rcout << "JT Component @ address " << this <<" is made of Nodes :";
    for( auto i : nodes )
        Rcout << " " << i;
    Rcout << '\n';

    Rcout << "  with Separator :";
    for( auto i : separator )
        Rcout << " " << i;
    Rcout << '\n';

    Rcout << "  Its Parent is @ " << parent.lock() << " and its Children are @:";
    for( auto i : childrens )
        Rcout << " " << i;
    Rcout << '\n' << '\n';
}

// ##########################################################################################
// ##########################################################################################
// ##########################################################################################

// create a JT with n disconnected components
JunctionTree::JunctionTree( const unsigned int n_ , const std::string type ) // type has a default of ""
{

    n = n_;

    if( type == "" || type == "empty" )
    {
        perfectEliminationOrder = std::vector<unsigned int>(n);
        perfectCliqueSequence = std::deque<std::shared_ptr<JTComponent>>(); // init empty list


        std::vector<unsigned int> node(1); 


        std::deque<std::shared_ptr<JTComponent>>::iterator pcs_it = perfectCliqueSequence.end();

        node[0] = 0;
        perfectCliqueSequence.insert( pcs_it, std::make_shared<JTComponent>( node ));
        // the above returns a shared ptr to a new component,
        // the component has one node, empty separator and no parent or child
        // and it adds it to the pcs_it position in the list

        pcs_it = perfectCliqueSequence.end(); //update iterator
        --pcs_it;

        perfectEliminationOrder[0] = 0;

        std::shared_ptr<JTComponent> tmpPointerToComponent;

        
        for( unsigned int i=1; i<n; ++i)
        {
            node[0] = i;
            perfectEliminationOrder[i] = i;

            // create a new component
            tmpPointerToComponent = std::make_shared<JTComponent>( node );
            // insert the new compoent as child in the previous component
            (*pcs_it) -> add1Children( tmpPointerToComponent );
            // add the previous component as parent of the new component
            tmpPointerToComponent -> setParent( *pcs_it );
            // increment iterator
            ++pcs_it;
            // add new component to the list
            perfectCliqueSequence.insert( pcs_it, tmpPointerToComponent );

            // update iterator since the previous insert deprecated it
            pcs_it = perfectCliqueSequence.end();
            --pcs_it;
        }

        adjacencyMatrix.zeros(n,n); // the matrix is empty, no edges

    }else{
        perfectEliminationOrder = std::vector<unsigned int>(n);
        perfectCliqueSequence = std::deque<std::shared_ptr<JTComponent>>();

        for( unsigned int i=0; i<n; ++i)
            perfectEliminationOrder[i] = i;

        perfectCliqueSequence.push_back( std::make_shared<JTComponent>(perfectEliminationOrder) );
        
        arma::umat tmp = arma::ones<arma::umat>(n,n);
        adjacencyMatrix = tmp - arma::eye<arma::umat>(n,n); // the matrix is full of edges

    }

}


JunctionTree::JunctionTree( const unsigned int n_, const std::deque<std::shared_ptr<JTComponent>>& PCS_)
{
    n = n_;
    perfectCliqueSequence = PCS_;
    updatePEO();
    updateAdjMat();
}

std::deque<std::shared_ptr<JTComponent>> JunctionTree::getPCS() const
{
    return perfectCliqueSequence;
}

std::vector<unsigned int> JunctionTree::getPEO() const
{
    return perfectEliminationOrder;
}

arma::sp_umat JunctionTree::getAdjMat() const
{
    return adjacencyMatrix;
}


unsigned int JunctionTree::getDimension() const
{
    return n;
}

void JunctionTree::print() const
{
    Rcout << '\n' << " ---------------------------------- " << '\n';
    for( auto i : perfectCliqueSequence )
        i->print();
    Rcout << " ---------------------------------- " << '\n' <<
        "The PEO for this JT is :" << '\n';

    for(auto i : perfectEliminationOrder )
        Rcout << i << " ";
    Rcout << '\n' << " ---------------------------------- " << '\n';

    arma::umat tmp(adjacencyMatrix);
    Rcout << "Graph's Adjacency Matrix: " << tmp << '\n' << '\n';
}

void JunctionTree::cloneRoot( std::shared_ptr<JTComponent>& newComponent , 
                                    std::shared_ptr<JTComponent>& oldComponent )
{
    newComponent->setNodes( oldComponent->getNodes() );
    newComponent->setSeparator( oldComponent->getSeparator() );
    newComponent->setParent( nullptr );

    unsigned int nChildrens = oldComponent->getChildrens().size();

    std::vector<std::shared_ptr<JTComponent>> newChildrens( nChildrens );
    std::vector<std::shared_ptr<JTComponent>> oldChildrens = oldComponent->getChildrens();

    for( unsigned int i=0; i<nChildrens; ++i )
    {
        newChildrens[i] = std::make_shared<JTComponent>();
        cloneComponent( newChildrens[i], oldChildrens[i], newComponent);
    }

    newComponent->setChildrens( newChildrens );
}


void JunctionTree::cloneComponent( std::shared_ptr<JTComponent>& newComponent , 
                                    std::shared_ptr<JTComponent>& oldComponent , 
                                    std::shared_ptr<JTComponent>& newParent )
{
    newComponent->setNodes( oldComponent->getNodes() );
    newComponent->setSeparator( oldComponent->getSeparator() );
    newComponent->setParent( newParent );

    unsigned int nChildrens = oldComponent->getChildrens().size();

    std::vector<std::shared_ptr<JTComponent>> newChildrens( nChildrens );
    std::vector<std::shared_ptr<JTComponent>> oldChildrens = oldComponent->getChildrens();

    for( unsigned int i=0; i<nChildrens; ++i )
    {
        newChildrens[i] = std::make_shared<JTComponent>();
        cloneComponent( newChildrens[i], oldChildrens[i], newComponent);
    }

    newComponent->setChildrens( newChildrens );
}


/* 
This assumes that all the elements before pos are correctly initialised
*/
void JunctionTree::buildNewPCS( std::deque<std::shared_ptr<JTComponent>>& newPCS, unsigned int& pos)
{

    unsigned int j=0; // index for the childrens inside each 'root' node
    std::vector<std::shared_ptr<JTComponent>> rootChildrens = newPCS[pos]->getChildrens();
    std::vector<std::shared_ptr<JTComponent>> childChildrens;

    while( ( rootChildrens.size() - j ) > 0 )
    {
        newPCS.insert( newPCS.begin() + (++pos) , rootChildrens[j] );  // I love how increments in c++ make things less readable..
        childChildrens = rootChildrens[j++]->getChildrens();

        if( childChildrens.size() > 0 )
        {
            buildNewPCS(newPCS, pos);
        }
    }
}

void JunctionTree::copyJT( JunctionTree& newJT ) // returns a copy of this, copying each component into newJT
{
  
    std::deque<std::shared_ptr<JTComponent>> newPCS;

    newPCS.insert( newPCS.end() , std::make_shared<JTComponent>() );

    cloneRoot(newPCS[0],perfectCliqueSequence[0]); //this clones the root component and ALL its childrens

    // build up again the newPCS (the O-th element is already set)
    unsigned int pcsPosition=0; // index for the PCS
    buildNewPCS( newPCS, pcsPosition);

    newJT = JunctionTree( this->n, newPCS );
}


void JunctionTree::updatePEO( )
{
    std::vector<unsigned int> tmpResidual;

    perfectEliminationOrder.clear();
    perfectEliminationOrder.reserve(n);

    std::vector<unsigned int> tmpNodes,tmpSeparator;

    for( auto c : perfectCliqueSequence )
    {
        tmpResidual.clear();
        
        tmpNodes = c->getNodes();
        tmpSeparator = c->getSeparator();

        std::set_difference(tmpNodes.begin(), tmpNodes.end(), 
                            tmpSeparator.begin(), tmpSeparator.end(), 
                        std::inserter(tmpResidual, tmpResidual.begin()));

        perfectEliminationOrder.insert(perfectEliminationOrder.end(),
                                    tmpResidual.begin(), tmpResidual.end()); 
    }
}


void JunctionTree::updateAdjMat( )
{
    adjacencyMatrix.zeros(n,n); //this maintains dimensions but reset it

    std::vector<unsigned int> componentNodes;

    for(auto i : perfectCliqueSequence)
    {
        componentNodes = i->getNodes();
        if( componentNodes.size() > 1 )
        {

            for( unsigned int j=0, nNodes=componentNodes.size() ; j<(nNodes-1) ; ++j )  // j < k -- nodes are ordered
            {
                for( unsigned int k=j+1 ; k<nNodes ; ++k )                          // so element (k,j) is in the lower tri
                {
                    if( !adjacencyMatrix(componentNodes[k],componentNodes[j]) )
                    {
                        adjacencyMatrix(componentNodes[k],componentNodes[j]) = 1;
                        adjacencyMatrix(componentNodes[j],componentNodes[k]) = 1;
                    }
                }
            }

        }
    }

    // adjacencyMatrix = adjacencyMatrix + adjacencyMatrix.t(); // to fill the upper tri part.
}


/*
JT Update proposal, from Green and Thomas 2013
The idea is to call this from a copied JT, as it will modify the JT structure (is this TOO costly?)
*/
std::pair<bool,double> JunctionTree::propose_single_edge_update( arma::uvec& update_idx )
{
    // We need to select randomly one separator -- in our structure is any component beside the first
    unsigned int numComponents = perfectCliqueSequence.size();
    unsigned int randomSep, randomComp;

    arma::uvec randomIndexes;

    std::shared_ptr<JTComponent> Cx, Cy, C, ClX, ClY, cLeft, cRight;
    std::vector<unsigned int> setCx, setCy, setCxlS, setCylS, setS, 
        setC, setNeighbour, xUS, yUS, setCLeft, setCRight;

    unsigned int x,y;

    std::vector<std::shared_ptr<JTComponent>> newChildrens;
    std::vector<unsigned int> newNodes;
    std::vector<unsigned int> newSeparator;
    std::shared_ptr<JTComponent> newParent;

    std::vector<std::shared_ptr<JTComponent>> neighbours, N, Nx, Ny, possibleComponents;
    bool definedCx=false, definedCy=false;
    std::deque<std::shared_ptr<JTComponent>> newPCS;
    unsigned int pos = 0;
    unsigned int countN = 0;

    double logP = 0;

    if( randU01() < 0.5 )
    {
        // propose addition
        // if there's one component only, return false and no update happened
        if( numComponents > 1 )
        {
            // get a Random separator
            randomSep = randIntUniform(1,numComponents-1);

            // populate Cx and Cy
            Cx = perfectCliqueSequence[randomSep]->getParent();
            Cy = perfectCliqueSequence[randomSep];

            setCx = Cx->getNodes();
            setCy = Cy->getNodes();
            setS = perfectCliqueSequence[randomSep]->getSeparator();
            
            // Populate Cx\S and Cy\S
            std::set_difference(setCx.begin(), setCx.end(), setS.begin(), setS.end(), 
                            std::inserter(setCxlS, setCxlS.begin()));

            std::set_difference(setCy.begin(), setCy.end(), setS.begin(), setS.end(), 
                            std::inserter(setCylS, setCylS.begin()));

            // now choose x and y from Cx\S and Cy\S
            x = setCxlS[ randIntUniform(0,setCxlS.size()-1) ];
            y = setCylS[ randIntUniform(0,setCylS.size()-1) ];

            // check whether or not Cx and Cy are superset of x U S and y U S and act accordingly
            // note there's no way for it to be the other way around by construction
            if( setCxlS.size() == 1 && setCylS.size() == 1 ) // a) this means that both Cx and Cy are exactly just x and S (and y and S) 
            {
                
                // Type a)

                // compute new children set
                newChildrens = Cy->getChildrens();
                for( auto c : Cx->getChildrens() )
                {
                    if( c != Cy )
                        newChildrens.push_back(c);
                }


                newNodes = setS; 
                newNodes.push_back(x);
                newNodes.push_back(y);

                // newSeparator = Cx->getSeparator();
                // newParent = Cx->getParent();

                // Modify Cx to C* (the new component we need)
                Cx->setChildrens(newChildrens);
                Cx->setNodes(newNodes);
                // Cx->setSeparator(newSeparator);
                // Cx->setParent(newParent);

                // make the new childrens point to C* (rather than Cy)
                for( auto c : newChildrens )
                {
                    if( c->getParent() != Cx )
                        c->setParent(Cx);
                }

                // Erase Cy
                perfectCliqueSequence.erase(perfectCliqueSequence.begin() + randomSep);

                // **** logP addition
                neighbours.clear();
                if( Cx->getParent() )
                    neighbours.push_back( Cx->getParent() );  // push parent only if it exists

                for( auto i : Cx->getChildrens() )
                    neighbours.push_back(i);

                for( auto i : neighbours )
                {
                    setNeighbour = i->getNodes();
                    if(std::find(setNeighbour.begin(), setNeighbour.end(), x) == setNeighbour.end()) // if x is NOT in there
                    {
                        if(std::find(setNeighbour.begin(), setNeighbour.end(), y) == setNeighbour.end()) // if y is NOT in there either
                        {
                            ++countN;
                        
                        }

                    }
                }

                logP -= countN * log(2.); // backward probability addition
                logP -= log((double)(numComponents-1.)) - log(2.) + ( log( (double)(Cx->getNodes().size()) ) + log( (double)(Cx->getNodes().size()-1.) ) ); // backward probability (-1 because we reduce the # component by 1)

            }else if( setCxlS.size() > 1 && setCylS.size() == 1 ) // b)
            {

                // Type b)

                Cy->add1Node(x);
                Cy->add1Separator(x);

                logP -= log((double)numComponents) - log(2.) + ( log( (double)(Cx->getNodes().size()) ) + log( (double)(Cx->getNodes().size()-1.) ) ); // backward probability

            }else if( setCxlS.size() == 1 && setCylS.size() > 1 ) //c)
            {
                // Type c)

                Cx->add1Node(y);
                Cy->add1Separator(y); // remember the separator is always the one from the Cy obj

                logP -= log((double)numComponents) - log(2.) + ( log( (double)(Cx->getNodes().size()) ) + log( (double)(Cx->getNodes().size()-1.) ) ); // backward probability

            }else if( setCxlS.size() > 1 && setCylS.size() > 1 ) // d) Cx and Cy contain more than just x and S (and y and S) 
            {
                // Type d)

                // Create C* (the new Component)
                newChildrens = std::vector<std::shared_ptr<JTComponent>>(1);
                newChildrens[0] = Cy;

                newNodes = setS; 
                newNodes.push_back(x);
                newNodes.push_back(y);

                newSeparator = setS;
                newSeparator.push_back(x);
                
                newParent = Cx;

                // insert it in the PCS
                perfectCliqueSequence.insert( perfectCliqueSequence.begin() + randomSep ,
                                std::make_shared<JTComponent>(newNodes,newSeparator,newChildrens,newParent) );

                // Modify Cx to point to C*
                newChildrens = Cx->getChildrens();
                for (auto itC = newChildrens.begin(); itC != newChildrens.end();  ++itC )
                {
                    if( *(itC) == Cy )
                        *(itC) = perfectCliqueSequence[randomSep];
                }
                Cx->setChildrens(newChildrens);

                // Modify Cy to have C* as parent and modify its separator
                Cy->setParent(perfectCliqueSequence[randomSep]);

                Cy->add1Separator(y);         


                logP -= log((double)(numComponents+1.)) - log(2.) +
                    ( log( (double)(perfectCliqueSequence[randomSep]->getNodes().size()) ) +
                        log( (double)(perfectCliqueSequence[randomSep]->getNodes().size()-1.) ) ); // backward probability (+1 because we inserta new component here)
                
            }

            logP += log((double)(numComponents-1.)) + log((double)(setCxlS.size())) + log((double)(setCylS.size())) ; // forward probability

        }else //if only one there's no option for edge addition
        {
            return std::make_pair(false,0.0); // and the graph is unchanged
        }
    
    }else         // propose deletion
    {
        // get a random component
        randomComp = randIntUniform(0,numComponents-1);
        C = perfectCliqueSequence[randomComp];

        setC = C->getNodes();

        if( setC.size() > 1 )
        {
            // deletion has its forward logP computed immediately as setC is likely to be modified after
            logP += log((double)numComponents) - log(2.) + ( log( (double)(setC.size()) ) + log( (double)(setC.size()-1.) ) ); // forward probability

            // partition C into 3 sets, x y and S (S might be empty)"
            randomIndexes = Distributions::randWeightedIndexSampleWithoutReplacement(setC.size(),2);

            x = setC[randomIndexes(0)];
            y = setC[randomIndexes(1)];

            std::copy(setC.begin(), setC.end(),
                std::back_inserter(setS));
            setS.erase(std::remove(setS.begin(), setS.end(), x), setS.end());
            setS.erase(std::remove(setS.begin(), setS.end(), y), setS.end());

            // now scan the neighbours of C and construct Nx, Ny and N"
            neighbours.clear();
            if( C->getParent() )
                neighbours.push_back( C->getParent() );  // push parent only if it exists

            for( auto i : C->getChildrens() )
                neighbours.push_back(i);

            for( auto i : neighbours )
            {
                setNeighbour = i->getNodes();
                if(std::find(setNeighbour.begin(), setNeighbour.end(), x) != setNeighbour.end()) // if x is in there
                {
                    if(std::find(setNeighbour.begin(), setNeighbour.end(), y) != setNeighbour.end()) // if y is in there as well
                    {
                        return std::make_pair(false,0.0); // exit, deletion is not possible
                    
                    }else{ // i contains only x

                        Nx.push_back(i);
                    }

                }else if(std::find(setNeighbour.begin(), setNeighbour.end(), y) != setNeighbour.end()) // if y is in there (but not x)
                {
                    Ny.push_back(i);
                }else{ // nor x nor y are in there
                    N.push_back(i);
                }
            }

            // search Nx for x U S and select Cx if possible
            possibleComponents.clear(); 

            std::copy(setS.begin(), setS.end(),
                std::back_inserter(xUS));
            xUS.push_back(x);
            std::sort(xUS.begin(),xUS.end());

            for( auto i : Nx )
            {
                setNeighbour = i->getNodes();
                if( std::includes(setNeighbour.begin(), setNeighbour.end(), xUS.begin(), xUS.end()) ) // note that includes only works on sorted ranges
                    possibleComponents.push_back(i);
            }
            if( possibleComponents.size() > 0 )
            {
                Cx = possibleComponents[ randIntUniform(0,possibleComponents.size()-1) ];
                definedCx = true;
            }

            // search Ny for y U S and select Cy if possible
            possibleComponents.clear();

            std::copy(setS.begin(), setS.end(),
                std::back_inserter(yUS));
            yUS.push_back(y);
            std::sort(yUS.begin(),yUS.end());

            for( auto i : Ny )
            {
                setNeighbour = i->getNodes();
                if( std::includes(setNeighbour.begin(), setNeighbour.end(), yUS.begin(), yUS.end()) ) // note that includes only works on sorted ranges
                    possibleComponents.push_back(i);
            }
            
            if( possibleComponents.size() > 0 )
            {
                Cy = possibleComponents[ randIntUniform(0,possibleComponents.size()-1) ];
                definedCy = true;
            }

            // Now onto the 4 cases

            if( !definedCx && !definedCy ) // if both are undefined
            {
                // Type a)

                // Create two new JTComponents, separated by S and made up by xUS and yUS
                ClX = std::make_shared<JTComponent>( setS );
                ClX->add1Node(y);

                ClY = std::make_shared<JTComponent>( setS );
                ClY->add1Node(x);

                // decide which (Clx or Cly) is the left and which the right clique
                if( std::find(Nx.begin(), Nx.end(), C->getParent() ) != Nx.end() ) // the parent was in Nx, hence I need ClY on the left
                {
                    cLeft = ClY;  // here I'm copying the pointer, so I can use them interchangeably
                    cRight = ClX;

                }else if( std::find(Ny.begin(), Ny.end(), C->getParent() ) != Ny.end() ) // the parent was in Ny, hence I need ClX on the left
                {
                    cLeft = ClX;  
                    cRight = ClY;

                }else{ // either C was the root or the parent is in N, so I can choose randomly

                    if( randU01() < 0.5 )
                    {
                        cLeft = ClY;  
                        cRight = ClX;
                    }else{
                        cLeft = ClX;  
                        cRight = ClY;                        
                    }
                }

                // regardless, add now cLeft to its childrens and remove C
                if( C->getParent() )
                {
                    newChildrens = C->getParent()->getChildrens();
                    newChildrens.erase(std::remove(newChildrens.begin(), newChildrens.end(), C), newChildrens.end());
                    newChildrens.push_back( cLeft );
                    C->getParent()->setChildrens( newChildrens );
                } // else it was the root so no changes

                // Add the original parent to cLeft                
                cLeft->setParent( C->getParent() );
                // and its separator, which might be empty
                cLeft->setSeparator( C->getSeparator() );

                // the parent of the right clique is then the left clique 
                // (and the right is a child for the left)
                // their separator is S
                cLeft->add1Children( cRight ); // note that this is just ONE children and there were none before
                cRight->setParent( cLeft );
                cRight->setSeparator( setS );    

                // Now connect all the neighbours to either cLeft or cRight (except the original parent which is already connected)
                for( auto i : Nx )
                {
                    if( i != C->getParent() )
                    {
                        ClY->add1Children( i );
                        i->setParent( ClY );
                    }
                }

                for( auto i : Ny )
                {
                    if( i != C->getParent() )
                    {
                        ClX->add1Children( i );
                        i->setParent( ClX );
                    }
                }

                for( auto i : N )
                {
                    if( i != C->getParent() )
                    {
                        if( randU01() < 0.5 )
                        {
                            ClX->add1Children( i );
                            i->setParent( ClX );

                        }else{
                            ClY->add1Children( i );
                            i->setParent( ClY );
                        }
                    }
                }

                // Finally erase C from PCS .. (this actually gets overwritten so no need..)
                // perfectCliqueSequence.erase(perfectCliqueSequence.begin() + randomComp);

                // .. and Re-build the PCS after this
                newPCS.clear();
                if( !( C->getParent() ) )
                {
                    // cLeft is the new root
                    newPCS.insert( newPCS.end() , cLeft );
                }else{
                    // the old root is still the root, but I've modified its childrens
                    newPCS.insert( newPCS.end() , perfectCliqueSequence[0] );
                }
                pos = 0;
                buildNewPCS( newPCS, pos );
                perfectCliqueSequence = newPCS; // substitute to the current one
                // PEO updated below

                logP += N.size() * log(2.); // forward probability addition
                logP -= log((double)numComponents) + log((double)(ClX->getNodes().size() - setS.size())) + log((double)(ClY->getNodes().size() - setS.size())) ; //backward probability (we added one component here so numComponents rather than numComponents-1)

            }else if( definedCx && !definedCy ) // if only Cx is defined
            {
                // Type b)

                // Separation is possible only if Nx contains only Cx
                if( Nx.size() == 1 )
                {
                    if( Nx[0] == Cx )
                    {
                        // remove x from C
                        setC.erase(std::remove(setC.begin(), setC.end(), x), setC.end());
                        C->setNodes( setC );
                        // find the separator that connects C to Cx
                        if( Cx == C->getParent() )
                        {
                            setS = C->getSeparator();
                            setS.erase(std::remove(setS.begin(), setS.end(), x), setS.end());
                            C->setSeparator( setS );
                        }else{ // it's a child
                            setS = Cx->getSeparator();
                            setS.erase(std::remove(setS.begin(), setS.end(), x), setS.end());
                            Cx->setSeparator( setS );
                        }
                    }else{
                        return std::make_pair(false,0.0); // and the graph is unchanged
                    }
                }else{
                    return std::make_pair(false,0.0); // and the graph is unchanged
                }

                logP -= log((double)(numComponents-1.)) + log((double)(Cx->getNodes().size() - setS.size())) + log((double)(setC.size() - setS.size())) ; //backward probability

            }else if( !definedCx && definedCy ) // if only Cy is defined
            {
                // Type c)

                // Separation is possible only if Ny contains only Cy
                if( Ny.size() == 1 )
                {
                    if( Ny[0] == Cy )
                    {
                        // remove y from C
                        setC.erase(std::remove(setC.begin(), setC.end(), y), setC.end());
                        C->setNodes( setC );
                        // find the separator that connects C to Cy
                        if( Cy == C->getParent() )
                        {
                            setS = C->getSeparator();
                            setS.erase(std::remove(setS.begin(), setS.end(), y), setS.end());
                            C->setSeparator( setS );
                        }else{ // it's a child
                            setS = Cy->getSeparator();
                            setS.erase(std::remove(setS.begin(), setS.end(), y), setS.end());
                            Cy->setSeparator( setS );
                        }
                    }else{
                        return std::make_pair(false,0.0); // and the graph is unchanged
                    }
                }else{
                    return std::make_pair(false,0.0); // and the graph is unchanged
                }
                
                logP -= log((double)(numComponents-1.)) + log((double)(Cy->getNodes().size() - setS.size()) ) + log((double)(setC.size() - setS.size())) ; //backward probability
                
            }else{ //if both are defined
                // Type d)

                // Separation is possible only if N* contains only C* [*=x,y] and N is empty
                if( Nx.size() == 1 && Ny.size() == 1 && N.size() == 0 )
                {
                    if( Nx[0]==Cx && Ny[0]==Cy)
                    {
                        // Find who's first between Cx and Cy
                        // adjust their links
                        // and change their separator to be a new S
                        // the separator for cLeft now is either its old separator (if it was the parent)
                        // or it needs to be set to the separator of C if cLeft was its child
                        if( Cx == C->getParent() )
                        {
                            cLeft = Cx;
                            cRight = Cy;

                        }else if( Cy == C->getParent() ) // if Cy is the parent
                        {
                            cLeft = Cy;
                            cRight = Cx;

                        }else if( !( C->getParent()) ) // C was their parent, so choose randomly
                        {                               // note as well that in this case (becase N is empty), C was the root of JT
                            if( randU01() < 0.5 )
                            {
                                cLeft = Cx;
                                cRight = Cy;
                            }else{
                                cLeft = Cy;
                                cRight = Cx;
                            }

                            // because C was the root of the tree, the clique chosen to replace him at the top get its separator emptied
                            newSeparator.clear();
                            cLeft->setSeparator( newSeparator );

                            // and its parent's emptied as well as it becomes the new root
                            cLeft->setParent( nullptr );

                            // DO WE NEED TO CHANGE THE PCS in this case? (it's safer to do so, but... do we need?)
                            // I'll do it because our code relies on the fact that perfectCliqueSequence[0] is the root

                        }

                        // Now fix links
                        // find C in cLeft's childs and substitute (or add) cRight
                        newChildrens = cLeft->getChildrens();
                        newChildrens.erase(std::remove(newChildrens.begin(), newChildrens.end(), C), newChildrens.end());
                        newChildrens.push_back( cRight );
                        cLeft->setChildrens( newChildrens );

                        // set cLeft as cRight's parent (cRight HAS to come from C's childs)
                        cRight->setParent( cLeft );

                        // set the separator between them as the intersection between them
                        newSeparator.clear();
                        setCLeft = cLeft->getNodes();
                        setCRight = cRight->getNodes();

                        std::set_intersection(setCLeft.begin(), setCLeft.end(),
                                                setCRight.begin(), setCRight.end(),
                                            std::back_inserter(newSeparator));

                        cRight->setSeparator( newSeparator );
                        
                        // Finally erase C
                        perfectCliqueSequence.erase(perfectCliqueSequence.begin() + randomComp);

                        if( !( C->getParent()) ) // if we erased the root node, we created a new one and we need to create a newPCS
                        {
                            newPCS.clear();
                            // cLeft is the new root
                            newPCS.insert( newPCS.end() , cLeft );
                            pos = 0;
                            buildNewPCS( newPCS, pos );
                            perfectCliqueSequence = newPCS; // substitute to the current one
                            // PEO updated below
                        }

                    }else{
                        return std::make_pair(false,0.0); // and the graph is unchanged
                    }

                }else{
                    return std::make_pair(false,0.0); // and the graph is unchanged
                }
                
                logP -= log((double)(numComponents-2.)) + log((double)(setCLeft.size()-newSeparator.size())) + log((double)(setCRight.size()-newSeparator.size())) ; //backward probability (-2 here cause we further deleted one component)
        
            }

        }else //if only one there's no option for edge deletion
        {
            return std::make_pair(false,0.0); // and the graph is unchanged
        }

    }

    updatePEO(); // this->updatePEO();
    updateAdjMat();

    update_idx.zeros(2); update_idx(0) = x; update_idx(1) = y;
    return std::make_pair(true,logP);

}

/*
JT Update proposal, from Green and Thomas 2013
The idea is to call this from a copied JT, as it will modify the JT structure (is this TOO costly?)
*/
std::pair<bool,double> JunctionTree::propose_single_edge_update( )
{
    // We need to select randomly one separator -- in our structure is any component beside the first
    unsigned int numComponents = perfectCliqueSequence.size();
    unsigned int randomSep, randomComp;

    arma::uvec randomIndexes;

    std::shared_ptr<JTComponent> Cx, Cy, C, ClX, ClY, cLeft, cRight;
    std::vector<unsigned int> setCx, setCy, setCxlS, setCylS, setS, 
        setC, setNeighbour, xUS, yUS, setCLeft, setCRight;

    unsigned int x,y;

    std::vector<std::shared_ptr<JTComponent>> newChildrens;
    std::vector<unsigned int> newNodes;
    std::vector<unsigned int> newSeparator;
    std::shared_ptr<JTComponent> newParent;

    std::vector<std::shared_ptr<JTComponent>> neighbours, N, Nx, Ny, possibleComponents;
    bool definedCx=false, definedCy=false;
    std::deque<std::shared_ptr<JTComponent>> newPCS;
    unsigned int pos = 0;
    unsigned int countN = 0;

    double logP = 0;

    if( randU01() < 0.5 )
    {
        // propose addition
        // if there's one component only, return false and no update happened
        if( numComponents > 1 )
        {
            // get a Random separator
            randomSep = randIntUniform(1,numComponents-1);

            // populate Cx and Cy
            Cx = perfectCliqueSequence[randomSep]->getParent();
            Cy = perfectCliqueSequence[randomSep];

            setCx = Cx->getNodes();
            setCy = Cy->getNodes();
            setS = perfectCliqueSequence[randomSep]->getSeparator();
            
            // Populate Cx\S and Cy\S
            std::set_difference(setCx.begin(), setCx.end(), setS.begin(), setS.end(), 
                            std::inserter(setCxlS, setCxlS.begin()));

            std::set_difference(setCy.begin(), setCy.end(), setS.begin(), setS.end(), 
                            std::inserter(setCylS, setCylS.begin()));

            // now choose x and y from Cx\S and Cy\S
            x = setCxlS[ randIntUniform(0,setCxlS.size()-1) ];
            y = setCylS[ randIntUniform(0,setCylS.size()-1) ];

            // check whether or not Cx and Cy are superset of x U S and y U S and act accordingly
            // note there's no way for it to be the other way around by construction
            if( setCxlS.size() == 1 && setCylS.size() == 1 ) // a) this means that both Cx and Cy are exactly just x and S (and y and S) 
            {
                
                // Type a)

                // compute new children set
                newChildrens = Cy->getChildrens();
                for( auto c : Cx->getChildrens() )
                {
                    if( c != Cy )
                        newChildrens.push_back(c);
                }


                newNodes = setS; 
                newNodes.push_back(x);
                newNodes.push_back(y);

                // newSeparator = Cx->getSeparator();
                // newParent = Cx->getParent();

                // Modify Cx to C* (the new component we need)
                Cx->setChildrens(newChildrens);
                Cx->setNodes(newNodes);
                // Cx->setSeparator(newSeparator);
                // Cx->setParent(newParent);

                // make the new childrens point to C* (rather than Cy)
                for( auto c : newChildrens )
                {
                    if( c->getParent() != Cx )
                        c->setParent(Cx);
                }

                // Erase Cy
                perfectCliqueSequence.erase(perfectCliqueSequence.begin() + randomSep);

                // **** logP addition
                neighbours.clear();
                if( Cx->getParent() )
                    neighbours.push_back( Cx->getParent() );  // push parent only if it exists

                for( auto i : Cx->getChildrens() )
                    neighbours.push_back(i);

                for( auto i : neighbours )
                {
                    setNeighbour = i->getNodes();
                    if(std::find(setNeighbour.begin(), setNeighbour.end(), x) == setNeighbour.end()) // if x is NOT in there
                    {
                        if(std::find(setNeighbour.begin(), setNeighbour.end(), y) == setNeighbour.end()) // if y is NOT in there either
                        {
                            ++countN;
                        }

                    }
                }

                logP -= countN * log(2.); // backward probability addition
                logP -= log((double)(numComponents-1.)) - log(2.) + ( log( (double)(Cx->getNodes().size() )) + log( (double)(Cx->getNodes().size()-1. )) ); // backward probability (-1 because we reduce the # component by 1)

           }else if( setCxlS.size() > 1 && setCylS.size() == 1 ) // b)
            {

                // Type b)

                Cy->add1Node(x);
                Cy->add1Separator(x);

                logP -= log((double)numComponents) - log(2.) + ( log( (double)(Cx->getNodes().size() )) + log( (double)(Cx->getNodes().size()-1.) ) ); // backward probability

            }else if( setCxlS.size() == 1 && setCylS.size() > 1 ) //c)
            {
                // Type c)

                Cx->add1Node(y);
                Cy->add1Separator(y); // remember the separator is always the one from the Cy obj

                logP -= log((double)numComponents) - log(2.) + ( log( (double)(Cx->getNodes().size()) ) + log( (double)(Cx->getNodes().size()-1.) ) ); // backward probability
            }else if( setCxlS.size() > 1 && setCylS.size() > 1 ) // d) Cx and Cy contain more than just x and S (and y and S) 
            {
                // Type d)

                // Create C* (the new Component)
                newChildrens = std::vector<std::shared_ptr<JTComponent>>(1);
                newChildrens[0] = Cy;

                newNodes = setS; 
                newNodes.push_back(x);
                newNodes.push_back(y);

                newSeparator = setS;
                newSeparator.push_back(x);
                
                newParent = Cx;

                // insert it in the PCS
                perfectCliqueSequence.insert( perfectCliqueSequence.begin() + randomSep ,
                                std::make_shared<JTComponent>(newNodes,newSeparator,newChildrens,newParent) );

                // Modify Cx to point to C*
                newChildrens = Cx->getChildrens();
                for (auto itC = newChildrens.begin(); itC != newChildrens.end();  ++itC )
                {
                    if( *(itC) == Cy )
                        *(itC) = perfectCliqueSequence[randomSep];
                }
                Cx->setChildrens(newChildrens);

                // Modify Cy to have C* as parent and modify its separator
                Cy->setParent(perfectCliqueSequence[randomSep]);

                Cy->add1Separator(y);         


                logP -= log((double)(numComponents+1.)) - log(2.) +
                    ( log( (double)(perfectCliqueSequence[randomSep]->getNodes().size()) ) +
                        log( (double)(perfectCliqueSequence[randomSep]->getNodes().size()-1.) ) ); // backward probability (+1 because we inserta new component here)
                
            }

            logP += log((double)(numComponents-1.)) + log((double)(setCxlS.size())) + log((double)(setCylS.size())) ; // forward probability

        }else //if only one there's no option for edge addition
        {
            return std::make_pair(false,0.0); // and the graph is unchanged
        }
    
    }else         // propose deletion
    {
        // get a random component
        randomComp = randIntUniform(0,numComponents-1);
        C = perfectCliqueSequence[randomComp];

        setC = C->getNodes();

        if( setC.size() > 1 )
        {
            // deletion has its forward logP computed immediately as setC is likely to be modified after
            logP += log((double)numComponents) - log(2.) + ( log( (double)(setC.size()) ) + log( (double)(setC.size()-1.) ) ); // forward probability

            // partition C into 3 sets, x y and S (S might be empty)"
            randomIndexes = Distributions::randWeightedIndexSampleWithoutReplacement(setC.size(),2);

            x = setC[randomIndexes(0)];
            y = setC[randomIndexes(1)];

            std::copy(setC.begin(), setC.end(),
                std::back_inserter(setS));
            setS.erase(std::remove(setS.begin(), setS.end(), x), setS.end());
            setS.erase(std::remove(setS.begin(), setS.end(), y), setS.end());

            // now scan the neighbours of C and construct Nx, Ny and N"
            neighbours.clear();
            if( C->getParent() )
                neighbours.push_back( C->getParent() );  // push parent only if it exists

            for( auto i : C->getChildrens() )
                neighbours.push_back(i);

            for( auto i : neighbours )
            {
                setNeighbour = i->getNodes();
                if(std::find(setNeighbour.begin(), setNeighbour.end(), x) != setNeighbour.end()) // if x is in there
                {
                    if(std::find(setNeighbour.begin(), setNeighbour.end(), y) != setNeighbour.end()) // if y is in there as well
                    {
                        return std::make_pair(false,0.0); // exit, deletion is not possible
                    
                    }else{ // i contains only x

                        Nx.push_back(i);
                    }

                }else if(std::find(setNeighbour.begin(), setNeighbour.end(), y) != setNeighbour.end()) // if y is in there (but not x)
                {
                    Ny.push_back(i);
                }else{ // nor x nor y are in there
                    N.push_back(i);
                }
            }

            // search Nx for x U S and select Cx if possible
            possibleComponents.clear(); 

            std::copy(setS.begin(), setS.end(),
                std::back_inserter(xUS));
            xUS.push_back(x);
            std::sort(xUS.begin(),xUS.end());

            for( auto i : Nx )
            {
                setNeighbour = i->getNodes();
                if( std::includes(setNeighbour.begin(), setNeighbour.end(), xUS.begin(), xUS.end()) ) // note that includes only works on sorted ranges
                    possibleComponents.push_back(i);
            }
            if( possibleComponents.size() > 0 )
            {
                Cx = possibleComponents[ randIntUniform(0,possibleComponents.size()-1) ];
                definedCx = true;
            }

            // search Ny for y U S and select Cy if possible
            possibleComponents.clear();

            std::copy(setS.begin(), setS.end(),
                std::back_inserter(yUS));
            yUS.push_back(y);
            std::sort(yUS.begin(),yUS.end());

            for( auto i : Ny )
            {
                setNeighbour = i->getNodes();
                if( std::includes(setNeighbour.begin(), setNeighbour.end(), yUS.begin(), yUS.end()) ) // note that includes only works on sorted ranges
                    possibleComponents.push_back(i);
            }
            
            if( possibleComponents.size() > 0 )
            {
                Cy = possibleComponents[ randIntUniform(0,possibleComponents.size()-1) ];
                definedCy = true;
            }

            // Now onto the 4 cases

            if( !definedCx && !definedCy ) // if both are undefined
            {
                // Type a)

                // Create two new JTComponents, separated by S and made up by xUS and yUS
                ClX = std::make_shared<JTComponent>( setS );
                ClX->add1Node(y);

                ClY = std::make_shared<JTComponent>( setS );
                ClY->add1Node(x);

                // decide which (Clx or Cly) is the left and which the right clique
                if( std::find(Nx.begin(), Nx.end(), C->getParent() ) != Nx.end() ) // the parent was in Nx, hence I need ClY on the left
                {
                    cLeft = ClY;  // here I'm copying the pointer, so I can use them interchangeably
                    cRight = ClX;

                }else if( std::find(Ny.begin(), Ny.end(), C->getParent() ) != Ny.end() ) // the parent was in Ny, hence I need ClX on the left
                {
                    cLeft = ClX;  
                    cRight = ClY;

                }else{ // either C was the root or the parent is in N, so I can choose randomly

                    if( randU01() < 0.5 )
                    {
                        cLeft = ClY;  
                        cRight = ClX;
                    }else{
                        cLeft = ClX;  
                        cRight = ClY;                        
                    }
                }

                // regardless, add now cLeft to its childrens and remove C
                if( C->getParent() )
                {
                    newChildrens = C->getParent()->getChildrens();
                    newChildrens.erase(std::remove(newChildrens.begin(), newChildrens.end(), C), newChildrens.end());
                    newChildrens.push_back( cLeft );
                    C->getParent()->setChildrens( newChildrens );
                } // else it was the root so no changes

                // Add the original parent to cLeft                
                cLeft->setParent( C->getParent() );
                // and its separator, which might be empty
                cLeft->setSeparator( C->getSeparator() );

                // the parent of the right clique is then the left clique 
                // (and the right is a child for the left)
                // their separator is S
                cLeft->add1Children( cRight ); // note that this is just ONE children and there were none before
                cRight->setParent( cLeft );
                cRight->setSeparator( setS );    

                // Now connect all the neighbours to either cLeft or cRight (except the original parent which is already connected)
                for( auto i : Nx )
                {
                    if( i != C->getParent() )
                    {
                        ClY->add1Children( i );
                        i->setParent( ClY );
                    }
                }

                for( auto i : Ny )
                {
                    if( i != C->getParent() )
                    {
                        ClX->add1Children( i );
                        i->setParent( ClX );
                    }
                }

                for( auto i : N )
                {
                    if( i != C->getParent() )
                    {
                        if( randU01() < 0.5 )
                        {
                            ClX->add1Children( i );
                            i->setParent( ClX );

                        }else{
                            ClY->add1Children( i );
                            i->setParent( ClY );
                        }
                    }
                }

                // Finally erase C from PCS .. (this actually gets overwritten so no need..)
                // perfectCliqueSequence.erase(perfectCliqueSequence.begin() + randomComp);

                // .. and Re-build the PCS after this
                newPCS.clear();
                if( !( C->getParent() ) )
                {
                    // cLeft is the new root
                    newPCS.insert( newPCS.end() , cLeft );
                }else{
                    // the old root is still the root, but I've modified its childrens
                    newPCS.insert( newPCS.end() , perfectCliqueSequence[0] );
                }
                pos = 0;
                buildNewPCS( newPCS, pos );
                perfectCliqueSequence = newPCS; // substitute to the current one
                // PEO updated below

                logP += N.size() * log(2.); // forward probability addition
                logP -= log((double)numComponents) + log((double)(ClX->getNodes().size() - setS.size())) + log((double)(ClY->getNodes().size() - setS.size())) ; //backward probability (we added one component here so numComponents rather than numComponents-1)

            }else if( definedCx && !definedCy ) // if only Cx is defined
            {
                // Type b)

                // Separation is possible only if Nx contains only Cx
                if( Nx.size() == 1 )
                {
                    if( Nx[0] == Cx )
                    {
                        // remove x from C
                        setC.erase(std::remove(setC.begin(), setC.end(), x), setC.end());
                        C->setNodes( setC );
                        // find the separator that connects C to Cx
                        if( Cx == C->getParent() )
                        {
                            setS = C->getSeparator();
                            setS.erase(std::remove(setS.begin(), setS.end(), x), setS.end());
                            C->setSeparator( setS );
                        }else{ // it's a child
                            setS = Cx->getSeparator();
                            setS.erase(std::remove(setS.begin(), setS.end(), x), setS.end());
                            Cx->setSeparator( setS );
                        }
                    }else{
                        return std::make_pair(false,0.0); // and the graph is unchanged
                    }
                }else{
                    return std::make_pair(false,0.0); // and the graph is unchanged
                }

                logP -= log((double)(numComponents-1.)) + log((double)(Cx->getNodes().size() - setS.size())) + log((double)(setC.size() - setS.size())) ; //backward probability

            }else if( !definedCx && definedCy ) // if only Cy is defined
            {
                // Type c)

                // Separation is possible only if Ny contains only Cy
                if( Ny.size() == 1 )
                {
                    if( Ny[0] == Cy )
                    {
                        // remove y from C
                        setC.erase(std::remove(setC.begin(), setC.end(), y), setC.end());
                        C->setNodes( setC );
                        // find the separator that connects C to Cy
                        if( Cy == C->getParent() )
                        {
                            setS = C->getSeparator();
                            setS.erase(std::remove(setS.begin(), setS.end(), y), setS.end());
                            C->setSeparator( setS );
                        }else{ // it's a child
                            setS = Cy->getSeparator();
                            setS.erase(std::remove(setS.begin(), setS.end(), y), setS.end());
                            Cy->setSeparator( setS );
                        }
                    }else{
                        return std::make_pair(false,0.0); // and the graph is unchanged
                    }
                }else{
                    return std::make_pair(false,0.0); // and the graph is unchanged
                }
                
                logP -= log((double)(numComponents-1.)) + log((double)(Cy->getNodes().size() - setS.size()) ) + log((double)(setC.size() - setS.size())) ; //backward probability
                
            }else{ //if both are defined
                // Type d)

                // Separation is possible only if N* contains only C* [*=x,y] and N is empty
                if( Nx.size() == 1 && Ny.size() == 1 && N.size() == 0 )
                {
                    if( Nx[0]==Cx && Ny[0]==Cy)
                    {
                        // Find who's first between Cx and Cy
                        // adjust their links
                        // and change their separator to be a new S
                        // the separator for cLeft now is either its old separator (if it was the parent)
                        // or it needs to be set to the separator of C if cLeft was its child
                        if( Cx == C->getParent() )
                        {
                            cLeft = Cx;
                            cRight = Cy;

                        }else if( Cy == C->getParent() ) // if Cy is the parent
                        {
                            cLeft = Cy;
                            cRight = Cx;

                        }else if( !( C->getParent()) ) // C was their parent, so choose randomly
                        {                               // note as well that in this case (becase N is empty), C was the root of JT
                            if( randU01() < 0.5 )
                            {
                                cLeft = Cx;
                                cRight = Cy;
                            }else{
                                cLeft = Cy;
                                cRight = Cx;
                            }

                            // because C was the root of the tree, the clique chosen to replace him at the top get its separator emptied
                            newSeparator.clear();
                            cLeft->setSeparator( newSeparator );

                            // and its parent's emptied as well as it becomes the new root
                            cLeft->setParent( nullptr );

                            // DO WE NEED TO CHANGE THE PCS in this case? (it's safer to do so, but... do we need?)
                            // I'll do it because our code relies on the fact that perfectCliqueSequence[0] is the root

                        }

                        // Now fix links
                        // find C in cLeft's childs and substitute (or add) cRight
                        newChildrens = cLeft->getChildrens();
                        newChildrens.erase(std::remove(newChildrens.begin(), newChildrens.end(), C), newChildrens.end());
                        newChildrens.push_back( cRight );
                        cLeft->setChildrens( newChildrens );

                        // set cLeft as cRight's parent (cRight HAS to come from C's childs)
                        cRight->setParent( cLeft );

                        // set the separator between them as the intersection between them
                        newSeparator.clear();
                        setCLeft = cLeft->getNodes();
                        setCRight = cRight->getNodes();

                        std::set_intersection(setCLeft.begin(), setCLeft.end(),
                                                setCRight.begin(), setCRight.end(),
                                            std::back_inserter(newSeparator));

                        cRight->setSeparator( newSeparator );
                        
                        // Finally erase C
                        perfectCliqueSequence.erase(perfectCliqueSequence.begin() + randomComp);

                        if( !( C->getParent()) ) // if we erased the root node, we created a new one and we need to create a newPCS
                        {
                            newPCS.clear();
                            // cLeft is the new root
                            newPCS.insert( newPCS.end() , cLeft );
                            pos = 0;
                            buildNewPCS( newPCS, pos );
                            perfectCliqueSequence = newPCS; // substitute to the current one
                            // PEO updated below
                        }

                    }else{
                        return std::make_pair(false,0.0); // and the graph is unchanged
                    }

                }else{
                    return std::make_pair(false,0.0); // and the graph is unchanged
                }
                
                logP -= log((double)(numComponents-2.)) + log((double)(setCLeft.size()-newSeparator.size())) + log((double)(setCRight.size()-newSeparator.size())) ; //backward probability (-2 here cause we further deleted one component)
        
            }

        }else //if only one there's no option for edge deletion
        {
            return std::make_pair(false,0.0); // and the graph is unchanged
        }

    }

    updatePEO(); // this->updatePEO();
    updateAdjMat();

    return std::make_pair(true,logP);

}

/*
JT Update proposal, from Green and Thomas 2013
MULTIPLE EDGES UPDATE
*/

std::pair<bool,double> JunctionTree::propose_multiple_edge_update( )
{
    // We need to select randomly one separator -- in our structure is any component beside the first
    unsigned int numComponents = perfectCliqueSequence.size();
    unsigned int randomSep, randomComp;

    arma::uvec randomIndexes;

    std::shared_ptr<JTComponent> Cx, Cy, C, ClX, ClY, cLeft, cRight;
    std::vector<unsigned int> setCx, setCy, setCxlS, setCylS, setS, 
        setC, setNeighbour, xUS, yUS, setCLeft, setCRight, tmpVec;

    std::vector<unsigned int> X,Y;
    unsigned int dimX, dimY; // dimension of the X and Y sets

    std::vector<std::shared_ptr<JTComponent>> newChildrens;
    std::vector<unsigned int> newNodes;
    std::vector<unsigned int> newSeparator;
    std::shared_ptr<JTComponent> newParent;

    std::vector<std::shared_ptr<JTComponent>> neighbours, N, Nx, Ny, possibleCxComponents, possibleCyComponents;
    bool definedCx=false, definedCy=false;
    std::deque<std::shared_ptr<JTComponent>> newPCS;
    bool XisRightYisLeft = false;
    unsigned int sizeXIntersect, sizeYIntersect;

    unsigned int pos = 0;
    unsigned int countN = 0;

    double logP = 0;

    if( randU01() < 0.5 )
    {
        // propose addition
        // if there's one component only, return false and no update happened
        if( numComponents > 1 )
        {
            // get a Random separator
            randomSep = randIntUniform(1,numComponents-1);

            // populate Cx and Cy
            Cx = perfectCliqueSequence[randomSep]->getParent();
            Cy = perfectCliqueSequence[randomSep];

            setCx = Cx->getNodes();
            setCy = Cy->getNodes();
            setS = perfectCliqueSequence[randomSep]->getSeparator();
            
            // Populate Cx\S and Cy\S
            std::set_difference(setCx.begin(), setCx.end(), setS.begin(), setS.end(), 
                            std::inserter(setCxlS, setCxlS.begin()));

            std::set_difference(setCy.begin(), setCy.end(), setS.begin(), setS.end(), 
                            std::inserter(setCylS, setCylS.begin()));

            // now choose x and y from Cx\S and Cy\S
            dimX = randIntUniform(1,setCxlS.size());
            dimY = randIntUniform(1,setCylS.size());

            X = Distributions::randSampleWithoutReplacement( setCxlS.size() , setCxlS , dimX ); 
            Y = Distributions::randSampleWithoutReplacement( setCylS.size() , setCylS , dimY ); 
                // this is uncorrect, we should sample
                // randomly from all possible subsets of CxlS

            // check whether or not Cx and Cy are superset of x U S and y U S and act accordingly
            // note there's no way for it to be the other way around by construction


            // forward probability (common for all merge moves)
            logP += log((double)(numComponents-1.)) +
                log((double)(setCxlS.size())) - std::lgamma((double)(dimX+1.)) - std::lgamma((double)(setCxlS.size()-dimX+1.)) + std::lgamma((double)(setCxlS.size()+1.)) +
                log((double)(setCylS.size())) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(setCylS.size()-dimY+1.)) + std::lgamma((double)(setCylS.size()+1.)) ;


            if( setCxlS.size() == X.size() && setCylS.size() == Y.size() ) // a) this means that both Cx and Cy are exactly just X and S (and Y and S) 
            {
                
                // Type a)

                // compute new children set
                newChildrens = Cy->getChildrens();
                for( auto c : Cx->getChildrens() )
                {
                    if( c != Cy )
                        newChildrens.push_back(c);
                }

                newNodes = setS; 
                for( auto i : X )
                    newNodes.push_back(i);
                for( auto i : Y )
                    newNodes.push_back(i);

                // newSeparator = Cx->getSeparator();
                // newParent = Cx->getParent();

                // Modify Cx to C* (the new component we need)
                Cx->setChildrens(newChildrens);
                Cx->setNodes(newNodes);
                // Cx->setSeparator(newSeparator);
                // Cx->setParent(newParent);

                // make the new childrens point to C* (rather than Cy)
                for( auto c : newChildrens )
                {
                    if( c->getParent() != Cx )
                        c->setParent(Cx);
                }

                // Erase Cy
                perfectCliqueSequence.erase(perfectCliqueSequence.begin() + randomSep);

                // **** logP addition
                neighbours.clear();
                if( Cx->getParent() )
                    neighbours.push_back( Cx->getParent() );  // push parent only if it exists

                for( auto i : Cx->getChildrens() )
                    neighbours.push_back(i);

                for( auto i : neighbours )
                {
                    setNeighbour = i->getNodes();

                    // // compute intersection between this neighbour and X
                    // tmpVec.clear();
                    // std::set_intersection(setNeighbour.begin(), setNeighbour.end(),
                    //         X.begin(), X.end(),
                    //         std::back_inserter(tmpVec));

                    // // if size is more than 0, then it intersects
                    // // if size is equal to dimX then the complete X is in "i" (neighbour)
                    // sizeXIntersect = tmpVec.size();

                    // // same for Y
                    // tmpVec.clear();
                    // std::set_intersection(setNeighbour.begin(), setNeighbour.end(),
                    //         Y.begin(), Y.end(),
                    //         std::back_inserter(tmpVec));

                    // // if size is more than 0, then it intersects
                    // // if size is equal to dimY then the complete Y is in "i" (neighbour)
                    // sizeYIntersect = tmpVec.size();
                    //
                    // if( sizeXIntersect == 0 && sizeYIntersect == 0 ) // if neither intersect i
                    //     ++countN;

                    // go to the next "i" as soon as you find one match in X
                    for( auto j : X )
                        if( std::find(setNeighbour.begin(), setNeighbour.end(), j ) != setNeighbour.end() ) // the element j in X is in neighbour i
                            continue;

                    // go to the next "i" as soon as you find one match in Y
                    for( auto j : Y )
                        if( std::find(setNeighbour.begin(), setNeighbour.end(), j ) != setNeighbour.end() ) // the element j in Y is in neighbour i
                            continue;

                    //  if not out of the loop yet
                    ++countN;

                }

                logP -= countN * log(2.); // backward probability addition
                logP -= // backward probability (-1 because we reduce the # component by 1)
                        log((double)(numComponents-1.)) - log(2.) + log( (double)(Cx->getNodes().size()-1.) ) + log( (double)(dimX+dimY-1.) ) -
                            std::lgamma((double)(dimX+1.)) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(Cx->getSeparator().size()+1.)) + std::lgamma((double)(Cx->getNodes().size()+1.));

            }else if( setCxlS.size() > X.size() && setCylS.size() == Y.size() ) // b)
            {

                // Type b)
                Cy->addNodes(X);
                Cy->addSeparators(X);

                logP -= // backward probability
                        log((double)numComponents) - log(2.) + ( log( (double)(Cy->getNodes().size()-1.) ) + log( (double)(dimX+dimY-1.) ) ) -
                            std::lgamma((double)(dimX+1.)) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(Cy->getSeparator().size()+1.)) + std::lgamma((double)(Cy->getNodes().size()+1.));

            }else if( setCxlS.size() == X.size() && setCylS.size() > Y.size() ) //c)
            {
                // Type c)


                Cx->addNodes(Y);
                Cy->addSeparators(Y); // remember the separator is always the one from the Cy obj

                logP -= // backward probability
                        log((double)numComponents) - log(2.) + ( log( (double)(Cx->getNodes().size()-1.) ) + log( (double)(dimX+dimY-1.) ) ) -
                            std::lgamma((double)(dimX+1.)) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(Cy->getSeparator().size()+1.)) + std::lgamma((double)(Cx->getNodes().size()+1.));

            }else if( setCxlS.size() > X.size() && setCylS.size() > Y.size() ) // d) Cx and Cy contain more than just x and S (and y and S) 
            {
                // Type d)

                // Create C* (the new Component)
                newChildrens = std::vector<std::shared_ptr<JTComponent>>(1);
                newChildrens[0] = Cy;

                newNodes = setS; 
                for( auto i : X )
                    newNodes.push_back(i);
                for( auto i : Y )
                    newNodes.push_back(i);

                newSeparator = setS;
                for( auto i : X ) 
                    newSeparator.push_back(i);
                
                newParent = Cx;

                // insert it in the PCS
                perfectCliqueSequence.insert( perfectCliqueSequence.begin() + randomSep ,
                                std::make_shared<JTComponent>(newNodes,newSeparator,newChildrens,newParent) );

                // Modify Cx to point to C*
                newChildrens = Cx->getChildrens();
                for (auto itC = newChildrens.begin(); itC != newChildrens.end();  ++itC )
                {
                    if( *(itC) == Cy )
                        *(itC) = perfectCliqueSequence[randomSep];
                }
                Cx->setChildrens(newChildrens);

                // Modify Cy to have C* as parent and modify its separator
                Cy->setParent(perfectCliqueSequence[randomSep]);

                Cy->addSeparators(Y);         


                logP -= // backward probability (+1 because we insert a new component here)
                    log((double)(numComponents+1.)) - log(2.) + log( (double)(perfectCliqueSequence[randomSep]->getNodes().size()-1.) ) + log( (double)(dimX+dimY-1.) ) -
                        std::lgamma((double)(dimX+1.)) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(perfectCliqueSequence[randomSep]->getSeparator().size()+1.)) + std::lgamma((double)(perfectCliqueSequence[randomSep]->getNodes().size()+1.));
                
            }

        }else //if only one there's no option for edge addition
        {
            return std::make_pair(false,0.0); // and the graph is unchanged
        }
    
    }else         // propose deletion
    {
        // get a random component
        randomComp = randIntUniform(0,numComponents-1);
        C = perfectCliqueSequence[randomComp];

        setC = C->getNodes();

        if( setC.size() > 1 )
        {
            // partition C into 3 sets, x y and S (S might be empty)"
            dimY = randIntUniform(2,setC.size());
            dimX = randIntUniform(1,dimY-1);
            dimY -= dimX;

            // Set S = Set C
            std::copy(setC.begin(), setC.end(),
                std::back_inserter(setS));
            // Select X
            X = Distributions::randSampleWithoutReplacement( setS.size() , setS , dimX ); 

            //remove X from setS
            tmpVec.clear();
            std::set_difference(std::make_move_iterator(setS.begin()), 
                                std::make_move_iterator(setS.end()), 
                                X.begin(), X.end(), 
                        std::inserter(tmpVec, tmpVec.begin()));
            setS.swap(tmpVec);

            // Select Y
            Y = Distributions::randSampleWithoutReplacement( setS.size() , setS , dimY ); 

            //remove Y from setS
            tmpVec.clear();
            std::set_difference(std::make_move_iterator(setS.begin()), 
                                std::make_move_iterator(setS.end()), 
                                Y.begin(), Y.end(), 
                        std::inserter(tmpVec, tmpVec.begin()));
            setS.swap(tmpVec);

            // what's left in setS is just S

            // forward logP computed immediately as setC is likely to be modified after
            logP += log((double)numComponents) - log(2.) + ( log( (double)(setC.size()-1.) ) + log( (double)(dimX+dimY-1.) ) ) -
                std::lgamma((double)(dimX+1.)) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(setS.size()+1.)) + std::lgamma((double)(setC.size()+1.)); // forward probability

            // now scan the neighbours of C and construct Nx, Ny and N"
            neighbours.clear();
            if( C->getParent() )
                neighbours.push_back( C->getParent() );  // push parent only if it exists

            for( auto i : C->getChildrens() )
                neighbours.push_back(i);


            possibleCxComponents.clear();
            possibleCyComponents.clear();

            for( auto i : neighbours )
            {
                setNeighbour = i->getNodes();

                // compute intersection between this neighbour and X
                tmpVec.clear();
                std::set_intersection(setNeighbour.begin(), setNeighbour.end(),
                          X.begin(), X.end(),
                          std::back_inserter(tmpVec));

                // if size is more than 0, then it intersects
                // if size is equal to dimX then the complete X is in "i" (neighbour)
                sizeXIntersect = tmpVec.size();

                // same for Y
                tmpVec.clear();
                std::set_intersection(setNeighbour.begin(), setNeighbour.end(),
                          Y.begin(), Y.end(),
                          std::back_inserter(tmpVec));

                // if size is more than 0, then it intersects
                // if size is equal to dimY then the complete Y is in "i" (neighbour)
                sizeYIntersect = tmpVec.size();

                if( sizeXIntersect > 0 && sizeYIntersect > 0 ) // if both intersect i
                {
                    return std::make_pair(false,0.0); // exit, deletion is not possible
                    
                }else if( sizeXIntersect > 0 && sizeYIntersect == 0 ) // i contains only x
                {
                    Nx.push_back(i);
                    if ( sizeXIntersect == dimX )
                    {
                        // possible components already have all X, so we need to be sure they contain all of S as well
                        if( std::includes(setNeighbour.begin(), setNeighbour.end(), setS.begin(), setS.end()) ) // note that includes only works on sorted ranges
                            possibleCxComponents.push_back(i);
                    }

                }else if( sizeXIntersect == 0 && sizeYIntersect > 0 ) // if y is in there (but not x)
                {
                    Ny.push_back(i);
                    if ( sizeYIntersect == dimY )
                    {
                        // possible components already have all X, so we need to be sure they contain all of S as well
                        if( std::includes(setNeighbour.begin(), setNeighbour.end(), setS.begin(), setS.end()) ) // note that includes only works on sorted ranges
                            possibleCyComponents.push_back(i);
                    }

                }else{ // nor x nor y are in there
                    N.push_back(i);
                }
            }

            // If we have candidates, select one at random as Cx and/or Cy
            if( possibleCxComponents.size() > 0 )
            {
                Cx = possibleCxComponents[ randIntUniform(0,possibleCxComponents.size()-1) ];
                definedCx = true;
            }

            if( possibleCyComponents.size() > 0 )
            {
                Cy = possibleCyComponents[ randIntUniform(0,possibleCyComponents.size()-1) ];
                definedCy = true;
            }

            // Now onto the 4 cases

            if( !definedCx && !definedCy ) // if both are undefined
            {

                // Type a)

                // Create two new JTComponents, separated by S and made up by xUS and yUS
                ClX = std::make_shared<JTComponent>( setS );
                ClX->addNodes(Y);

                ClY = std::make_shared<JTComponent>( setS );
                ClY->addNodes(X);

                // decide which (Clx or Cly) is the left and which the right clique
                if( std::find(Nx.begin(), Nx.end(), C->getParent() ) != Nx.end() ) // the parent was in Nx, hence I need ClY on the left
                {
                    cLeft = ClY;  // here I'm copying the pointer, so I can use them interchangeably
                    cRight = ClX;

                }else if( std::find(Ny.begin(), Ny.end(), C->getParent() ) != Ny.end() ) // the parent was in Ny, hence I need ClX on the left
                {
                    cLeft = ClX;  
                    cRight = ClY;

                }else{ // either C was the root or the parent is in N, so I can choose randomly

                    if( randU01() < 0.5 )
                    {
                        cLeft = ClY;  
                        cRight = ClX;
                    }else{
                        cLeft = ClX;  
                        cRight = ClY;                        
                    }
                }

                // regardless, add now cLeft to its childrens and remove C
                if( C->getParent() )
                {
                    newChildrens = C->getParent()->getChildrens();
                    newChildrens.erase(std::remove(newChildrens.begin(), newChildrens.end(), C), newChildrens.end());
                    newChildrens.push_back( cLeft );
                    C->getParent()->setChildrens( newChildrens );
                } // else it was the root so no changes

                // Add the original parent to cLeft                
                cLeft->setParent( C->getParent() );
                // and its separator, which might be empty
                cLeft->setSeparator( C->getSeparator() );

                // the parent of the right clique is then the left clique 
                // (and the right is a child for the left)
                // their separator is S
                cLeft->add1Children( cRight ); // note that this is just ONE children and there were none before
                cRight->setParent( cLeft );
                cRight->setSeparator( setS );    

                // Now connect all the neighbours to either cLeft or cRight (except the original parent which is already connected)
                for( auto i : Nx )
                {
                    if( i != C->getParent() )
                    {
                        ClY->add1Children( i );
                        i->setParent( ClY );
                    }
                }

                for( auto i : Ny )
                {
                    if( i != C->getParent() )
                    {
                        ClX->add1Children( i );
                        i->setParent( ClX );
                    }
                }

                for( auto i : N )
                {
                    if( i != C->getParent() )
                    {
                        if( randU01() < 0.5 )
                        {
                            ClX->add1Children( i );
                            i->setParent( ClX );

                        }else{
                            ClY->add1Children( i );
                            i->setParent( ClY );
                        }
                    }
                }

                // Finally erase C from PCS .. (this actually gets overwritten so no need..)
                // perfectCliqueSequence.erase(perfectCliqueSequence.begin() + randomComp);

                // .. and Re-build the PCS after this
                newPCS.clear();
                if( !( C->getParent() ) )
                {
                    // cLeft is the new root
                    newPCS.insert( newPCS.end() , cLeft );
                }else{
                    // the old root is still the root, but I've modified its childrens
                    newPCS.insert( newPCS.end() , perfectCliqueSequence[0] );
                }
                pos = 0;
                buildNewPCS( newPCS, pos );
                perfectCliqueSequence = newPCS; // substitute to the current one
                // PEO updated below

                logP += N.size() * log(2.); // forward probability addition
                logP -= //backward probability (we added one component here so numComponents rather than numComponents-1)
                        log((double)numComponents) +
                        log((double)(ClY->getNodes().size() - setS.size())) - std::lgamma((double)(dimX+1.)) - std::lgamma((double)(ClY->getNodes().size()-setS.size()-dimX+1.)) + std::lgamma((double)(ClY->getNodes().size()-setS.size()+1.)) +
                        log((double)(ClX->getNodes().size() - setS.size())) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(ClX->getNodes().size()-setS.size()-dimY+1.)) + std::lgamma((double)(ClX->getNodes().size()-setS.size()+1.)) ;

            }else if( definedCx && !definedCy ) // if only Cx is defined
            {
                // Type b)

                // Separation is possible only if Nx contains only Cx
                if( Nx.size() == 1 )
                {
                    if( Nx[0] == Cx )
                    {
                        // remove x from C
                        tmpVec.clear();
                        std::set_difference(std::make_move_iterator(setC.begin()), 
                                            std::make_move_iterator(setC.end()), 
                                            X.begin(), X.end(), 
                                    std::inserter(tmpVec, tmpVec.begin()));
                        setC.swap(tmpVec); // move operator, move semantic and swap...ah
                    
                        C->setNodes( setC );

                        // find the separator that connects C to Cx
                        if( Cx == C->getParent() )
                        {
                            setS = C->getSeparator();

                            tmpVec.clear();
                            std::set_difference(std::make_move_iterator(setS.begin()), 
                                            std::make_move_iterator(setS.end()), 
                                            X.begin(), X.end(), 
                                    std::inserter(tmpVec, tmpVec.begin()));
                            setS.swap(tmpVec);

                            C->setSeparator( setS );

                        }else{ // it's a child

                            setS = Cx->getSeparator();
                            
                            tmpVec.clear();
                            std::set_difference(std::make_move_iterator(setS.begin()), 
                                            std::make_move_iterator(setS.end()), 
                                            X.begin(), X.end(), 
                                    std::inserter(tmpVec, tmpVec.begin()));
                            setS.swap(tmpVec);
                            
                            Cx->setSeparator( setS );
                        }
                    }else{
                        return std::make_pair(false,0.0); // and the graph is unchanged
                    }
                }else{
                    return std::make_pair(false,0.0); // and the graph is unchanged
                }

                logP -= //backward probability
                        log((double)(numComponents-1.)) +
                        log((double)(Cx->getNodes().size() - setS.size())) - std::lgamma((double)(dimX+1.)) - std::lgamma((double)(Cx->getNodes().size()-setS.size()-dimX+1.)) + std::lgamma((double)(Cx->getNodes().size()-setS.size()+1.)) +
                        log((double)(setC.size() - setS.size())) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(setC.size()-setS.size()-dimY+1.)) + std::lgamma((double)(setC.size()-setS.size()+1.)) ;

            }else if( !definedCx && definedCy ) // if only Cy is defined
            {
                // Type c)

                // Separation is possible only if Ny contains only Cy
                if( Ny.size() == 1 )
                {
                    if( Ny[0] == Cy )
                    {
                        // remove y from C
                        tmpVec.clear();
                        std::set_difference(std::make_move_iterator(setC.begin()), 
                                            std::make_move_iterator(setC.end()), 
                                            Y.begin(), Y.end(), 
                                    std::inserter(tmpVec, tmpVec.begin()));
                        setC.swap(tmpVec); // move operator, move semantic and swap...ah

                        C->setNodes( setC );
                        // find the separator that connects C to Cy
                        if( Cy == C->getParent() )
                        {
                            setS = C->getSeparator();
                            
                            tmpVec.clear();
                            std::set_difference(std::make_move_iterator(setS.begin()), 
                                            std::make_move_iterator(setS.end()), 
                                            Y.begin(), Y.end(), 
                                    std::inserter(tmpVec, tmpVec.begin()));
                            setS.swap(tmpVec);

                            C->setSeparator( setS );
                        }else{ // it's a child
                            setS = Cy->getSeparator();

                            tmpVec.clear();
                            std::set_difference(std::make_move_iterator(setS.begin()), 
                                            std::make_move_iterator(setS.end()), 
                                            Y.begin(), Y.end(), 
                                    std::inserter(tmpVec, tmpVec.begin()));
                            setS.swap(tmpVec);

                            Cy->setSeparator( setS );
                        }
                    }else{
                        return std::make_pair(false,0.0); // and the graph is unchanged
                    }
                }else{
                    return std::make_pair(false,0.0); // and the graph is unchanged
                }
                
                logP -= //backward probability
                    log((double)(numComponents-1.)) +
                    log((double)(setC.size() - setS.size())) - std::lgamma((double)(dimX+1.)) - std::lgamma((double)(setC.size()-setS.size()-dimX+1.)) + std::lgamma((double)(setC.size()-setS.size()+1.)) +
                    log((double)(Cy->getNodes().size() - setS.size())) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(Cy->getNodes().size()-setS.size()-dimY+1.)) + std::lgamma((double)(Cy->getNodes().size()-setS.size()+1.));


            }else{ //if both are defined
                // Type d)

                // Separation is possible only if N* contains only C* [*=x,y] and N is empty
                if( Nx.size() == 1 && Ny.size() == 1 && N.size() == 0 )
                {
                    if( Nx[0]==Cx && Ny[0]==Cy)
                    {
                        // Find who's first between Cx and Cy
                        // adjust their links
                        // and change their separator to be a new S
                        // the separator for cLeft now is either its old separator (if it was the parent)
                        // or it needs to be set to the separator of C if cLeft was its child
                        if( Cx == C->getParent() )
                        {
                            cLeft = Cx;
                            cRight = Cy;

                            XisRightYisLeft = false;

                        }else if( Cy == C->getParent() ) // if Cy is the parent
                        {
                            cLeft = Cy;
                            cRight = Cx;

                            XisRightYisLeft = true;

                        }else if( !( C->getParent()) ) // C was their parent, so choose randomly
                        {                               // note as well that in this case (becase N is empty), C was the root of JT
                            if( randU01() < 0.5 )
                            {
                                cLeft = Cx;
                                cRight = Cy;
                                XisRightYisLeft = false;

                            }else{
                                cLeft = Cy;
                                cRight = Cx;
                                XisRightYisLeft = true;
                            }

                            // because C was the root of the tree, the clique chosen to replace him at the top get its separator emptied
                            newSeparator.clear();
                            cLeft->setSeparator( newSeparator );

                            // and its parent's emptied as well as it becomes the new root
                            cLeft->setParent( nullptr );

                            // DO WE NEED TO CHANGE THE PCS in this case? (it's safer to do so, but... do we need?)
                            // I'll do it because our code relies on the fact that perfectCliqueSequence[0] is the root

                        }

                        // Now fix links
                        // find C in cLeft's childs and substitute (or add) cRight
                        newChildrens = cLeft->getChildrens();
                        newChildrens.erase(std::remove(newChildrens.begin(), newChildrens.end(), C), newChildrens.end());
                        newChildrens.push_back( cRight );
                        cLeft->setChildrens( newChildrens );

                        // set cLeft as cRight's parent (cRight HAS to come from C's childs)
                        cRight->setParent( cLeft );

                        // set the separator between them as the intersection between them
                        newSeparator.clear();
                        setCLeft = cLeft->getNodes();
                        setCRight = cRight->getNodes();

                        std::set_intersection(setCLeft.begin(), setCLeft.end(),
                                                setCRight.begin(), setCRight.end(),
                                            std::back_inserter(newSeparator));

                        cRight->setSeparator( newSeparator );
                        
                        // Finally erase C
                        perfectCliqueSequence.erase(perfectCliqueSequence.begin() + randomComp);

                        if( !( C->getParent()) ) // if we erased the root node, we created a new one and we need to create a newPCS
                        {
                            newPCS.clear();
                            // cLeft is the new root
                            newPCS.insert( newPCS.end() , cLeft );
                            pos = 0;
                            buildNewPCS( newPCS, pos );
                            perfectCliqueSequence = newPCS; // substitute to the current one
                            // PEO updated below
                        }

                        if( XisRightYisLeft )
                        {
                            logP -= //backward probability (-2 here cause we further deleted one component)
                                log((double)(numComponents-2.)) +
                                log((double)(setCLeft.size() - newSeparator.size())) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(setCLeft.size()-newSeparator.size()-dimY+1.)) + std::lgamma((double)(setCLeft.size()-newSeparator.size()+1.)) +
                                log((double)(setCRight.size() - newSeparator.size())) - std::lgamma((double)(dimX+1.)) - std::lgamma((double)(setCRight.size()-newSeparator.size()-dimX+1.)) + std::lgamma((double)(setCRight.size()-newSeparator.size()+1.));
                        }else{
                            logP -= //backward probability (-2 here cause we further deleted one component)
                                log((double)(numComponents-2.)) +
                                log((double)(setCLeft.size() - newSeparator.size())) - std::lgamma((double)(dimX+1.)) - std::lgamma((double)(setCLeft.size()-newSeparator.size()-dimX+1.)) + std::lgamma((double)(setCLeft.size()-newSeparator.size()+1.)) +
                                log((double)(setCRight.size() - newSeparator.size())) - std::lgamma((double)(dimY+1.)) - std::lgamma((double)(setCRight.size()-newSeparator.size()-dimY+1.)) + std::lgamma((double)(setCRight.size()-newSeparator.size()+1.));
                        }
                        

                    }else{
                        return std::make_pair(false,0.0); // and the graph is unchanged
                    }

                }else{
                    return std::make_pair(false,0.0); // and the graph is unchanged
                }
                

            }

        }else //if only one there's no option for edge deletion
        {
            return std::make_pair(false,0.0); // and the graph is unchanged
        }

    }

    updatePEO(); // this->updatePEO();
    updateAdjMat();

    return std::make_pair(true,logP);

}



void JunctionTree::swapParentChild( std::shared_ptr<JTComponent>& parent , std::shared_ptr<JTComponent>& child )
{

    std::vector<std::shared_ptr<JTComponent>> childrenVec;
    std::shared_ptr<JTComponent> hyperParent;

    if( parent->getParent() )
    {
        hyperParent = parent->getParent();
        swapParentChild( hyperParent , parent );

    }
    
    // parent is now the root of the tree, so now swap
    
    child -> add1Children(parent);
    child -> setParent(nullptr);

    parent->setParent(child);
    
    childrenVec = parent->getChildrens();
    childrenVec.erase(std::remove(childrenVec.begin(), childrenVec.end(), child), childrenVec.end());
    parent -> setChildrens( childrenVec );

    parent -> setSeparator( child->getSeparator() );
    child -> clearSeparator();

}

void JunctionTree::reRoot()
{
    unsigned int pos = 0;
    std::deque<std::shared_ptr<JTComponent>> newPCS;
    std::shared_ptr<JTComponent> parent;

    // select at random one component as the NEW root ( sampling from 1 excludes current root )
    newPCS.insert( newPCS.end() , perfectCliqueSequence[ randIntUniform( 1 , perfectCliqueSequence.size() -1 ) ] );

    // get its parent
    parent = newPCS[pos]->getParent();

    // make the selected component recursively the root by swapping parents and childs
    swapParentChild( parent , newPCS[pos] );

    // now build the new PCS, starting from the root and its (new and old) childrens ..
    buildNewPCS( newPCS, pos );
    // ... and substitute to the current one
    perfectCliqueSequence = newPCS; 

    // updateAdjMat(); // this shouldn't be necessary !
    updatePEO();
    
}

bool JunctionTree::isChild( std::shared_ptr<JTComponent>& parent , std::shared_ptr<JTComponent>& node )
{
    std::vector<std::shared_ptr<JTComponent>> parentsChildrens = parent -> getChildrens();
    const unsigned int nChildrens = parentsChildrens.size();

    if( nChildrens == 0 )
        return false;

    std::vector<bool> res;
    res.resize(nChildrens,false); // init and set all to false

    for( unsigned int i=0; i<nChildrens; ++i )
    {
        if( parentsChildrens[i] == node )
            return true;
        else{
            res[i] = isChild( parentsChildrens[i] , node );
        }
    }

    for( auto const& b : res)
        if(b) return true;
    return false;

}

void JunctionTree::randomJTPermutation()
{
    unsigned int numComponents = perfectCliqueSequence.size();

    // this->print();

    // reroot the tree, this happens anyway
    if( numComponents > 1 ) // otherwise is at best the reverse of the above reRoot move
        this->reRoot();

    // then see if we can shuffle something
    if( numComponents > 2 ) // otherwise is at best the reverse of thea bove reRoot move
    {
        // Select at random one component (not the root)
        std::shared_ptr<JTComponent> thisComponent = perfectCliqueSequence[ randIntUniform( 1 , numComponents -1 ) ];
        std::shared_ptr<JTComponent> itsParent = thisComponent -> getParent();
        std::vector<std::shared_ptr<JTComponent>> itsChildrens = thisComponent -> getChildrens();
        std::vector<unsigned int> setI;
        std::vector<unsigned int> itsSeparator = thisComponent->getSeparator();
        std::vector<std::shared_ptr<JTComponent>> parentsChildrens;


        std::vector<std::shared_ptr<JTComponent>> possibleNewParent;
        unsigned int idx;

        // search in the PCS if there are other possible fathers (not in its childrens)
        bool isSelf, isOGParent, isItsChild, hasRightSep;
        for( unsigned int i=0; i<numComponents; ++i )
        {
            isSelf = ( perfectCliqueSequence[i] == thisComponent );

            if( !isSelf )
            {
                isOGParent = ( perfectCliqueSequence[i] == itsParent );

                if( !isOGParent )
                {
                    isItsChild = this->isChild( thisComponent , perfectCliqueSequence[i] );

                    if( !isItsChild )
                    {
                        setI = perfectCliqueSequence[i]->getNodes();
                        hasRightSep = std::includes(setI.begin(), setI.end(), itsSeparator.begin(), itsSeparator.end()); 
                        // should I check if the intersection between setI and setThis is just setS? no right?
                        // because if that was the case, the RIProperty would be violated from the start

                        if( hasRightSep )
                        {
                            // select this as a candidate
                            possibleNewParent.push_back( perfectCliqueSequence[i] );
                        }
                    }
                }
            }
        }

        // now select a new parent
        if( possibleNewParent.size() > 0 )
        {
            idx = randIntUniform(0,possibleNewParent.size()-1);

            //swap its parent with this new one, break and rebuild PCS
            thisComponent->setParent( possibleNewParent[idx] );
            possibleNewParent[idx] -> add1Children( thisComponent );
            
            parentsChildrens = itsParent -> getChildrens();
            parentsChildrens.erase(std::remove(parentsChildrens.begin(), parentsChildrens.end(), thisComponent), parentsChildrens.end());
            itsParent -> setChildrens( parentsChildrens );

            // Rebuild the PCS
            unsigned int pos = 0;
            std::deque<std::shared_ptr<JTComponent>> newPCS;

            newPCS.insert( newPCS.end() , perfectCliqueSequence[ 0 ] );

            buildNewPCS( newPCS , pos );

            // updateAdjMat(); // this shouldn't be necessary !
            updatePEO();

        } //else nothing to do  
    }
}
