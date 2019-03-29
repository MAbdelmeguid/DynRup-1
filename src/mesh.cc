#include "mesh.h"
#include "EigLab.h"

extern "C"
{
  #include <parmetis.h>
}

using namespace std;

/*****************************************************************/
/**                     Default constructor                     **/
/*****************************************************************/
Mesh::Mesh()
{
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &nprocs);
}

/*****************************************************************/
/**                     Read in Abaqus mesh                     **/
/*****************************************************************/
PetscErrorCode Mesh::Abaqus_IO(std::string &fname)
{
  PetscErrorCode err = 0;
  if (myid == 0)
  {
    std::string s;
    std::ifstream _in(fname);
    while (true)
    {
      std::getline(_in,s);
      // Convert s to uppercase
      std::string upper(s);
      std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
      // 0.) Look for the "*Part" Section
//      if (upper.find("*PART")== static_cast<std::string::size_type>(0))
//      {
//        std::cout<<"Find *PART"<<std::endl;
//      }
      // 1.) Loop for the "*Nodes" section
      if (upper.find("*NODE")== static_cast<std::string::size_type>(0))
      {
        //std::string nset_name = s ;
        //std::cout<<s<<std::endl;
        // Temperatry variables for parsing lines of text
        char c ;
        std::string line;
        int i=0;
        while (_in.peek()!='*'&&_in.peek()!=EOF)
        {
          // Read an entire line which corresponds to a single points's id and (x,y) value
          std::getline(_in,line);
          // Revomie all whitesspaces characters from the line
          line.erase(std::remove_if(line.begin(),line.end(),::isspace),line.end());
          // Make a stream out of the modified line so we can stream values from it in the usaly way
          std::stringstream ss(line);
          int abaqus_Node_id = 0;
          double x = 0 , y = 0;
          ss >> abaqus_Node_id >> c >> x >>c >> y ;
          Node.conservativeResize(i+1, 2);
          Node.row(i)<< x , y ;
          i++;
        }
        //  std::cout<< Nodes<< "\n"<< Nodes.rows()<<std::endl;
      }
      else if (upper.find("*ELEMENT,")==static_cast<std::string::size_type>(0))
      {
        //std::string elset_name = s;
        //std::cout<<s <<std::endl;
        char c ;
        std::string line;
        int i = 0;
        while (_in.peek()!='*'&&_in.peek()!=EOF)
        {
          // Read an entire line which corresponds to a single Element's id and connectivity value (Q4_only)
          std::getline(_in,line);
          // Revomie all whitesspaces characters from the line
          line.erase(std::remove_if(line.begin(),line.end(),::isspace),line.end());
          // Make a stream out of the modified line so we can stream values from it in the usaly way
          std::stringstream ss(line);
          int abaqus_el_id = 0;
          int Node1 =0 , Node2 =0 , Node3 = 0 , Node4 =0;
          ss >> abaqus_el_id >> c >> Node1 >> c >>Node2 >> c >> Node3 >> c >> Node4 ;
          // add -1 here becasue we are using a zero Node numbering zero is the first Node
          Element.conservativeResize(i+1, 4);
          Element.row(i)<<Node1-1,Node2-1,Node3-1,Node4-1;
          i++;
        }
      }
      if (_in.eof())
        break;
    }

    // Set distribution variables
    nLocElem = Element.rows();
    nLocNode = Node.rows();
    elDist.setConstant(nprocs+1, nLocElem);
    ndDist.setConstant(nprocs+1, nLocNode);
    elDist(0) = 0;
    ndDist(0) = 0;
  }
  else
  {
    // Set distribution variables
    nLocElem = 0;
    nLocNode = 0;
    elDist.resize(nprocs+1);
    ndDist.resize(nprocs+1);
  }

  // Distribute Element and Node distribution arrays
  err = MPI_Bcast(elDist.data(), nprocs+1, MPI_INT, 0, comm);
  err = MPI_Bcast(ndDist.data(), nprocs+1, MPI_INT, 0, comm);
  nElem = elDist(nprocs);
  nNode = ndDist(nprocs);
  nDims = 2; // TODO generalize this to ND meshes

  return err;
}

/*****************************************************************/
/**                    Redistribute Elements                    **/
/*****************************************************************/
PetscErrorCode Mesh::Redistribute()
{
  if (nprocs == 1)
    return 0;

  PetscErrorCode err = 0;
  if (elDist(1) == elDist(nprocs))
  {
    // Everything on processor 1, use METIS
    err = ReorderMETIS(); CHKERRQ(err);
  }
  else if ((elDist.segment(1, nprocs) - elDist.segment(0, nprocs)).minCoeff() == 0)
  {
    // Some processors are empty, ParMETIS will fail
    err = MPI_Abort(comm, 404);
  }
  else
  {
    // We can use ParMETIS
    err = ReorderParMETIS(); CHKERRQ(err);
  }

  err = NodeDist(); CHKERRQ(err);
  gElem = ArrayXI::LinSpaced(this->nLocElem, elDist(myid), elDist(myid+1)-1);
  gNode = ArrayXI::LinSpaced(this->nLocNode, ndDist(myid), ndDist(myid+1)-1);
  err = Expand_Node(); CHKERRQ(err);
  err = Initialize_Vectors(); CHKERRQ(err);

 /// Local Element Numbering
  err = Localize(); CHKERRQ(err);

  return err;
}

/*****************************************************************/
/**            Get element partitioning in serial               **/
/*****************************************************************/
PetscErrorCode Mesh::ReorderMETIS(PetscInt nparts, PetscInt nCommonNodes,
		                  PetscScalar *tpwgts, PetscInt *elmwgt,
				  PetscInt *opts)
{
  PetscErrorCode err = 0;
  ArrayXI ePart = ArrayXI::Ones(nLocElem);
  ArrayXI nPart = ArrayXI::Ones(nLocNode);

  if (myid == 0)
  {
    ArrayXI eptr = ArrayXI::LinSpaced(nElem+1, 0, Element.size());

    // Verify Inputs
    if (nparts <= 0)
      nparts = nprocs;
    if (nCommonNodes <= 0)
      nCommonNodes = pow(2, nDims-1);

    /*if (opts == NULL)                         //0 for default options
    { opts = new PetscInt; opts[0] = 0;}*/

    real_t *tpwgts_r;
    if (tpwgts == NULL)                //Vertex weight in each subdomain
    {
      tpwgts_r = new real_t[nparts];
      for (int i = 0; i < nparts; i++)
      {
        tpwgts_r[i] = 1.0/nparts;
      }
    }
    else
    {
      tpwgts_r = new real_t[nparts];
      for (int i = 0; i < nparts; i++)
      {
        tpwgts_r[i] = tpwgts[i];
      }
    }

    // Call METIS
    PetscInt METIS, edgecut;
    if (sizeof(idx_t) != sizeof(PetscInt))
    {
      cout << "WARNING, PetscInt and ParMETIS int (PetscInt) are of different " <<
              "sizes, skipping reordering with ParMETIS.\n";
      METIS = METIS_OK;
    }
    else
      METIS = METIS_PartMeshDual(&nElem, &nNode, eptr.data(),
              Element.data(), elmwgt, NULL, &nCommonNodes, &nparts,
	      tpwgts_r, opts, &edgecut, ePart.data(), nPart.data());

    delete[] tpwgts_r;
    //delete opts;

    if (METIS != METIS_OK)
    {
      std::cerr << "Error partitioning matrix! Error code: " << METIS << "\n";
      return METIS;
    }
  }

  err = ElemDist(ePart);

  return err;
}

/*****************************************************************/
/**           Get element partitioning in parallel              **/
/*****************************************************************/
PetscErrorCode Mesh::ReorderParMETIS(PetscInt nparts, PetscInt nCommonNodes, PetscScalar *tpwgts,
                                     PetscScalar *ubvec, PetscInt *opts, PetscInt ncon,
                                     PetscInt *elmwgt, PetscInt wgtflag, PetscInt numflag)
{
  PetscErrorCode err = 0;

  /// Verify Inputs
  if (nparts <= 0)
    nparts = nprocs;
  if (nCommonNodes <= 0)
    nCommonNodes = pow(2, nDims-1);

  // Local Element Descriptions - Element contains the Nodes,
  // eptr specifies where each Element starts
  short elementSize = pow(2,nDims);
  Eigen::Array<PetscInt, -1, -1> eptr =
    Eigen::Array<PetscInt, -1, 1>::LinSpaced(nLocElem+1,0,nLocElem*elementSize);
  Eigen::Array<PetscInt, -1, 1> partition = myid*Eigen::Array<PetscInt, -1, 1>::Ones(nLocElem);

  /// Initialize ParMETIS Variables
  real_t *ubvec_r;
  if (ubvec == NULL)                            //Imbalance tolerance
  {
    ubvec_r = new real_t[ncon];
    for (int i = 0; i < ncon; i++)
      ubvec_r[i] = 1.05+(double)nparts/nElem;
  }
  else
  {
    ubvec_r = new real_t[ncon];
    for (int i = 0; i < ncon; i++)
      ubvec_r[i] = ubvec[i];
  }

  idx_t *opts_r;
  if (opts == NULL)                         //0 for default options
  { opts_r = new PetscInt; opts_r[0] = 0; }
  else
  { opts_r = new PetscInt; opts_r[0] = opts[0]; }

  real_t *tpwgts_r;
  if (tpwgts == NULL)                //Vertex weight in each subdomain
  {
    tpwgts_r = new real_t[ncon*nparts];
    for (int i = 0; i < nparts; i++)
    {
      for (int j = 0; j < ncon; j++)
        tpwgts_r[i*ncon+j] = 1.0/nparts;
    }
  }
  else
  {
    tpwgts_r = new real_t[ncon*nparts];
    for (int i = 0; i < nparts; i++)
    {
      for (int j = 0; j < ncon; j++)
        tpwgts_r[i*ncon+j] = tpwgts[i*ncon+j];
    }
  }

  PetscInt edgecut;

  // Call ParMETIS
  PetscInt METIS;
  if (sizeof(idx_t) != sizeof(PetscInt))
  {
    cout << "WARNING, PetscInt and ParMETIS int (PetscInt) are of different " <<
            "sizes, skipping reordering with ParMETIS.\n";
    partition.setConstant(myid); METIS = METIS_OK;
  }
  else
    METIS = ParMETIS_V3_PartMeshKway(elDist.data(), eptr.data(),
            Element.data(), elmwgt, &wgtflag, &numflag, &ncon,
            &nCommonNodes, &nparts, tpwgts_r, ubvec_r, opts_r,
            &edgecut, partition.data(), &comm);
  
  delete[] ubvec_r;
  delete[] tpwgts_r;
  delete opts_r;

  if (METIS != METIS_OK)
  {
    std::cerr << "Error partitioning matrix! Error code: " << METIS << "\n";
    return METIS;
  }

  err = ElemDist(partition);

  return err;
}

/*****************************************************************/
/**                    Redistribute Elements                    **/
/*****************************************************************/
PetscErrorCode Mesh::ElemDist(Eigen::Array<PetscInt, -1, 1> &partition)
{
  PetscErrorCode ierr = 0;
  /// Reallocate Elements
  /// Note abbreviations: senddisp = first Element in array sent to each process
  /// sendcnt = how many Elements sent to each process - TO BE REMOVED
  /// transferSize = how many Elements each process is sending to the other processes
  /// recvcnt = how many Elements received from each process
  /// recvdsp = beginning location of buffer to receive Elements from each process
  /// elmcpy = a copy of Element reordered for continguous send buffers
  /// where = after initial sorting, the local number of each Element
  /// permute = permutation vector for filter matrix (global)
  // Initialize transfer Variables
  short elementSize = pow(2,nDims);
  ArrayXI where = EigLab::gensort(partition).cast<PetscInt>();
  ArrayXXI transferSize = ArrayXXI::Zero(nprocs,nprocs);
  ArrayXXIRM elmcpy(Element.rows(),Element.cols());
  for (PetscInt i = 0; i < partition.rows(); i++)
  {
    elmcpy.row(i) = Element.row(where(i));
    transferSize(partition(i),myid)++;
  }

  // How many Elements are transferred between each pair of processes
  MPI_Allgather(MPI_IN_PLACE, 0, MPIU_INT, transferSize.data(),
                nprocs, MPIU_INT, comm);
  Eigen::ArrayXi sendcnt = elementSize*transferSize.col(myid).cast<int>();
  Eigen::ArrayXi recvcnt = elementSize*transferSize.row(myid).cast<int>();

  // Offsets in sent messages
  Eigen::ArrayXi senddsp = Eigen::ArrayXi::Zero(nprocs);
  for (short i = 1; i < nprocs; i++)
    senddsp(i) = sendcnt(i-1) + senddsp(i-1);

  // Offsets in received messages
  Eigen::ArrayXi recvdsp = Eigen::ArrayXi::Zero(nprocs);
  for (short i = 1; i < nprocs; i++)
    recvdsp(i) = recvcnt(i-1) + recvdsp(i-1);

  // The Element transfer
  Element.resize(recvcnt.sum()/elementSize, elementSize);
  MPI_Alltoallv(elmcpy.data(), sendcnt.data(), senddsp.data(),
                MPIU_INT, Element.data(), recvcnt.data(),
                recvdsp.data(), MPIU_INT, comm);

  // Update distribution across processes
  elDist(myid+1) = Element.rows();
  nLocElem = Element.rows();
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, elDist.data()+1,
                1, MPIU_INT, comm);
  for (short i = 1; i <= nprocs; i++)
    elDist(i) += elDist(i-1);

  // Create global permutation array after sharing Elements
  // This is currently assembling a global vector on all processes and reducing
  // it.  The performance could possibly be improved by sharing the local parts
  // and then assembling after transfer, thereby reducing communications.
  // Permute = permutation vector, permute(i) = newi
  // Indices = vector indicating where this process can start assigning Elements
  //            on each process (i.e. global locations in the permute vector)
  ArrayXI permute = ArrayXI::Zero(nElem);
  ArrayXI indices = ArrayXI::Zero(nprocs);
  indices.segment(1,nprocs-1) = transferSize.block(0, 0, nprocs-1, nprocs)
                                .rowwise().sum();
  partial_sum(indices.data(), indices.data()+nprocs, indices.data());
  indices += transferSize.block(0, 0, nprocs, myid).rowwise().sum();
  int permuteStart = transferSize.block(0, 0, nprocs, myid).sum();
  for (PetscInt i = 0; i < partition.rows(); i++)
  {
    permute(where(i)+permuteStart) = indices(partition(i))++;
  }
  MPI_Allreduce(MPI_IN_PLACE, permute.data(), nElem, MPIU_INT, MPI_SUM, comm);

  return ierr;
}

/*****************************************************************/
/**                      Node Distribution                      **/
/*****************************************************************/
PetscErrorCode Mesh::NodeDist()
{
    PetscErrorCode ierr = 0;

    /// Find which Nodes each process interacts with
    Eigen::Array<short,-1,1> pckproc = Eigen::Array<short,-1,1>::Zero(nNode);
    for (PetscInt el = 0; el < Element.rows(); el++)
    {
        for (short nd = 0; nd < pow(2, nDims); nd++)
            pckproc(Element(el,nd)) = myid;
    }

    // Assign Nodes to the highest numbered processor that uses them
    MPI_Allreduce(MPI_IN_PLACE, pckproc.data(), nNode, MPI_SHORT, MPI_MAX, comm);

    /// Sort Nodes into chunks to go to each process
    Eigen::Array<short,-1,1> locpart = pckproc.segment(ndDist(myid),nLocNode);
    ArrayXI reorder = EigLab::gensort(locpart).cast<PetscInt>();
    /// Package the Nodes into a new array for sending to each process
    /// And track how many are being sent to each process
    ArrayXXSRM ndcpy(Node.rows(),Node.cols());
    Eigen::ArrayXi sendcnt = Eigen::ArrayXi::Zero(nprocs);
    for (PetscInt i = 0; i < locpart.rows(); i++)
    {
        ndcpy.row(i) = Node.row(reorder(i)).transpose();
        sendcnt(locpart(i))++;
    }

    // How much to receive from every process
    Eigen::ArrayXi recvcnt(nprocs);
    MPI_Alltoall(sendcnt.data(), 1, MPI_INT, recvcnt.data(), 1, MPI_INT, comm);

    /// Offsets in sent messages
    Eigen::ArrayXi senddsp = Eigen::ArrayXi::Zero(nprocs);
    for (short i = 1; i < nprocs; i++)
        senddsp(i) = nDims*sendcnt(i-1) + senddsp(i-1);

    /// Offsets in recieved messages
    Eigen::ArrayXi recvdsp = Eigen::ArrayXi::Zero(nprocs);
    for (short i = 1; i < nprocs; i++)
        recvdsp(i) = nDims*recvcnt(i-1) + recvdsp(i-1);

    /// The Node transfer
    Node.resize(recvcnt.sum(),nDims);
    recvcnt *= nDims; sendcnt *= nDims;
    MPI_Alltoallv(ndcpy.data(), sendcnt.data(), senddsp.data(),
                  MPI_DOUBLE, Node.data(), recvcnt.data(),
                  recvdsp.data(), MPI_DOUBLE, comm);

    /// Update the distribution of Nodes
    ndDist.setZero(nprocs+1);
    nLocNode = Node.rows();
    ndDist(myid+1) = Node.rows();
    MPI_Allreduce(MPI_IN_PLACE, ndDist.data()+1, nprocs, MPIU_INT, MPI_MAX, comm);
    for (short i = 1; i <= nprocs; i++)
        ndDist(i) += ndDist(i-1);

    /// Renumber Nodes in Element array
    reorder = EigLab::gensort(pckproc).cast<PetscInt>();
    ArrayXI invreorder(reorder.rows());
    for (PetscInt i = 0; i < reorder.rows(); i++)
      invreorder(reorder(i)) = i;
    for (PetscInt el = 0; el < Element.rows(); el++)
    {
        for (short nd = 0; nd < pow(2,nDims); nd++)
        {
            Element(el,nd) = invreorder(Element(el,nd));
        }
    }

    return ierr;
}

/*****************************************************************/
/**       Capture surrounding Nodes on other processes          **/
/*****************************************************************/
PetscErrorCode Mesh::Expand_Node()
{
  PetscErrorCode ierr = 0;

  /// List of all the Nodes the local Elements need
  ArrayXI ndlist = Eigen::Map<ArrayXI>(Element.data(),Element.size());
  EigLab::Unique(ndlist, 1);

  /// Pull out already owned Nodes
  PetscInt ind = 0, nind = 0;
  for (PetscInt i = 0; i < ndlist.rows(); i++)
  {
      if (ind == gNode.size())
      {
          ndlist.segment(nind, ndlist.rows()-i) = ndlist.segment(i, ndlist.rows()-i);
          nind += ndlist.rows()-i;
          break;
      }
      if (gNode(ind) != ndlist(i))
          ndlist(nind++) = ndlist(i);
      else
          ind++;
  }
  ndlist.conservativeResize(nind);

  /// Find where all those Nodes are and how many Nodes are neede from each process
  ArrayXI where( ndlist.rows() );
  Eigen::ArrayXi perproc = Eigen::ArrayXi::Zero( nprocs );
  short proc = 0;
  for (PetscInt i = 0; i < ndlist.rows(); i++)
  {
      while( ndlist(i) >= ndDist(proc+1) )
          proc++;
      where(i) = proc;
      perproc(proc)++;
  }

  /// Tell each process how many Nodes you need sent over
  Eigen::ArrayXi sendcnt(nprocs);
  MPI_Alltoall(perproc.data(), 1, MPI_INT, sendcnt.data(), 1, MPI_INT, comm);

  /// Offsets in recieved messages regarding which Nodes are requested
  Eigen::ArrayXi senddsp = Eigen::ArrayXi::Zero(nprocs);
  for (short i = 1; i < nprocs; i++)
      senddsp(i) = sendcnt(i-1) + senddsp(i-1);

  /// Get offsets in sent messages requesting Nodes
  Eigen::ArrayXi perprocdisp = Eigen::ArrayXi::Zero(nprocs);
  for (short i = 1; i < nprocs; i++)
      perprocdisp(i) = perprocdisp(i-1)+perproc(i-1);

  /// Send the Nodes you want to each process
  ArrayXI sendnd(sendcnt.sum());
  MPI_Alltoallv(ndlist.data(), perproc.data(), perprocdisp.data(),
                MPIU_INT, sendnd.data(), sendcnt.data(),
                senddsp.data(), MPIU_INT, comm);

  /// Pack up all the Nodes for sending
  ArrayXXSRM ndpack(sendcnt.sum(),nDims);
  for (int i = 0; i < sendcnt.sum(); i++)
      ndpack.row(i) = Node.row(sendnd(i)-ndDist(myid));

  /// Ship the Nodes
  Node.conservativeResize(nLocNode+perproc.sum(), nDims);
  perprocdisp += nLocNode; perprocdisp *= nDims;
  sendcnt *= nDims; senddsp *= nDims; perproc *= nDims;
  MPI_Alltoallv(ndpack.data(), sendcnt.data(), senddsp.data(),
                MPI_DOUBLE, Node.data(), perproc.data(),
                perprocdisp.data(), MPI_DOUBLE, comm);

  /// Update the global Node list
  gNode.conservativeResize(Node.rows());
  gNode.segment(nLocNode, ndlist.rows()) = ndlist;

  return ierr;
}

/*****************************************************************/
/**       Set up ghost communications for PEtSc vectors         **/
/*****************************************************************/
PetscErrorCode Mesh::Initialize_Vectors()
{
    PetscErrorCode ierr = 0;
    // Nodal ghost info
    ArrayXXIRM ghosts( gNode.size()-nLocNode, nDims );
    ghosts.col(0) = nDims*gNode.segment( nLocNode,gNode.size()-nLocNode );
    for (short i = 1; i < nDims; i++)
      ghosts.col(i) = ghosts.col(i-1) + 1;

    ierr = VecCreateGhost(comm, nDims*nLocNode, nDims*nNode,
                          ghosts.size(), ghosts.data(), &U); CHKERRQ(ierr);
    ierr = VecSet(U, 0.0); CHKERRQ(ierr);
    ierr = VecDuplicate(U, &F); CHKERRQ(ierr);
    ierr = VecSet(F, 0.0); CHKERRQ(ierr);

    return ierr;
}

/*****************************************************************/
/**       Convert global numberings to local numberings         **/
/*****************************************************************/
PetscErrorCode Mesh::Localize()
{
    PetscErrorCode ierr = 0;

    /// Convert Elements to local Node numbers
    PetscInt *start = gNode.data(), *finish = gNode.data()+gNode.size();
    for (PetscInt i = 0; i < Element.rows(); i++)
    {
        for (short j = 0; j < Element.cols(); j++)
            Element(i,j) = std::find(start, finish, Element(i,j)) - start;
    }

    return ierr;
}
