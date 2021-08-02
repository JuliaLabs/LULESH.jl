
function sortRegions(regReps_h, regSorted_h, ::AbstractDomain)
end

function createRegionIndexSets(nr, balance, ::AbstractDomain)
end

# access elements for comms
get_delv_xi(idx::IndexT, dom::AbstractDomain) = dom.d_delv_xi[idx]
get_delv_eta(idx::IndexT, dom::AbstractDomain) = dom.d_delv_eta[idx]
get_delv_zeta(idx::IndexT, dom::AbstractDomain) = dom.d_delv_zeta[idx]

get_x(idx::IndexT, dom::AbstractDomain) = dom.d_x[idx]
get_y(idx::IndexT, dom::AbstractDomain) = dom.d_y[idx]
get_z(idx::IndexT, dom::AbstractDomain) = dom.d_z[idx]

get_xd(idx::IndexT, dom::AbstractDomain) = dom.d_xd[idx]
get_yd(idx::IndexT, dom::AbstractDomain) = dom.d_yd[idx]
get_zd(idx::IndexT, dom::AbstractDomain) = dom.d_zd[idx]

get_fx(idx::IndexT, dom::AbstractDomain) = dom.d_fx[idx]
get_fy(idx::IndexT, dom::AbstractDomain) = dom.d_fy[idx]
get_fz(idx::IndexT, dom::AbstractDomain) = dom.d_fz[idx]

# host access
get_nodalMass(idx::IndexT, dom::AbstractDomain) = dom.h_nodalMass[idx]

colLoc(dom::AbstractDomain) = dom.m_colLoc
rowLoc(dom::AbstractDomain) = dom.m_rowLoc
planeLoc(dom::AbstractDomain) = dom.m_planeLoc
tp(dom::AbstractDomain) = dom.m_tp

function AllocateNodalPersistent(prob, domNodes)
	T = prob.devicetype{prob.floattype}
	x = T(undef, domNodes)   # coordinates
	y = T(undef, domNodes)
	z = T(undef, domNodes)

	xd = T(undef, domNodes)  # velocities
	yd = T(undef, domNodes)
	zd = T(undef, domNodes)

	xdd = T(undef, domNodes) # accelerations
	ydd = T(undef, domNodes) # accelerations
	zdd = T(undef, domNodes) # accelerations

	fx = T(undef, domNodes)   # forces
	fy = T(undef, domNodes)
	fz = T(undef, domNodes)

 	dfx = T(undef, domNodes)  # AD derivative of the forces
 	dfy = T(undef, domNodes)
 	dfz = T(undef, domNodes)

	nodalMass = T(undef, domNodes)  # mass
	return (x,y,z,xd,yd,zd,xdd,ydd,zdd,fx,fy,fz,dfx,dfy,dfz,nodalMass)
end

function AllocateElemPersistent(prob, domElems, padded_domElems)
	VDF = prob.devicetype{prob.floattype}
	VDI = prob.devicetype{IndexT}
	VDInt = prob.devicetype{Int}
	matElemlist = VDI(undef, domElems) ;  # material indexset */
	nodelist = VDI(undef, 8*padded_domElems) ;   # elemToNode connectivity */

	lxim = VDI(undef, domElems)  # elem connectivity through face g
	lxip = VDI(undef, domElems)
	letam = VDI(undef, domElems)
	letap = VDI(undef, domElems)
	lzetam = VDI(undef, domElems)
	lzetap = VDI(undef, domElems)

	elemBC = VDInt(undef, domElems)   # elem face symm/free-surf flag g

	e = VDF(undef, domElems)    # energy g
	p = VDF(undef, domElems)    # pressure g

	d_e = VDF(undef, domElems)  # AD derivative of energy E g

	q = VDF(undef, domElems)    # q g
	ql = VDF(undef, domElems)   # linear term for q g
	qq = VDF(undef, domElems)   # quadratic term for q g
	v = VDF(undef, domElems)      # relative volume g

	volo = VDF(undef, domElems)   # reference volume g
	delv = VDF(undef, domElems)   # m_vnew - m_v g
	vdov = VDF(undef, domElems)   # volume derivative over volume g

	arealg = VDF(undef, domElems)   # elem characteristic length g

	ss = VDF(undef, domElems)       # "sound speed" g

	elemMass = VDF(undef, domElems)   # mass g
	return (matElemlist, nodelist, lxim, lxip, letam, letap, lzetam, lzetap, elemBC, e, p, d_e, q, ql, qq, v, volo, delv, vdov, arealg, ss, elemMass)
end

function NewDomain(prob::LuleshProblem)
	numRanks = getNumRanks(prob.comm)
	colLoc = prob.col
	rowLoc = prob.row
	planeLoc = prob.plane
	nx = prob.nx
	tp = prob.side
	structured = prob.structured
	nr = prob.nr
	balance = prob.balance
	cost = prob.cost


	max_streams = 32;
	# domain->streams.resize(domain->max_streams);
	# TODO: CUDA stream stuff goes here
	streams = nothing

# #   for (Int_t i=0;i<domain->max_streams;i++)
# #     cudaStreamCreate(&(domain->streams[i]));

# #   cudaEventCreateWithFlags(&domain->time_constraint_computed,cudaEventDisableTiming);

#   Index_t domElems;
#   Index_t domNodes;
#   Index_t padded_domElems;

#   Vector_h<Index_t> nodelist_h;
#   Vector_h<Real_t> x_h;
#   Vector_h<Real_t> y_h;
#   Vector_h<Real_t> z_h;

  if structured
		m_tp       = tp
		m_numRanks = numRanks

	    m_colLoc   =   colLoc
	    m_rowLoc   =   rowLoc
	    m_planeLoc = planeLoc

		edgeElems = nx
		edgeNodes = edgeElems+1

		sizeX = edgeElems
		sizeY = edgeElems
		sizeZ = edgeElems

		numElem = sizeX*sizeY*sizeZ ;
		@show typeof(numElem)
		padded_numElem = PAD(numElem,32);
		@show typeof(padded_numElem)

		numNode = (sizeX+1)*(sizeY+1)*(sizeZ+1)
		padded_numNode = PAD(numNode,32);

		domElems = numElem
		domNodes = numNode
		padded_domElems = padded_numElem

		(matElemlist, nodelist, lxim, lxip, letam, letap, lzetam, lzetap,
		 elemBC, e, p, d_e, q, ql, qq, v, volo, delv, vdov, arealg, ss, elemMass) = AllocateElemPersistent(prob,domElems,padded_domElems);
		(x,y,z,xd,yd,zd,xdd,ydd,zdd,fx,fy,fz,dfx,dfy,dfz,nodalMass) = AllocateNodalPersistent(prob,domNodes);

#     domain->SetupCommBuffers(edgeNodes);

#     InitializeFields(domain);

#     domain->BuildMesh(nx, edgeNodes, edgeElems, domNodes, padded_domElems, x_h, y_h, z_h, nodelist_h);

#     domain->numSymmX = domain->numSymmY = domain->numSymmZ = 0;

#     if (domain->m_colLoc == 0)
#       domain->numSymmX = (edgeElems+1)*(edgeElems+1) ;
#     if (domain->m_rowLoc == 0)
#       domain->numSymmY = (edgeElems+1)*(edgeElems+1) ;
#     if (domain->m_planeLoc == 0)
#       domain->numSymmZ = (edgeElems+1)*(edgeElems+1) ;

#     AllocateSymmX(domain,edgeNodes*edgeNodes);
#     AllocateSymmY(domain,edgeNodes*edgeNodes);
#     AllocateSymmZ(domain,edgeNodes*edgeNodes);

#     /* set up symmetry nodesets */

#     Vector_h<Index_t> symmX_h(domain->symmX.size());
#     Vector_h<Index_t> symmY_h(domain->symmY.size());
#     Vector_h<Index_t> symmZ_h(domain->symmZ.size());

#     Int_t nidx = 0 ;
#     for (Index_t i=0; i<edgeNodes; ++i) {
#        Index_t planeInc = i*edgeNodes*edgeNodes ;
#        Index_t rowInc   = i*edgeNodes ;
#        for (Index_t j=0; j<edgeNodes; ++j) {
#          if (domain->m_planeLoc == 0) {
#            symmZ_h[nidx] = rowInc   + j ;
#          }
#          if (domain->m_rowLoc == 0) {
#            symmY_h[nidx] = planeInc + j ;
#          }
#          if (domain->m_colLoc == 0) {
#            symmX_h[nidx] = planeInc + j*edgeNodes ;
#          }
#         ++nidx ;
#        }
#     }

#     if (domain->m_planeLoc == 0)
#       domain->symmZ = symmZ_h;
#     if (domain->m_rowLoc == 0)
#       domain->symmY = symmY_h;
#     if (domain->m_colLoc == 0)
#       domain->symmX = symmX_h;

#     SetupConnectivityBC(domain, edgeElems);
#   }
  else
	error("Reading unstructured mesh is currently missing in the Julia version of LULESH.")
  end

#   /* set up node-centered indexing of elements */
#   Vector_h<Index_t> nodeElemCount_h(domNodes);

#   for (Index_t i=0; i<domNodes; ++i) {
#      nodeElemCount_h[i] = 0 ;
#   }

#   for (Index_t i=0; i<domElems; ++i) {
#      for (Index_t j=0; j < 8; ++j) {
#         ++(nodeElemCount_h[nodelist_h[j*padded_domElems+i]]);
#      }
#   }

#   Vector_h<Index_t> nodeElemStart_h(domNodes);

#   nodeElemStart_h[0] = 0;
#   for (Index_t i=1; i < domNodes; ++i) {
#      nodeElemStart_h[i] =
#         nodeElemStart_h[i-1] + nodeElemCount_h[i-1] ;
#   }

#   Vector_h<Index_t> nodeElemCornerList_h(nodeElemStart_h[domNodes-1] +
#                  nodeElemCount_h[domNodes-1] );

#   for (Index_t i=0; i < domNodes; ++i) {
#      nodeElemCount_h[i] = 0;
#   }

#   for (Index_t j=0; j < 8; ++j) {
#     for (Index_t i=0; i < domElems; ++i) {
#         Index_t m = nodelist_h[padded_domElems*j+i];
#         Index_t k = padded_domElems*j + i ;
#         Index_t offset = nodeElemStart_h[m] +
#                          nodeElemCount_h[m] ;
#         nodeElemCornerList_h[offset] = k;
#         ++(nodeElemCount_h[m]) ;
#      }
#   }

#   Index_t clSize = nodeElemStart_h[domNodes-1] +
#                    nodeElemCount_h[domNodes-1] ;
#   for (Index_t i=0; i < clSize; ++i) {
#      Index_t clv = nodeElemCornerList_h[i] ;
#      if ((clv < 0) || (clv > padded_domElems*8)) {
#           fprintf(stderr,
#    "AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!\n");
#           exit(1);
#      }
#   }

#   domain->nodeElemStart = nodeElemStart_h;
#   domain->nodeElemCount = nodeElemCount_h;
#   domain->nodeElemCornerList = nodeElemCornerList_h;

#   /* Create a material IndexSet (entire domain same material for now) */
#   Vector_h<Index_t> matElemlist_h(domElems);
#   for (Index_t i=0; i<domElems; ++i) {
#      matElemlist_h[i] = i ;
#   }
#   domain->matElemlist = matElemlist_h;

#   cudaMallocHost(&domain->dtcourant_h,sizeof(Real_t),0);
#   cudaMallocHost(&domain->dthydro_h,sizeof(Real_t),0);
#   cudaMallocHost(&domain->bad_vol_h,sizeof(Index_t),0);
#   cudaMallocHost(&domain->bad_q_h,sizeof(Index_t),0);

#   *(domain->bad_vol_h)=-1;
#   *(domain->bad_q_h)=-1;
#   *(domain->dthydro_h)=1e20;
#   *(domain->dtcourant_h)=1e20;

#   /* initialize material parameters */
#   domain->time_h      = Real_t(0.) ;
#   domain->dtfixed = Real_t(-1.0e-6) ;
#   domain->deltatimemultlb = Real_t(1.1) ;
#   domain->deltatimemultub = Real_t(1.2) ;
#   domain->stoptime  = Real_t(1.0e-2) ;
#   domain->dtmax     = Real_t(1.0e-2) ;
#   domain->cycle   = 0 ;

#   domain->e_cut = Real_t(1.0e-7) ;
#   domain->p_cut = Real_t(1.0e-7) ;
#   domain->q_cut = Real_t(1.0e-7) ;
#   domain->u_cut = Real_t(1.0e-7) ;
#   domain->v_cut = Real_t(1.0e-10) ;

#   domain->hgcoef      = Real_t(3.0) ;
#   domain->ss4o3       = Real_t(4.0)/Real_t(3.0) ;

#   domain->qstop              =  Real_t(1.0e+12) ;
#   domain->monoq_max_slope    =  Real_t(1.0) ;
#   domain->monoq_limiter_mult =  Real_t(2.0) ;
#   domain->qlc_monoq          = Real_t(0.5) ;
#   domain->qqc_monoq          = Real_t(2.0)/Real_t(3.0) ;
#   domain->qqc                = Real_t(2.0) ;

#   domain->pmin =  Real_t(0.) ;
#   domain->emin = Real_t(-1.0e+15) ;

#   domain->dvovmax =  Real_t(0.1) ;

#   domain->eosvmax =  Real_t(1.0e+9) ;
#   domain->eosvmin =  Real_t(1.0e-9) ;

#   domain->refdens =  Real_t(1.0) ;

#   /* initialize field data */
#   Vector_h<Real_t> nodalMass_h(domNodes);
#   Vector_h<Real_t> volo_h(domElems);
#   Vector_h<Real_t> elemMass_h(domElems);

#   for (Index_t i=0; i<domElems; ++i) {
#      Real_t x_local[8], y_local[8], z_local[8] ;
#      for( Index_t lnode=0 ; lnode<8 ; ++lnode )
#      {
#        Index_t gnode = nodelist_h[lnode*padded_domElems+i];
#        x_local[lnode] = x_h[gnode];
#        y_local[lnode] = y_h[gnode];
#        z_local[lnode] = z_h[gnode];
#      }

#      // volume calculations
#      Real_t volume = CalcElemVolume(x_local, y_local, z_local );
#      volo_h[i] = volume ;
#      elemMass_h[i] = volume ;
#      for (Index_t j=0; j<8; ++j) {
#         Index_t gnode = nodelist_h[j*padded_domElems+i];
#         nodalMass_h[gnode] += volume / Real_t(8.0) ;
#      }
#   }

#   domain->nodalMass = nodalMass_h;
#   domain->volo = volo_h;
#   domain->elemMass= elemMass_h;

#    /* deposit energy */
#    domain->octantCorner = 0;
#   // deposit initial energy
#   // An energy of 3.948746e+7 is correct for a problem with
#   // 45 zones along a side - we need to scale it
#   const Real_t ebase = 3.948746e+7;
#   Real_t scale = (nx*domain->m_tp)/45.0;
#   Real_t einit = ebase*scale*scale*scale;
#   //Real_t einit = ebase;
#   if (domain->m_rowLoc + domain->m_colLoc + domain->m_planeLoc == 0) {
#      // Dump into the first zone (which we know is in the corner)
#      // of the domain that sits at the origin
#        domain->e[0] = einit;
#   }

#   //set initial deltatime base on analytic CFL calculation
#   domain->deltatime_h = (.5*cbrt(domain->volo[0]))/sqrt(2*einit);

#   domain->cost = cost;
#   domain->regNumList.resize(domain->numElem) ;  // material indexset
#   domain->regElemlist.resize(domain->numElem) ;  // material indexset
#   domain->regCSR.resize(nr);
#   domain->regReps.resize(nr);
#   domain->regSorted.resize(nr);

#   // Setup region index sets. For now, these are constant sized
#   // throughout the run, but could be changed every cycle to
#   // simulate effects of ALE on the lagrange solver

#   domain->CreateRegionIndexSets(nr, balance);

	# return domain ;
end