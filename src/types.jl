@with_kw mutable struct Domain{FT}
    # Something CUDA
    max_streams::IndexT
    streams::Nothing

    # Elem-centered
    matElemlist::VD{IndexT} # material indexset
    nodelist::VD{IndexT}    # elemToNode connectivity

    lxim::VD{IndexT}         # element connectivity through face g
    lxip::VD{IndexT}
    letam::VD{IndexT}
    letap::VD{IndexT}
    lzetam::VD{IndexT}
    lzetap::VD{IndexT}

    elemBC::VD{Int}         # elem face symm/free-surf flag g

    e::VD{FT}             # energy g
    d_e::VD{FT}           # change in energy g

    p::VD{FT}             # pressure g

    q::VD{FT}             # q g
    ql::VD{FT}            # linear term for q g
    qq::VD{FT}            # quadratic term for q g

    v::VD{FT}             # relative volume g

    volo::VD{FT}          # reference volume g
    delv::VD{FT}          # m_vnew - m_v g
    vdov::VD{FT}          # volume derivative over volume g

    arealg::VD{FT}        # char length of an element g

    ss::VD{FT}            # "sound speed" g

    elemMass::VD{FT}      # mass g

    # TODO: Not making distinction between pointers to arrays down here

    vnew::VD{FT}          # new relative volume -- temporary g

    delv_xi::VD{FT}       # velocity gradient -- temporary g
    delv_eta::VD{FT}
    delv_zeta::VD{FT}

    delx_xi::VD{FT}       # coordinate gradient -- temporary g
    delx_eta::VD{FT}
    delx_zeta::VD{FT}

    dxx::VD{FT}           # principal strains -- temporary g
    dyy::VD{FT}
    dzz::VD{FT}

    # Node-centered g

    x::VD{FT}             # coordinates g
    y::VD{FT}
    z::VD{FT}

    xd::VD{FT}            # velocities g
    yd::VD{FT}
    zd::VD{FT}

    xdd::VD{FT}           # accelerations g
    ydd::VD{FT}
    zdd::VD{FT}


    fx::VD{FT}            # forces g
    fy::VD{FT}
    fz::VD{FT}

    dfx::VD{FT}          # AD of the forces g
    dfy::VD{FT}
    dfz::VD{FT}

    nodalMass::VD{FT}     # mass g
    h_nodalMass::Vector{FT}     # mass - host g

    # device pointers for comms g
    # TODO: Not sure how to store the pointers or if even necessary
    # Real_t *d_delv_xi       # velocity gradient -- temporary g
    # Real_t *d_delv_eta
    # Real_t *d_delv_zeta

    # Real_t *d_x             # coordinates g
    # Real_t *d_y
    # Real_t *d_z

    # Real_t *d_xd            # velocities g
    # Real_t *d_yd
    # Real_t *d_zd

    # Real_t *d_fx            # forces g
    # Real_t *d_fy
    # Real_t *d_fz
    # Boundary nodesets

    symmX::VD{IndexT}      # symmetry plane nodesets
    symmY::VD{IndexT}
    symmZ::VD{IndexT}

    nodeElemCount::VD{Int}
    nodeElemStart::VD{Int}
    nodeElemCornerList::VD{IndexT}

    # Parameters

    dtfixed::FT                # fixed time increment g
    deltatimemultlb::FT
    deltatimemultub::FT
    stoptime::FT               # end time for simulation g
    dtmax::FT                  # maximum allowable time increment g
    cycle::Int                 # iteration count for simulation g

    dthydro_h::FT             # hydro time constraint g
    d_dthydro_h::FT           # AD change of the hydro time constraint g
    dtcourant_h::FT           # courant time constraint g
    d_dtcourant_h::FT         # AD of the courant time constraint g
    bad_q_h::IndexT              # flag to indicate Q error g
    bad_vol_h::IndexT            # flag to indicate volume error g

    # cuda Events to indicate completion of certain kernels
    # TODO Will check later how this works with KA
    # cudaEvent_t time_constraint_computed;
    time_h::FT                # current time g
    deltatime_h::FT           # variable time increment g

    u_cut::FT                 # velocity tolerance g
    hgcoef::FT                # hourglass control g
    qstop::FT                 # excessive q indicator g
    monoq_max_slope::FT
    monoq_limiter_mult::FT
    e_cut::FT                 # energy tolerance g
    p_cut::FT                 # pressure tolerance g
    ss4o3::FT
    q_cut::FT                 # q tolerance g
    v_cut::FT                 # relative volume tolerance g
    qlc_monoq::FT             # linear term coef for q g
    qqc_monoq::FT             # quadratic term coef for q g
    qqc::FT
    eosvmax::FT
    eosvmin::FT
    pmin::FT                  # pressure floor g
    emin::FT                  # energy floor g
    dvovmax::FT               # maximum allowable volume change g
    refdens::FT               # reference density g

    m_colLoc::IndexT
    m_rowLoc::IndexT
    m_planeLoc::IndexT
    m_tp::IndexT
    sizeX::IndexT
    sizeY::IndexT
    sizeZ::IndexT
    maxPlaneSize::IndexT
    maxEdgeSize::IndexT

    numElem::IndexT
    # padded_numElem::IndexT

    numNode::IndexT
    # padded_numNode::IndexT

    numSymmX::IndexT
    numSymmY::IndexT
    numSymmZ::IndexT

    octantCorner::IndexT

    # Region information
    numReg::Int                     # number of regions (def:11)
    balance::Int                    # Load balance between regions of a domain (def: 1)
    cost::Int                       # imbalance cost (def: 1)
    regElemSize::Vector{Int}        # Size of region sets
    regCSR::VD{Int}                 # records the begining and end of each region
    regReps::VD{Int}                # records the rep number per region
    regNumList::VD{IndexT}         # Region number per domain element
    regElemlist::VD{IndexT}        # region indexset
    regSorted::VD{IndexT}          # keeps index of sorted regions


    # MPI-Related additional data

    # TODO I think we can handle this differently

    # IndexT m_numRanks;
    # IndexT& numRanks() { return m_numRanks ; }


    # # Used in setup
    m_rowMin::IndexT
    m_rowMax::IndexT
    m_colMin::IndexT
    m_colMax::IndexT
    m_planeMin::IndexT
    m_planeMax::IndexT

    # # Communication Work space
    commDataSend::VD{FT}
    commDataRecv::VD{FT}

    # Maximum number of block neighbors
    recvRequest::Vector{MPI.Request} # 6 faces + 12 edges + 8 corners
    sendRequest::Vector{MPI.Request} # 6 faces + 12 edges + 8 corners
end
