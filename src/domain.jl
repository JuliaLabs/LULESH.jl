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

# assume communication to 6 neighbors by default
m_rowMin(domain::Domain) = (domain.m_rowLoc == 0)             ? false : true
m_rowMax(domain::Domain) = (domain.m_rowLoc == domain.m_tp-1) ? false : true
m_colMin(domain::Domain) = (domain.m_colLoc == 0)             ? false : true
m_colMax(domain::Domain) = (domain.m_colLoc == domain.m_tp-1) ? false : true
m_planeMin(domain::Domain) = (domain.m_planeLoc == 0)         ? false : true
m_planeMax(domain::Domain) = (domain.m_planeLoc == domain.m_tp-1) ? false : true

# host access
get_nodalMass(idx::IndexT, dom::AbstractDomain) = dom.nodalMass[idx]

colLoc(dom::AbstractDomain) = dom.m_colLoc
rowLoc(dom::AbstractDomain) = dom.m_rowLoc
planeLoc(dom::AbstractDomain) = dom.m_planeLoc
tp(dom::AbstractDomain) = dom.m_tp

function allocateNodalPersistent!(domain, domNodes)
    fill!(resize!(domain.x, domNodes),0)   # coordinates
    fill!(resize!(domain.y, domNodes),0)
    fill!(resize!(domain.z, domNodes),0)

    fill!(resize!(domain.xd, domNodes),0)  # velocities
    fill!(resize!(domain.yd, domNodes),0)
    fill!(resize!(domain.zd, domNodes),0)

    fill!(resize!(domain.xdd, domNodes),0) # accelerations
    fill!(resize!(domain.ydd, domNodes),0) # accelerations
    fill!(resize!(domain.zdd, domNodes),0) # accelerations

    fill!(resize!(domain.fx, domNodes),0)   # forces
    fill!(resize!(domain.fy, domNodes),0)
    fill!(resize!(domain.fz, domNodes),0)

    fill!(resize!(domain.dfx, domNodes),0)  # AD derivative of the forces
    fill!(resize!(domain.dfy, domNodes),0)
    fill!(resize!(domain.dfz, domNodes),0)

    fill!(resize!(domain.nodalMass, domNodes),0)  # mass
    return nothing
end

function allocateElemPersistent!(domain, domElems)
    resize!(domain.matElemlist, domElems)  # material indexset
    resize!(domain.nodelist, 8*domElems)   # elemToNode connectivity
    fill!(domain.nodelist, 0)

    resize!(domain.lxim, domElems)  # elem connectivity through face g
    resize!(domain.lxip, domElems)
    resize!(domain.letam, domElems)
    resize!(domain.letap, domElems)
    resize!(domain.lzetam, domElems)
    resize!(domain.lzetap, domElems)

    resize!(domain.elemBC, domElems)   # elem face symm/free-surf flag g

    resize!(domain.e, domElems)    # energy g
    resize!(domain.p, domElems)    # pressure g

    resize!(domain.d_e, domElems)  # AD derivative of energy E g

    resize!(domain.q, domElems)    # q g
    resize!(domain.ql, domElems)   # linear term for q g
    resize!(domain.qq, domElems)   # quadratic term for q g
    resize!(domain.v, domElems)      # relative volume g

    resize!(domain.volo, domElems)   # reference volume g
    resize!(domain.delv, domElems)   # m_vnew - m_v g
    resize!(domain.vdov, domElems)   # volume derivative over volume g

    resize!(domain.arealg, domElems)   # elem characteristic length g

    resize!(domain.ss, domElems)       # "sound speed" g

    resize!(domain.elemMass, domElems)   # mass g
    return nothing
end

function initializeFields!(domain)
    # Basic Field Initialization

    fill!(domain.ss,0.0);
    fill!(domain.e,0.0)
    fill!(domain.p,0.0)
    fill!(domain.q,0.0)
    fill!(domain.v,1.0)

    fill!(domain.d_e,0.0)

    fill!(domain.xd,0.0)
    fill!(domain.yd,0.0)
    fill!(domain.zd,0.0)

    fill!(domain.xdd,0.0)
    fill!(domain.ydd,0.0)
    fill!(domain.zdd,0.0)

    fill!(domain.nodalMass,0.0)
end

function buildMesh!(domain, nx, edgeNodes, edgeElems, domNodes, domElems, x, y, z, nodelist)
    meshEdgeElems = domain.m_tp*nx

    resize!(x, domNodes)
    resize!(y, domNodes)
    resize!(z, domNodes)
    # initialize nodal coordinates
    # INDEXING
    nidx::IndexT = 1
    tz = 1.125*(domain.m_planeLoc*nx)/meshEdgeElems
    for plane in 1:edgeNodes
        ty = 1.125*(domain.m_rowLoc*nx)/meshEdgeElems
        for row in 1:edgeNodes
        tx = 1.125*(domain.m_colLoc*nx)/meshEdgeElems
            for col in 1:edgeNodes
                x[nidx] = tx
                y[nidx] = ty
                z[nidx] = tz
                nidx+=1
                # tx += ds ; // may accumulate roundoff...
                tx = 1.125*(domain.m_colLoc*nx+col)/meshEdgeElems
            end
        #// ty += ds ;  // may accumulate roundoff...
        ty = 1.125*(domain.m_rowLoc*nx+row)/meshEdgeElems
        end
        #// tz += ds ;  // may accumulate roundoff...
        tz = 1.125*(domain.m_planeLoc*nx+plane)/meshEdgeElems
    end

    copyto!(domain.x, x)
    copyto!(domain.y, y)
    copyto!(domain.z, z)
    resize!(nodelist, domElems*8);

    # embed hexehedral elements in nodal point lattice
    # INDEXING
    zidx::IndexT = 0
    nidx = 0
    for plane in 1:edgeElems
        for row in 1:edgeElems
            for col in 1:edgeElems
                nodelist[8*zidx+1] = nidx
                nodelist[8*zidx+2] = nidx                                   + 1
                nodelist[8*zidx+3] = nidx                       + edgeNodes + 1
                nodelist[8*zidx+4] = nidx                       + edgeNodes
                nodelist[8*zidx+5] = nidx + edgeNodes*edgeNodes
                nodelist[8*zidx+6] = nidx + edgeNodes*edgeNodes             + 1
                nodelist[8*zidx+7] = nidx + edgeNodes*edgeNodes + edgeNodes + 1
                nodelist[8*zidx+8] = nidx + edgeNodes*edgeNodes + edgeNodes
                zidx+=1
                nidx+=1
            end
        nidx+=1
        end
    nidx+=edgeNodes
    end
    copyto!(domain.nodelist, nodelist)
end

function setupConnectivityBC!(domain::Domain, edgeElems)
    domElems = domain.numElem;

    lxim = Vector{IndexT}(undef, domElems)
    lxip = Vector{IndexT}(undef, domElems)
    letam = Vector{IndexT}(undef, domElems)
    letap = Vector{IndexT}(undef, domElems)
    lzetam = Vector{IndexT}(undef, domElems)
    lzetap = Vector{IndexT}(undef, domElems)

    # set up elemement connectivity information
    lxim[1] = 0 ;
    for i in 2:domElems
       lxim[i]   = i-2
       lxip[i-1] = i-1
    end
    # MAYBE
    lxip[domElems] = domElems-1

    # INDEXING
    for i in 1:edgeElems
       letam[i] = i-1
       letap[domElems-edgeElems+i] = domElems-edgeElems+i-1
    end

    for i in (edgeElems+1):domElems
       letam[i] = i-edgeElems-1
       letap[i-edgeElems] = i-1
    end

    for i in 1:edgeElems*edgeElems
       lzetam[i] = i-1
       lzetap[domElems-edgeElems*edgeElems+i] = domElems-edgeElems*edgeElems+i-1
    end

    for i in (edgeElems*edgeElems+1):domElems
       lzetam[i] = i - edgeElems*edgeElems-1
       lzetap[i-edgeElems*edgeElems] = i-1
    end


    # set up boundary condition information
    elemBC = Vector{IndexT}(undef, domElems)
    for i in 1:domElems
        elemBC[i] = 0   # clear BCs by default
    end

    ghostIdx = [typemin(IndexT) for i in 1:6]::Vector{IndexT} # offsets to ghost locations

    pidx = domElems
    if m_planeMin(domain) != 0
        ghostIdx[1] = pidx
        pidx += domain.sizeX*domain.sizeY
    end

    if m_planeMax(domain) != 0
        ghostIdx[2] = pidx
        pidx += domain.sizeX*domain.sizeY
    end

    if m_rowMin(domain) != 0
        ghostIdx[3] = pidx
        pidx += domain.sizeX*domain.sizeZ
    end

    if m_rowMax(domain) != 0
        ghostIdx[4] = pidx
        pidx += domain.sizeX*domain.sizeZ
    end

    if m_colMin(domain) != 0
        ghostIdx[5] = pidx
        pidx += domain.sizeY*domain.sizeZ
    end

    if m_colMax(domain) != 0
        ghostIdx[6] = pidx
    end

    # symmetry plane or free surface BCs
    for i in 1:edgeElems
        planeInc = (i-1)*edgeElems*edgeElems
        rowInc   = (i-1)*edgeElems
        for j in 1:edgeElems
            if domain.m_planeLoc == 0
                elemBC[rowInc+j] |= ZETA_M_SYMM
            else
                elemBC[rowInc+j] |= ZETA_M_COMM
                lzetam[rowInc+j] = ghostIdx[1] + rowInc + (j-1)
            end

            if domain.m_planeLoc == domain.m_tp-1
                elemBC[rowInc+j+domElems-edgeElems*edgeElems] |= ZETA_P_FREE
            else
                elemBC[rowInc+j+domElems-edgeElems*edgeElems] |= ZETA_P_COMM
                lzetap[rowInc+j+domElems-edgeElems*edgeElems] = ghostIdx[2] + rowInc + (j-1)
            end

            if domain.m_rowLoc == 0
                elemBC[planeInc+j] |= ETA_M_SYMM
            else
                elemBC[planeInc+j] |= ETA_M_COMM
                letam[planeInc+j] = ghostIdx[3] + rowInc + (j-1)
            end

            if domain.m_rowLoc == domain.m_tp-1
                elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_FREE
            else
                elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_COMM
                letap[planeInc+j+edgeElems*edgeElems-edgeElems] = ghostIdx[4] +  rowInc + (j-1)
            end

            if domain.m_colLoc == 0
                elemBC[planeInc+(j-1)*edgeElems+1] |= XI_M_SYMM
            else
                elemBC[planeInc+(j-1)*edgeElems+1] |= XI_M_COMM
                lxim[planeInc+(j-1)*edgeElems+1] = ghostIdx[5] + rowInc + (j-1)
            end

            if domain.m_colLoc == domain.m_tp-1
                elemBC[planeInc+(j-1)*edgeElems+edgeElems] |= XI_P_FREE
            else
                elemBC[planeInc+(j-1)*edgeElems+edgeElems] |= XI_P_COMM
                lxip[planeInc+(j-1)*edgeElems+edgeElems] = ghostIdx[6] + rowInc + (j-1)
            end
        end
    end

    copyto!(domain.elemBC, elemBC)
    copyto!(domain.lxim, lxim)
    copyto!(domain.lxip, lxip)
    copyto!(domain.letam, letam)
    copyto!(domain.letap, letap)
    copyto!(domain.lzetam, lzetam)
    copyto!(domain.lzetap, lzetap)
end

function sortRegions(regReps::Vector{IndexT}, regSorted::Vector{IndexT}, regElemSize, numReg)
    regIndex = [v for v in 1:numReg]::Vector{IndexT}

    for i in 1:numReg-1
        for j in 1:numReg-i-1
            if regReps[j] < regReps[j+1]
                temp = regReps[j]
                regReps[j] = regReps[j+1]
                regReps[j+1] = temp

                temp = regElemSize[j]
                regElemSize[j] = regElemSize[j+1]
                regElemSize[j+1] = temp

                temp = regIndex[j]
                regIndex[j] = regIndex[j+1]
                regIndex[j+1] = temp
            end
        end
    end
    for i in 1:numReg
        regSorted[regIndex[i]] = i
    end
end

function createRegionIndexSets!(domain::Domain, nr::Int, b::Int, comm::Union{MPI.Comm, Nothing})
    domain.numReg = nr
    domain.balance = b
    @unpack_Domain domain
    myRank = getMyRank(comm)
    Random.seed!(myRank)

    regElemSize = Vector{Int}(undef, numReg)
    nextIndex::IndexT = 0

    regCSR = convert(Vector{Int}, regCSR) # records the begining and end of each region
    regReps = convert(Vector{Int}, regReps) # records the rep number per region
    regNumList = convert(Vector{IndexT}, regNumList) # Region number per domain element
    regElemlist = convert(Vector{IndexT}, regElemlist) # region indexset
    regSorted = convert(Vector{IndexT}, regSorted) # keeps index of sorted regions

    # if we only have one region just fill it
    # Fill out the regNumList with material numbers, which are always
    # the region index plus one
    if numReg == 1
        while nextIndex < numElem
            regNumList[nextIndex+1] = 1
            nextIndex+=1
        end
        regElemSize[1] = 0
    # If we have more than one region distribute the elements.
    else
        lastReg::Int = -1
        runto::IndexT = 0
        costDenominator::Int = 0
        regBinEnd = Vector{Int}(undef, numReg)
        # Determine the relative weights of all the regions.
        for i in 1:numReg
            regElemSize[i] = 0
            # INDEXING
            costDenominator += i^balance  # Total cost of all regions
            regBinEnd[i] = costDenominator  # Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
        end
        # Until all elements are assigned
        while nextIndex < numElem
            # pick the region
            regionVar = rand(Int) % costDenominator
            i = 0
            # INDEXING
            while regionVar >= regBinEnd[i+1]
                i += 1
            end
            # rotate the regions based on MPI rank.  Rotation is Rank % NumRegions
            regionNum = ((i + myRank) % numReg) + 1
            # make sure we don't pick the same region twice in a row
            while regionNum == lastReg
                regionVar = rand(Int) % costDenominator
                i = 0
                while regionVar >= regBinEnd[i+1]
                    i += 1
                end
                regionNum = ((i + myRank) % numReg) + 1
            end
            # Pick the bin size of the region and determine the number of elements.
            binSize = rand(Int) % 1000
            if binSize < 773
                elements = rand(Int) % 15 + 1
            elseif binSize < 937
                elements = rand(Int) % 16 + 16
            elseif binSize < 970
                elements = rand(Int) % 32 + 32
            elseif binSize < 974
                elements = rand(Int) % 64 + 64
            elseif binSize < 978
                elements = rand(Int) % 128 + 128
            elseif binSize < 981
                elements = rand(Int) % 256 + 256
            else
                elements = rand(Int) % 1537 + 512
            end
            runto = elements + nextIndex
            # Store the elements.  If we hit the end before we run out of elements then just stop.
            while nextIndex < runto && nextIndex < numElem
                # INDEXING
                regNumList[nextIndex+1] = regionNum
                nextIndex += 1
            end
            lastReg = regionNum
        end
    end
    # Convert regNumList to region index sets
    # First, count size of each region
    for i in 1:numElem
        # INDEXING
        r = regNumList[i] # region index == regnum-1
        regElemSize[r]+=1
    end

    # Second, allocate each region index set
    for r in 1:numReg
        if r < div(numReg, 2)
            rep = 1
        elseif r < (numReg - div((numReg+15),20))
            rep = 1 + cost;
        else
            rep = 10 * (1 + cost)
        end
        regReps[r] = rep
    end
    sortRegions(regReps, regSorted, regElemSize, numReg);

    regCSR[1] = 0;
    # Second, allocate each region index set
    for i in 2:numReg
        regCSR[i] = regCSR[i-1] + regElemSize[i-1];
    end

    # Third, fill index sets
    for i in 1:numElem
        # INDEXING
        r = regSorted[regNumList[i]] # region index == regnum-1
        regElemlist[regCSR[r]+1] = i
        regCSR[r] += 1
    end

    # Copy to device
    copyto!(regCSR, regCSR) # records the begining and end of each region
    copyto!(regReps, regReps) # records the rep number per region
    copyto!(regNumList, regNumList) # Region number per domain element
    copyto!(regElemlist, regElemlist) # region indexset
    copyto!(regSorted, regSorted) # keeps index of sorted regions
    @pack_Domain! domain
end

function Domain(prob::LuleshProblem)
    VDF = prob.devicetype{prob.floattype}
    VDI = prob.devicetype{IndexT}
    VDInt = prob.devicetype{Int}
    colLoc = prob.col
    rowLoc = prob.row
    planeLoc = prob.plane
    nx = prob.nx
    tp = prob.side
    structured = prob.structured
    nr = prob.nr
    balance = prob.balance
    cost = prob.cost
    domain = Domain{prob.floattype}(
        prob.comm,
        VDI(), VDI(),
        VDI(), VDI(), VDI(), VDI(), VDI(), VDI(),
        VDInt(),
        VDF(), VDF(),
        VDF(),
        VDF(), VDF(), VDF(),
        VDF(),
        VDF(), VDF(), VDF(), # volo
        VDF(),
        VDF(),
        VDF(), # elemMass
        VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        # FIXIT This is wrong
        VDF(),
        VDI(), VDI(), VDI(),
        VDInt(), VDInt(), VDI(),
        0.0, 0.0, 0.0, 0.0, 0.0, 0,
        0.0, 0.0, 0.0, 0.0, 0, 0,
        0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0,0,0,0,0,0,0,0,0,
        0,
        0,
        0,0,0,
        0,
        0,0,0, Vector{Int}(), VDInt(), VDInt(), VDI(), VDI(), VDI(),
        0,0,0,0,0,0,
        VDF(),VDF(),
        Vector{MPI.Request}(undef, 26), Vector{MPI.Request}(undef, 26)
    )

    nodelist = Vector{IndexT}()
    x = Vector{prob.floattype}()
    y = Vector{prob.floattype}()
    z = Vector{prob.floattype}()

    if structured
        domain.m_tp       = tp

        domain.m_colLoc   =   colLoc
        domain.m_rowLoc   =   rowLoc
        domain.m_planeLoc = planeLoc

        edgeElems = nx
        edgeNodes = edgeElems+1

        domain.sizeX = edgeElems
        domain.sizeY = edgeElems
        domain.sizeZ = edgeElems

        domain.numElem = domain.sizeX*domain.sizeY*domain.sizeZ ;
        # domain.padded_numElem = PAD(domain.numElem,32);

        domain.numNode = (domain.sizeX+1)*(domain.sizeY+1)*(domain.sizeZ+1)
        # domain.padded_numNode = PAD(domain.numNode,32);

        domElems = domain.numElem
        domNodes = domain.numNode
        # padded_domElems = domain.padded_numElem

        # Build domain object here. Not nice.

        allocateElemPersistent!(domain, domElems)
        allocateNodalPersistent!(domain, domNodes)

        setupCommBuffers!(domain, edgeNodes)

        initializeFields!(domain)

        buildMesh!(domain, nx, edgeNodes, edgeElems, domNodes, domElems, x, y, z, nodelist)

        domain.numSymmX = domain.numSymmY = domain.numSymmZ = 0

        if domain.m_colLoc == 0
            domain.numSymmX = (edgeElems+1)*(edgeElems+1)
        end
        if domain.m_rowLoc == 0
            domain.numSymmY = (edgeElems+1)*(edgeElems+1)
        end
        if domain.m_planeLoc == 0
            domain.numSymmZ = (edgeElems+1)*(edgeElems+1)
        end
        # Set up symmetry nodesets

        symmX = convert(Vector, domain.symmX)
        symmY = convert(Vector, domain.symmY)
        symmZ = convert(Vector, domain.symmZ)

        fill!(symmX, 0)
        fill!(symmY, 0)
        fill!(symmZ, 0)

        nidx = 1
        # INDEXING
        for i in 1:edgeNodes
            planeInc = (i-1)*edgeNodes*edgeNodes
            rowInc   = (i-1)*edgeNodes
            for j in 1:edgeNodes
                if domain.m_planeLoc == 0
                    symmZ[nidx] = rowInc   + j-1
                end
                if domain.m_rowLoc == 0
                    symmY[nidx] = planeInc + j-1
                end
                if domain.m_colLoc == 0
                    symmX[nidx] = planeInc + (j-1)*edgeNodes
                end
                nidx+=1
            end
        end
        if domain.m_planeLoc == 0
            domain.symmZ = symmZ
        end
        if domain.m_rowLoc == 0
            domain.symmY = symmY
        end
        if domain.m_colLoc == 0
            domain.symmX = symmX
        end

        setupConnectivityBC!(domain, edgeElems)
    else
        error("Reading unstructured mesh is currently missing in the Julia version of LULESH.")
    end
    # set up node-centered indexing of elements */
    nodeElemCount = zeros(IndexT, domNodes)
    # INDEXING
    for i in 1:domElems
        for j in 1:8
            nodeElemCount[nodelist[8*(i-1)+j]+1] += 1
        end
    end

    nodeElemStart = zeros(IndexT, domNodes)
    nodeElemStart[1] = 0
    for i in 2:domNodes
        nodeElemStart[i] = nodeElemStart[i-1] + nodeElemCount[i-1]
    end
    nodeElemCornerList = Vector{IndexT}(undef, nodeElemStart[domNodes] + nodeElemCount[domNodes] )

    nodeElemCount .= 0

    for i in 1:domElems
        for j in 1:8
            # @show i,j
            m = nodelist[8*(i-1)+j] + 1
            k = 8*(i-1) + j
            # INDEXING
            offset = nodeElemStart[m] + nodeElemCount[m]
            nodeElemCornerList[offset+1] = k
            nodeElemCount[m] += 1
        end
    end

    clSize = nodeElemStart[domNodes] + nodeElemCount[domNodes]
    for i in 1:clSize
        clv = nodeElemCornerList[i] ;
        if (clv < 0) || (clv > domElems*8)
            error("AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!")
        end
    end

    domain.nodeElemStart = convert(VDI, nodeElemStart)
    domain.nodeElemCount = convert(VDI, nodeElemCount)
    domain.nodeElemCornerList = convert(VDI, nodeElemCornerList)

    # Create a material IndexSet (entire domain same material for now)
    matElemlist = Vector{IndexT}(undef, domElems)
    for i in 1:domElems
        matElemlist[i] = i
    end
    copyto!(domain.matElemlist, matElemlist)

    domain.bad_vol = -1
    domain.bad_q = -1
    domain.dthydro = 1e20
    domain.dtcourant = 1e20

    # initialize material parameters
    domain.time      = 0.
    domain.dtfixed = -1.0e-6
    domain.deltatimemultlb = 1.1
    domain.deltatimemultub = 1.2
    domain.stoptime  = 1.0e-2
    domain.dtmax     = 1.0e-2
    domain.cycle   = 0

    domain.e_cut = 1.0e-7
    domain.p_cut = 1.0e-7
    domain.q_cut = 1.0e-7
    domain.u_cut = 1.0e-7
    domain.v_cut = 1.0e-10

    domain.hgcoef      = 3.0
    domain.ss4o3       = 4.0/3.0

    domain.qstop              =  1.0e+12
    domain.monoq_max_slope    =  1.0
    domain.monoq_limiter_mult =  2.0
    domain.qlc_monoq          = 0.5
    domain.qqc_monoq          = 2.0/3.0
    domain.qqc                = 2.0

    domain.pmin =  0.
    domain.emin = -1.0e+15

    domain.dvovmax =  0.1

    domain.eosvmax =  1.0e+9
    domain.eosvmin =  1.0e-9

    domain.refdens =  1.0

    # initialize field data
    nodalMass = Vector{prob.floattype}(undef, domNodes)
    volo = Vector{prob.floattype}(undef, domElems)
    elemMass = Vector{prob.floattype}(undef, domElems)
    fill!(nodalMass, 0)

    for i in 1:domElems
        x_local = Vector{prob.floattype}(undef, 8)
        y_local = Vector{prob.floattype}(undef, 8)
        z_local = Vector{prob.floattype}(undef, 8)
        for lnode in 1:8
            gnode = nodelist[(i-1)*8+lnode]+1
            x_local[lnode] = x[gnode]
            y_local[lnode] = y[gnode]
            z_local[lnode] = z[gnode]
        end
        # volume calculations
        volume = calcElemVolume(x_local, y_local, z_local )
        volo[i] = volume
        elemMass[i] = volume
        for j in 1:8
            gnode = nodelist[(i-1)*8+j]+1
            nodalMass[gnode] += volume / 8.0
        end
    end

    copyto!(domain.nodalMass, nodalMass)
    copyto!(domain.volo, volo)
    copyto!(domain.elemMass, elemMass)

    # deposit energy
    domain.octantCorner = 0;
    # deposit initial energy
    # An energy of 3.948746e+7 is correct for a problem with
    # 45 zones along a side - we need to scale it
    ebase = 3.948746e+7
    scale = (nx*domain.m_tp)/45.0;
    einit = ebase*scale*scale*scale;
    if domain.m_rowLoc + domain.m_colLoc + domain.m_planeLoc == 0
        # Dump into the first zone (which we know is in the corner)
        # of the domain that sits at the origin
        # TODO This only works for CUDA
        domain.e[1] = einit;
    end

    # set initial deltatime base on analytic CFL calculation
    domain.deltatime = (.5*cbrt(domain.volo[1]))/sqrt(2.0*einit)

    domain.cost = cost
    resize!(domain.regNumList, domain.numElem)  # material indexset
    resize!(domain.regElemlist, domain.numElem)  # material indexset
    resize!(domain.regCSR, nr)
    resize!(domain.regReps, nr)
    resize!(domain.regSorted, nr)

    # Setup region index sets. For now, these are constant sized
    # throughout the run, but could be changed every cycle to
    # simulate effects of ALE on the lagrange solver
    createRegionIndexSets!(domain, nr, balance, prob.comm)
    # Setup symmetry nodesets
    setupSymmetryPlanes(domain, edgeNodes)

    # Setup element connectivities
    setupElementConnectivities(domain, edgeElems)

    # Setup symmetry planes and free surface boundary arrays
    setupBoundaryConditions(domain, edgeElems)
    return domain
end

function  setupSymmetryPlanes(domain, edgeNodes)
    nidx = 1
    for i in 0:(edgeNodes-1)
        planeInc = i*edgeNodes*edgeNodes
        rowInc   = i*edgeNodes
        for j in 0:(edgeNodes-1)
            if domain.m_planeLoc == 0
                domain.symmZ[nidx] = rowInc   + j
            end
            if domain.m_rowLoc == 0
                domain.symmY[nidx] = planeInc + j
            end
            if domain.m_colLoc == 0
                domain.symmX[nidx] = planeInc + j*edgeNodes
            end
            nidx += 1
        end
    end
end

function setupElementConnectivities(domain::Domain, edgeElems)
    domain.lxim[1] = 1
    for i in 1:(domain.numElem - 1)
        domain.lxim[i+1] = i-1
        domain.lxip[i] = i
    end
    domain.lxip[domain.numElem] = domain.numElem - 1

    for i in 0:(edgeElems - 1)
        domain.letam[i+1] = i
        domain.letap[domain.numElem - edgeElems + i + 1] = i
    end

    for i in edgeElems:(domain.numElem - 1)
        domain.letam[i+1] = i - edgeElems
        domain.letap[i-edgeElems + 1] = i
    end

    for i in 0:(edgeElems*edgeElems - 1)
        domain.lzetam[i+1] = i
        domain.lzetap[domain.numElem - edgeElems*edgeElems + i + 1] = domain.numElem - edgeElems*edgeElems+i
    end

    for i in edgeElems*edgeElems:(domain.numElem - 1)
        domain.lzetam[i+1] = i - edgeElems * edgeElems
        domain.lzetap[i - edgeElems*edgeElems+1] = i
    end
end

function setupBoundaryConditions(domain::Domain, edgeElems)
  ghostIdx = Vector{IndexT}(undef, 6) # offsets to ghost locations

  # set up boundary condition information
  for i in 0:domain.numElem - 1
    domain.elemBC[i+1] = 0
  end

  for i in 1:6
    ghostIdx[i] = typemin(IndexT)
  end

  pidx = domain.numElem

  if domain.m_planeMin != 0
    ghostIdx[1] = pidx
    pidx += domain.sizeX*domain.sizeY
  end

  if m_planeMax != 0
    ghostIdx[2] = pidx
    pidx += domain.sizeX*domain.sizeY
  end

  if m_rowMin != 0
    ghostIdx[3] = pidx
    pidx += domain.sizeX*domain.sizeZ
  end

  if m_rowMax != 0
    ghostIdx[4] = pidx
    pidx += domain.sizeX*domain.sizeZ
  end

  if m_colMin != 0
    ghostIdx[5] = pidx
    pidx += domain.sizeY*domain.sizeZ
  end

  if m_colMax != 0
    ghostIdx[6] = pidx
  end


  # symmetry plane or free surface BCs

    for i in 1:edgeElems
        planeInc = (i-1)*edgeElems*edgeElems
        rowInc   = (i-1)*edgeElems
        for j in 1:edgeElems
            if domain.m_planeLoc == 0
                domain.elemBC[rowInc+j] |= ZETA_M_SYMM
            else
                domain.elemBC[rowInc+j] |= ZETA_M_COMM
                domain.lzetam[rowInc+j] = ghostIdx[1] + rowInc + (j-1)
            end

            if domain.m_planeLoc == domain.m_tp-1
                domain.elemBC[rowInc+j+domain.numElem-edgeElems*edgeElems] |= ZETA_P_FREE
            else
                domain.elemBC[rowInc+j+domain.numElem-edgeElems*edgeElems] |= ZETA_P_COMM
                domain.lzetap[rowInc+j+domain.numElem-edgeElems*edgeElems] = ghostIdx[2] + rowInc + (j-1)
            end

            if domain.m_rowLoc == 0
                domain.elemBC[planeInc+j] |= ETA_M_SYMM
            else
                domain.elemBC[planeInc+j] |= ETA_M_COMM
                domain.letam[planeInc+j] = ghostIdx[3] + rowInc + (j-1)
            end


            if domain.m_rowLoc == domain.m_tp-1
                domain.elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_FREE
            else
                domain.elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_COMM
                domain.letap[planeInc+j+edgeElems*edgeElems-edgeElems] = ghostIdx[4] +  rowInc + (j-1)
            end

            if domain.m_colLoc == 0
                domain.elemBC[planeInc+(j-1)*edgeElems+1] |= XI_M_SYMM
            else
                domain.elemBC[planeInc+(j-1)*edgeElems+1] |= XI_M_COMM
                domain.lxim[planeInc+(j-1)*edgeElems+1] = ghostIdx[5] + rowInc + (j-1)
            end

            if domain.m_colLoc == domain.m_tp-1
                domain.elemBC[planeInc+(j-1)*edgeElems+edgeElems] |= XI_P_FREE
            else
                domain.elemBC[planeInc+(j-1)*edgeElems+edgeElems] |= XI_P_COMM
                domain.lxip[planeInc+(j-1)*edgeElems+edgeElems] = ghostIdx[6] + rowInc + (j-1)
            end
        end
    end
end



function timeIncrement!(domain::Domain)
    targetdt = domain.stoptime - domain.time
    if domain.dtfixed <= 0.0 && domain.cycle != 0
        olddt = domain.deltatime

        gnewdt = typemax(Float64)
        newdt = 0.0

        if domain.dtcourant < gnewdt
            gnewdt = domain.dtcourant / 2.0
        end

        if domain.dthydro < gnewdt
            gnewdt = domain.dthydro * 2.0 / 3.0
        end
        newdt = comm_min(gnewdt, domain.comm)

        ratio = newdt / olddt
        if ratio >= 1.0
            if ratio < domain.deltatimemultlb
                newdt = olddt
            elseif ratio > domain.deltatimemultub
                newdt = olddt * domain.deltatimemultub
            end
        end

        newdt = min(newdt, domain.dtmax)

        domain.deltatime = newdt
    end

    # try to prevent very small scaling on the next cycle
    if domain.deltatime < targetdt < 4.0 * domain.deltatime / 3.0
        targetdt = 2.0 * domain.deltatime / 3.0
    end

    if targetdt < domain.deltatime
        domain.deltatime = targetdt
    end
    domain.time += domain.deltatime
    domain.cycle += 1
end

function initStressTermsForElems(domain::Domain, sigxx, sigyy, sigzz)
    # Based on FORTRAN implementation down from here
    @assert axes(sigxx) == axes(sigyy) == axes(sigzz) == axes(domain.p) == axes(domain.q)
    for i in 1:domain.numElem
        sigxx[i] = sigyy[i] = sigzz[i] = - domain.p[i] - domain.q[i]
    end
end

function sumElemFaceNormal(x0,  y0,  z0,
                           x1,  y1,  z1,
                           x2,  y2,  z2,
                           x3,  y3,  z3)
  RHALF = 0.5
  RQTR = 0.25

  bisectX0 = RHALF * (x3 + x2 - x1 - x0)
  bisectY0 = RHALF * (y3 + y2 - y1 - y0)
  bisectZ0 = RHALF * (z3 + z2 - z1 - z0)
  bisectX1 = RHALF * (x2 + x1 - x3 - x0)
  bisectY1 = RHALF * (y2 + y1 - y3 - y0)
  bisectZ1 = RHALF * (z2 + z1 - z3 - z0)
  areaX = RQTR * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1)
  areaY = RQTR * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1)
  areaZ = RQTR * (bisectX0 * bisectY1 - bisectY0 * bisectX1)

  return areaX, areaY, areaZ
end

@inline function calcElemNodeNormals(x, y, z)
    @inbounds @views begin
    pf = zeros(MMatrix{8, 3, Float64})
    pfx = view(pf, :, 1)
    pfy = view(pf, :, 2)
    pfz = view(pf, :, 3)

    # evaluate face one: nodes 1, 2, 3, 4
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[1], y[1], z[1], x[2], y[2], z[2],
                            x[3], y[3], z[3], x[4], y[4], z[4])
    pfx[1:4] .+= areaX
    pfy[1:4] .+= areaY
    pfz[1:4] .+= areaZ

    # evaluate face two: nodes 1, 5, 6, 2
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[1], y[1], z[1], x[5], y[5], z[5],
                            x[6], y[6], z[6], x[2], y[2], z[2])

    pfx[1:2] .+= areaX
    pfx[5:6] .+= areaX
    pfy[1:2] .+= areaY
    pfy[5:6] .+= areaY
    pfz[1:2] .+= areaZ
    pfz[5:6] .+= areaZ

    #evaluate face three: nodes 2, 6, 7, 3
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[2], y[2], z[2], x[6], y[6], z[6],
                            x[7], y[7], z[7], x[3], y[3], z[3])

    pfx[2:3] .+= areaX
    pfx[6:7] .+= areaX
    pfy[2:3] .+= areaY
    pfy[6:7] .+= areaY
    pfz[2:3] .+= areaZ
    pfz[6:7] .+= areaZ

    #evaluate face four: nodes 3, 7, 8, 4
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[3], y[3], z[3], x[7], y[7], z[7],
                            x[8], y[8], z[8], x[4], y[4], z[4])

    pfx[3:4] .+= areaX
    pfx[7:8] .+= areaX
    pfy[3:4] .+= areaY
    pfy[7:8] .+= areaY
    pfz[3:4] .+= areaZ
    pfz[7:8] .+= areaZ

    # evaluate face five: nodes 4, 8, 5, 1
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[4], y[4], z[4], x[8], y[8], z[8],
                            x[5], y[5], z[5], x[1], y[1], z[1])

    pfx[1]    += areaX
    pfx[4:5] .+= areaX
    pfx[8]    += areaX
    pfy[1]    += areaY
    pfy[4:5] .+= areaY
    pfy[8]    += areaY
    pfy[1]    += areaY
    pfz[4:5] .+= areaZ
    pfz[8]    += areaZ

    # evaluate face six: nodes 5, 8, 7, 6
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[5], y[5], z[5], x[8], y[8], z[8],
                            x[7], y[7], z[7], x[6], y[6], z[6])
    pfx[5:8] .+= areaX
    pfy[5:8] .+= areaY
    pfz[5:8] .+= areaZ

    return SMatrix(pf)
    end #@inbounds
end

@inline function sumElemStressesToNodeForces(B, sig_xx, sig_yy, sig_zz,  fx_elem,  fy_elem,  fz_elem, k)

  @inbounds begin
    stress_xx = sig_xx[k]
    stress_yy = sig_yy[k]
    stress_zz = sig_zz[k]

    fx = -stress_xx .* B[:, 1]
    fy = -stress_yy .* B[:, 2]
    fz = -stress_zz .* B[:, 3]

    fx_elem[(k-1)*8+1:k*8] = fx
    fy_elem[(k-1)*8+1:k*8] = fy
    fz_elem[(k-1)*8+1:k*8] = fz
  end
end


function integrateStressForElems(domain::Domain, sigxx, sigyy, sigzz, determ)
    # Based on FORTRAN implementation down from here
    # loop over all elements
    numElem8 = domain.numElem*8
    T = typeof(domain.x)
    fx_elem = T(undef, numElem8)
    fy_elem = T(undef, numElem8)
    fz_elem = T(undef, numElem8)
    # FIXIT. This has to be device type
    nodelist = domain.nodelist
    @inbounds for k in 1:domain.numElem
        x_local = collectNodal(nodelist, domain.x, (k-1)*8)
        y_local = collectNodal(nodelist, domain.y, (k-1)*8)
        z_local = collectNodal(nodelist, domain.z, (k-1)*8)
        _, detJ = calcElemShapeFunctionDerivatives(x_local, y_local, z_local)
        determ[k] = detJ
        if determ[k] <= 0.0
            error("Early Volume Error")
        end
        B = calcElemNodeNormals(x_local, y_local, z_local)
        sumElemStressesToNodeForces(B, sigxx, sigyy, sigzz, fx_elem, fy_elem, fz_elem, k)
    end

    numNode = domain.numNode

    @inbounds for gnode in 1:numNode
        count = domain.nodeElemCount[gnode]
        start = domain.nodeElemStart[gnode]
        fx = zero(eltype(fx_elem))
        fy = zero(eltype(fy_elem))
        fz = zero(eltype(fz_elem))
        @simd for i in 1:count
            elem = domain.nodeElemCornerList[start+i]
            fx = fx + fx_elem[elem]
            fy = fy + fy_elem[elem]
            fz = fz + fz_elem[elem]
        end
        domain.fx[gnode] = fx
        domain.fy[gnode] = fy
        domain.fz[gnode] = fz
    end
end

@inline function collectDomainNodesToElemNodes(domain::Domain, i)
    @inbounds begin
    nd0i = domain.nodelist[i]
    nd1i = domain.nodelist[i+1]
    nd2i = domain.nodelist[i+2]
    nd3i = domain.nodelist[i+3]
    nd4i = domain.nodelist[i+4]
    nd5i = domain.nodelist[i+5]
    nd6i = domain.nodelist[i+6]
    nd7i = domain.nodelist[i+7]

    elemX = SVector(
        domain.x[nd0i+1],
        domain.x[nd1i+1],
        domain.x[nd2i+1],
        domain.x[nd3i+1],
        domain.x[nd4i+1],
        domain.x[nd5i+1],
        domain.x[nd6i+1],
        domain.x[nd7i+1],
    )

    elemY = SVector(
        domain.y[nd0i+1],
        domain.y[nd1i+1],
        domain.y[nd2i+1],
        domain.y[nd3i+1],
        domain.y[nd4i+1],
        domain.y[nd5i+1],
        domain.y[nd6i+1],
        domain.y[nd7i+1],
    )

    elemZ = SVector(
       domain.z[nd0i+1],
       domain.z[nd1i+1],
       domain.z[nd2i+1],
       domain.z[nd3i+1],
       domain.z[nd4i+1],
       domain.z[nd5i+1],
       domain.z[nd6i+1],
       domain.z[nd7i+1],
    )

    return elemX, elemY, elemZ
    end
end

@inline function voluDer(x0, x1, x2,
                   x3, x4, x5,
                   y0, y1, y2,
                   y3, y4, y5,
                   z0, z1, z2,
                   z3, z4, z5)

    dvdx =
        (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
        (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
        (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5)

    dvdy =
        - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
        (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
        (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5)

    dvdz =
        - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
        (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
        (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5)

    twelfth = 1.0 / 12.0

    dvdx = dvdx * twelfth
    dvdy = dvdy * twelfth
    dvdz = dvdz * twelfth

    return dvdx, dvdy, dvdz
end

@inline function calcElemVolumeDerivative(x, y, z)

    dvdx1, dvdy1, dvdz1 = voluDer(
        x[2], x[3], x[4], x[5], x[6], x[8],
        y[2], y[3], y[4], y[5], y[6], y[8],
        z[2], z[3], z[4], z[5], z[6], z[8])

    dvdx2, dvdy2, dvdz2 = voluDer(
        x[1], x[2], x[3], x[8], x[5], x[7],
        y[1], y[2], y[3], y[8], y[5], y[7],
        z[1], z[2], z[3], z[8], z[5], z[7])

    dvdx3, dvdy3, dvdz3 = voluDer(
        x[4], x[1], x[2], x[7], x[8], x[6],
        y[4], y[1], y[2], y[7], y[8], y[6],
        z[4], z[1], z[2], z[7], z[8], z[6])

    dvdx4, dvdy4, dvdz4 = voluDer(
        x[3], x[4], x[1], x[6], x[7], x[5],
        y[3], y[4], y[1], y[6], y[7], y[5],
        z[3], z[4], z[1], z[6], z[7], z[5])

    dvdx5, dvdy5, dvdz5 = voluDer(
        x[8], x[7], x[6], x[1], x[4], x[2],
        y[8], y[7], y[6], y[1], y[4], y[2],
        z[8], z[7], z[6], z[1], z[4], z[2])

    dvdx6, dvdy6, dvdz6 = voluDer(
        x[5], x[8], x[7], x[2], x[1], x[3],
        y[5], y[8], y[7], y[2], y[1], y[3],
        z[5], z[8], z[7], z[2], z[1], z[3])

    dvdx7, dvdy7, dvdz7 = voluDer(
        x[6], x[5], x[8], x[3], x[2], x[4],
        y[6], y[5], y[8], y[3], y[2], y[4],
        z[6], z[5], z[8], z[3], z[2], z[4])

    dvdx8, dvdy8, dvdz8 = voluDer(
        x[7], x[6], x[5], x[4], x[3], x[1],
        y[7], y[6], y[5], y[4], y[3], y[1],
        z[7], z[6], z[5], z[4], z[3], z[1])

    dvdx = SVector(dvdx1, dvdx2, dvdx3, dvdx4, dvdx5, dvdx6, dvdx7, dvdx8)
    dvdy = SVector(dvdy1, dvdy2, dvdy3, dvdy4, dvdy5, dvdy6, dvdy7, dvdy8)
    dvdz = SVector(dvdz1, dvdz2, dvdz3, dvdz4, dvdz5, dvdz6, dvdz7, dvdz8)

    return dvdx, dvdy, dvdz
end

@inline function calcElemFBHourglassForce(xd, yd, zd,
                                  hourgam0, hourgam1,
                                  hourgam2, hourgam3,
                                  hourgam4, hourgam5,
                                  hourgam6, hourgam7,
                                  coefficient)
    i00 = 1
    i01 = 2
    i02 = 3
    i03 = 4
    @inbounds begin
    h00 = (
      hourgam0[i00] * xd[1] + hourgam1[i00] * xd[2] +
      hourgam2[i00] * xd[3] + hourgam3[i00] * xd[4] +
      hourgam4[i00] * xd[5] + hourgam5[i00] * xd[6] +
      hourgam6[i00] * xd[7] + hourgam7[i00] * xd[8]
    )

    h01 = (
      hourgam0[i01] * xd[1] + hourgam1[i01] * xd[2] +
      hourgam2[i01] * xd[3] + hourgam3[i01] * xd[4] +
      hourgam4[i01] * xd[5] + hourgam5[i01] * xd[6] +
      hourgam6[i01] * xd[7] + hourgam7[i01] * xd[8]
    )

    h02 = (
      hourgam0[i02] * xd[1] + hourgam1[i02] * xd[2] +
      hourgam2[i02] * xd[3] + hourgam3[i02] * xd[4] +
      hourgam4[i02] * xd[5] + hourgam5[i02] * xd[6] +
      hourgam6[i02] * xd[7] + hourgam7[i02] * xd[8]
    )

    h03 = (
      hourgam0[i03] * xd[1] + hourgam1[i03] * xd[2] +
      hourgam2[i03] * xd[3] + hourgam3[i03] * xd[4] +
      hourgam4[i03] * xd[5] + hourgam5[i03] * xd[6] +
      hourgam6[i03] * xd[7] + hourgam7[i03] * xd[8]
    )

    hgfx1 = coefficient *
     (hourgam0[i00] * h00 + hourgam0[i01] * h01 +
      hourgam0[i02] * h02 + hourgam0[i03] * h03)

    hgfx2 = coefficient *
     (hourgam1[i00] * h00 + hourgam1[i01] * h01 +
      hourgam1[i02] * h02 + hourgam1[i03] * h03)

    hgfx3 = coefficient *
     (hourgam2[i00] * h00 + hourgam2[i01] * h01 +
      hourgam2[i02] * h02 + hourgam2[i03] * h03)

    hgfx4 = coefficient *
     (hourgam3[i00] * h00 + hourgam3[i01] * h01 +
      hourgam3[i02] * h02 + hourgam3[i03] * h03)

    hgfx5 = coefficient *
     (hourgam4[i00] * h00 + hourgam4[i01] * h01 +
      hourgam4[i02] * h02 + hourgam4[i03] * h03)

    hgfx6 = coefficient *
     (hourgam5[i00] * h00 + hourgam5[i01] * h01 +
      hourgam5[i02] * h02 + hourgam5[i03] * h03)

    hgfx7 = coefficient *
     (hourgam6[i00] * h00 + hourgam6[i01] * h01 +
      hourgam6[i02] * h02 + hourgam6[i03] * h03)

    hgfx8 = coefficient *
     (hourgam7[i00] * h00 + hourgam7[i01] * h01 +
      hourgam7[i02] * h02 + hourgam7[i03] * h03)

    hgfx = SVector(hgfx1, hgfx2, hgfx3, hgfx4, hgfx5, hgfx6, hgfx7, hgfx8)

    h00 = (
      hourgam0[i00] * yd[1] + hourgam1[i00] * yd[2] +
      hourgam2[i00] * yd[3] + hourgam3[i00] * yd[4] +
      hourgam4[i00] * yd[5] + hourgam5[i00] * yd[6] +
      hourgam6[i00] * yd[7] + hourgam7[i00] * yd[8]
    )

    h01 = (
      hourgam0[i01] * yd[1] + hourgam1[i01] * yd[2] +
      hourgam2[i01] * yd[3] + hourgam3[i01] * yd[4] +
      hourgam4[i01] * yd[5] + hourgam5[i01] * yd[6] +
      hourgam6[i01] * yd[7] + hourgam7[i01] * yd[8]
    )

    h02 = (
      hourgam0[i02] * yd[1] + hourgam1[i02] * yd[2]+
      hourgam2[i02] * yd[3] + hourgam3[i02] * yd[4]+
      hourgam4[i02] * yd[5] + hourgam5[i02] * yd[6]+
      hourgam6[i02] * yd[7] + hourgam7[i02] * yd[8]
    )

    h03 = (
      hourgam0[i03] * yd[1] + hourgam1[i03] * yd[2] +
      hourgam2[i03] * yd[3] + hourgam3[i03] * yd[4] +
      hourgam4[i03] * yd[5] + hourgam5[i03] * yd[6] +
      hourgam6[i03] * yd[7] + hourgam7[i03] * yd[8]
    )


    hgfy1 = coefficient *
     (hourgam0[i00] * h00 + hourgam0[i01] * h01 +
      hourgam0[i02] * h02 + hourgam0[i03] * h03)

    hgfy2 = coefficient *
     (hourgam1[i00] * h00 + hourgam1[i01] * h01 +
      hourgam1[i02] * h02 + hourgam1[i03] * h03)

    hgfy3 = coefficient *
     (hourgam2[i00] * h00 + hourgam2[i01] * h01 +
      hourgam2[i02] * h02 + hourgam2[i03] * h03)

    hgfy4 = coefficient *
     (hourgam3[i00] * h00 + hourgam3[i01] * h01 +
      hourgam3[i02] * h02 + hourgam3[i03] * h03)

    hgfy5 = coefficient *
     (hourgam4[i00] * h00 + hourgam4[i01] * h01 +
      hourgam4[i02] * h02 + hourgam4[i03] * h03)

    hgfy6 = coefficient *
     (hourgam5[i00] * h00 + hourgam5[i01] * h01 +
      hourgam5[i02] * h02 + hourgam5[i03] * h03)

    hgfy7 = coefficient *
     (hourgam6[i00] * h00 + hourgam6[i01] * h01 +
      hourgam6[i02] * h02 + hourgam6[i03] * h03)

    hgfy8 = coefficient *
     (hourgam7[i00] * h00 + hourgam7[i01] * h01 +
      hourgam7[i02] * h02 + hourgam7[i03] * h03)

    hgfy = SVector(hgfy1, hgfy2, hgfy3, hgfy4, hgfy5, hgfy6, hgfy7, hgfy8)

    h00 = (
      hourgam0[i00] * zd[1] + hourgam1[i00] * zd[2] +
      hourgam2[i00] * zd[3] + hourgam3[i00] * zd[4] +
      hourgam4[i00] * zd[5] + hourgam5[i00] * zd[6] +
      hourgam6[i00] * zd[7] + hourgam7[i00] * zd[8]
    )

    h01 = (
      hourgam0[i01] * zd[1] + hourgam1[i01] * zd[2] +
      hourgam2[i01] * zd[3] + hourgam3[i01] * zd[4] +
      hourgam4[i01] * zd[5] + hourgam5[i01] * zd[6] +
      hourgam6[i01] * zd[7] + hourgam7[i01] * zd[8]
    )

    h02 =(
      hourgam0[i02] * zd[1] + hourgam1[i02] * zd[2]+
      hourgam2[i02] * zd[3] + hourgam3[i02] * zd[4]+
      hourgam4[i02] * zd[5] + hourgam5[i02] * zd[6]+
      hourgam6[i02] * zd[7] + hourgam7[i02] * zd[8]
    )

    h03 = (
      hourgam0[i03] * zd[1] + hourgam1[i03] * zd[2] +
      hourgam2[i03] * zd[3] + hourgam3[i03] * zd[4] +
      hourgam4[i03] * zd[5] + hourgam5[i03] * zd[6] +
      hourgam6[i03] * zd[7] + hourgam7[i03] * zd[8]
    )


    hgfz1 = coefficient *
     (hourgam0[i00] * h00 + hourgam0[i01] * h01 +
      hourgam0[i02] * h02 + hourgam0[i03] * h03)

    hgfz2 = coefficient *
     (hourgam1[i00] * h00 + hourgam1[i01] * h01 +
      hourgam1[i02] * h02 + hourgam1[i03] * h03)

    hgfz3 = coefficient *
     (hourgam2[i00] * h00 + hourgam2[i01] * h01 +
      hourgam2[i02] * h02 + hourgam2[i03] * h03)

    hgfz4 = coefficient *
     (hourgam3[i00] * h00 + hourgam3[i01] * h01 +
      hourgam3[i02] * h02 + hourgam3[i03] * h03)

    hgfz5 = coefficient *
     (hourgam4[i00] * h00 + hourgam4[i01] * h01 +
      hourgam4[i02] * h02 + hourgam4[i03] * h03)

    hgfz6 = coefficient *
     (hourgam5[i00] * h00 + hourgam5[i01] * h01 +
      hourgam5[i02] * h02 + hourgam5[i03] * h03)

    hgfz7 = coefficient *
     (hourgam6[i00] * h00 + hourgam6[i01] * h01 +
      hourgam6[i02] * h02 + hourgam6[i03] * h03)

    hgfz8 = coefficient *
     (hourgam7[i00] * h00 + hourgam7[i01] * h01 +
      hourgam7[i02] * h02 + hourgam7[i03] * h03)

    hgfz = SVector(hgfz1, hgfz2, hgfz3, hgfz4, hgfz5, hgfz6, hgfz7, hgfz8)

    end # inbounds

    return hgfx, hgfy, hgfz
end

@inline function calcFBHourglassForceForElems(domain, determ,
                                        x8n, y8n, z8n,
                                        dvdx, dvdy, dvdz,
                                        hourg             )

    # *************************************************
    # *
    # *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
    # *               force.
    # *
    # *************************************************

    numElem = domain.numElem
    numElem8 = numElem * 8

    fx_elem = Vector{Float64}(undef, numElem8)
    fy_elem = Vector{Float64}(undef, numElem8)
    fz_elem = Vector{Float64}(undef, numElem8)

    hourgam0 = MVector{4, Float64}(undef)
    hourgam1 = MVector{4, Float64}(undef)
    hourgam2 = MVector{4, Float64}(undef)
    hourgam3 = MVector{4, Float64}(undef)
    hourgam4 = MVector{4, Float64}(undef)
    hourgam5 = MVector{4, Float64}(undef)
    hourgam6 = MVector{4, Float64}(undef)
    hourgam7 = MVector{4, Float64}(undef)

    gamma = @SMatrix [
        1.0   1.0   1.0  -1.0
        1.0  -1.0  -1.0   1.0
       -1.0  -1.0   1.0  -1.0
       -1.0   1.0  -1.0   1.0
       -1.0  -1.0   1.0   1.0
       -1.0   1.0  -1.0  -1.0
        1.0   1.0   1.0   1.0
        1.0  -1.0  -1.0  -1.0
    ]

    # *************************************************
    # compute the hourglass modes


    @inbounds for i2 in 1:numElem

        i3=8*(i2-1)+1
        volinv= 1.0/determ[i2]

        for i1 in 1:4

            hourmodx =
                x8n[i3]   * gamma[1,i1] + x8n[i3+1] * gamma[2,i1] +
                x8n[i3+2] * gamma[3,i1] + x8n[i3+3] * gamma[4,i1] +
                x8n[i3+4] * gamma[5,i1] + x8n[i3+5] * gamma[6,i1] +
                x8n[i3+6] * gamma[7,i1] + x8n[i3+7] * gamma[8,i1]

            hourmody =
                y8n[i3]   * gamma[1,i1] + y8n[i3+1] * gamma[2,i1] +
                y8n[i3+2] * gamma[3,i1] + y8n[i3+3] * gamma[4,i1] +
                y8n[i3+4] * gamma[5,i1] + y8n[i3+5] * gamma[6,i1] +
                y8n[i3+6] * gamma[7,i1] + y8n[i3+7] * gamma[8,i1]

            hourmodz =
                z8n[i3]   * gamma[1,i1] + z8n[i3+1] * gamma[2,i1] +
                z8n[i3+2] * gamma[3,i1] + z8n[i3+3] * gamma[4,i1] +
                z8n[i3+4] * gamma[5,i1] + z8n[i3+5] * gamma[6,i1] +
                z8n[i3+6] * gamma[7,i1] + z8n[i3+7] * gamma[8,i1]

            hourgam0[i1] = gamma[1,i1] -  volinv*(dvdx[i3  ] * hourmodx +
                            dvdy[i3  ] * hourmody + dvdz[i3  ] * hourmodz )

            hourgam1[i1] = gamma[2,i1] -  volinv*(dvdx[i3+1] * hourmodx +
                            dvdy[i3+1] * hourmody + dvdz[i3+1] * hourmodz )

            hourgam2[i1] = gamma[3,i1] -  volinv*(dvdx[i3+2] * hourmodx +
                            dvdy[i3+2] * hourmody + dvdz[i3+2] * hourmodz )

            hourgam3[i1] = gamma[4,i1] -  volinv*(dvdx[i3+3] * hourmodx +
                            dvdy[i3+3] * hourmody + dvdz[i3+3] * hourmodz )

            hourgam4[i1] = gamma[5,i1] -  volinv*(dvdx[i3+4] * hourmodx +
                            dvdy[i3+4] * hourmody + dvdz[i3+4] * hourmodz )

            hourgam5[i1] = gamma[6,i1] -  volinv*(dvdx[i3+5] * hourmodx +
                            dvdy[i3+5] * hourmody + dvdz[i3+5] * hourmodz )

            hourgam6[i1] = gamma[7,i1] -  volinv*(dvdx[i3+6] * hourmodx +
                            dvdy[i3+6] * hourmody + dvdz[i3+6] * hourmodz )

            hourgam7[i1] = gamma[8,i1] -  volinv*(dvdx[i3+7] * hourmodx +
                            dvdy[i3+7] * hourmody + dvdz[i3+7] * hourmodz )
        end


        #   compute forces
        #   store forces into h arrays (force arrays)

        ss1 = domain.ss[i2]
        mass1 = domain.elemMass[i2]
        volume13=cbrt(determ[i2])

        n0si2 = domain.nodelist[(i2-1)*8+1]
        n1si2 = domain.nodelist[(i2-1)*8+2]
        n2si2 = domain.nodelist[(i2-1)*8+3]
        n3si2 = domain.nodelist[(i2-1)*8+4]
        n4si2 = domain.nodelist[(i2-1)*8+5]
        n5si2 = domain.nodelist[(i2-1)*8+6]
        n6si2 = domain.nodelist[(i2-1)*8+7]
        n7si2 = domain.nodelist[(i2-1)*8+8]

        xd1 = SVector(
            domain.xd[n0si2+1],
            domain.xd[n1si2+1],
            domain.xd[n2si2+1],
            domain.xd[n3si2+1],
            domain.xd[n4si2+1],
            domain.xd[n5si2+1],
            domain.xd[n6si2+1],
            domain.xd[n7si2+1],
        )

        yd1 = SVector(
            domain.yd[n0si2+1],
            domain.yd[n1si2+1],
            domain.yd[n2si2+1],
            domain.yd[n3si2+1],
            domain.yd[n4si2+1],
            domain.yd[n5si2+1],
            domain.yd[n6si2+1],
            domain.yd[n7si2+1],
        )

        zd1 = SVector(
            domain.zd[n0si2+1],
            domain.zd[n1si2+1],
            domain.zd[n2si2+1],
            domain.zd[n3si2+1],
            domain.zd[n4si2+1],
            domain.zd[n5si2+1],
            domain.zd[n6si2+1],
            domain.zd[n7si2+1],
        )

        coefficient = - hourg * 0.01 * ss1 * mass1 / volume13

        hgfx, hgfy, hgfz = calcElemFBHourglassForce(xd1,yd1,zd1,
                                 hourgam0,hourgam1,hourgam2,hourgam3,
                                 hourgam4,hourgam5,hourgam6,hourgam7,
                                 coefficient)

        @inbounds begin
            fx_elem[i3] = hgfx[1]
            fx_elem[i3+1] = hgfx[2]
            fx_elem[i3+2] = hgfx[3]
            fx_elem[i3+3] = hgfx[4]
            fx_elem[i3+4] = hgfx[5]
            fx_elem[i3+5] = hgfx[6]
            fx_elem[i3+6] = hgfx[7]
            fx_elem[i3+7] = hgfx[8]

            fy_elem[i3] = hgfy[1]
            fy_elem[i3+1] = hgfy[2]
            fy_elem[i3+2] = hgfy[3]
            fy_elem[i3+3] = hgfy[4]
            fy_elem[i3+4] = hgfy[5]
            fy_elem[i3+5] = hgfy[6]
            fy_elem[i3+6] = hgfy[7]
            fy_elem[i3+7] = hgfy[8]

            fz_elem[i3] = hgfz[1]
            fz_elem[i3+1] = hgfz[2]
            fz_elem[i3+2] = hgfz[3]
            fz_elem[i3+3] = hgfz[4]
            fz_elem[i3+4] = hgfz[5]
            fz_elem[i3+5] = hgfz[6]
            fz_elem[i3+6] = hgfz[7]
            fz_elem[i3+7] = hgfz[8]
        end

    # #if 0
    #     domain%m_fx(n0si2) = domain%m_fx(n0si2) + hgfx(0)
    #     domain%m_fy(n0si2) = domain%m_fy(n0si2) + hgfy(0)
    #     domain%m_fz(n0si2) = domain%m_fz(n0si2) + hgfz(0)

    #     domain%m_fx(n1si2) = domain%m_fx(n1si2) + hgfx(1)
    #     domain%m_fy(n1si2) = domain%m_fy(n1si2) + hgfy(1)
    #     domain%m_fz(n1si2) = domain%m_fz(n1si2) + hgfz(1)

    #     domain%m_fx(n2si2) = domain%m_fx(n2si2) + hgfx(2)
    #     domain%m_fy(n2si2) = domain%m_fy(n2si2) + hgfy(2)
    #     domain%m_fz(n2si2) = domain%m_fz(n2si2) + hgfz(2)

    #     domain%m_fx(n3si2) = domain%m_fx(n3si2) + hgfx(3)
    #     domain%m_fy(n3si2) = domain%m_fy(n3si2) + hgfy(3)
    #     domain%m_fz(n3si2) = domain%m_fz(n3si2) + hgfz(3)

    #     domain%m_fx(n4si2) = domain%m_fx(n4si2) + hgfx(4)
    #     domain%m_fy(n4si2) = domain%m_fy(n4si2) + hgfy(4)
    #     domain%m_fz(n4si2) = domain%m_fz(n4si2) + hgfz(4)

    #     domain%m_fx(n5si2) = domain%m_fx(n5si2) + hgfx(5)
    #     domain%m_fy(n5si2) = domain%m_fy(n5si2) + hgfy(5)
    #     domain%m_fz(n5si2) = domain%m_fz(n5si2) + hgfz(5)

    #     domain%m_fx(n6si2) = domain%m_fx(n6si2) + hgfx(6)
    #     domain%m_fy(n6si2) = domain%m_fy(n6si2) + hgfy(6)
    #     domain%m_fz(n6si2) = domain%m_fz(n6si2) + hgfz(6)

    #     domain%m_fx(n7si2) = domain%m_fx(n7si2) + hgfx(7)
    #     domain%m_fy(n7si2) = domain%m_fy(n7si2) + hgfy(7)
    #     domain%m_fz(n7si2) = domain%m_fz(n7si2) + hgfz(7)
    # #endif
    end

    numNode = domain.numNode

    @inbounds for gnode in 1:numNode
        count = domain.nodeElemCount[gnode]
        start = domain.nodeElemStart[gnode]
        fx = 0.0
        fy = 0.0
        fz = 0.0
        @simd for i in 1:count
            elem = domain.nodeElemCornerList[start+i]
            fx += fx_elem[elem]
            fy += fy_elem[elem]
            fz += fz_elem[elem]
        end
        domain.fx[gnode] += fx
        domain.fy[gnode] += fy
        domain.fz[gnode] += fz
    end
end

function calcHourglassControlForElems(domain::Domain, determ, hgcoef)
    numElem = domain.numElem
    numElem8 = numElem * 8
    dvdx = Vector{Float64}(undef, numElem8)
    dvdy = Vector{Float64}(undef, numElem8)
    dvdz = Vector{Float64}(undef, numElem8)
    x8n = Vector{Float64}(undef, numElem8)
    y8n = Vector{Float64}(undef, numElem8)
    z8n = Vector{Float64}(undef, numElem8)

    #start loop over elements
    @inbounds for i in 1:numElem
        x1, y1, z1    = collectDomainNodesToElemNodes(domain, (i-1)*8+1)
        pfx, pfy, pfz = calcElemVolumeDerivative(x1, y1, z1)

        #   load into temporary storage for FB Hour Glass control
        for ii in 1:8
            jj = 8*(i-1) + ii

            dvdx[jj] = pfx[ii]
            dvdy[jj] = pfy[ii]
            dvdz[jj] = pfz[ii]
            x8n[jj]  = x1[ii]
            y8n[jj]  = y1[ii]
            z8n[jj]  = z1[ii]
        end

        determ[i] = domain.volo[i] * domain.v[i]

        #  Do a check for negative volumes
        if domain.v[i] <= 0.0
            error("Volume Error: Volume is negative")
        end
    end

    if hgcoef > 0.0
        calcFBHourglassForceForElems(domain,determ,x8n,y8n,z8n,dvdx,dvdy,dvdz,hgcoef)
    end

    return nothing
end

function calcVolumeForceForElems(domain::Domain)
    # Based on FORTRAN implementation down from here
    hgcoef = domain.hgcoef
    numElem = domain.numElem
    VTD = typeof(domain.x)
    sigxx = VTD(undef, numElem)
    sigyy = VTD(undef, numElem)
    sigzz = VTD(undef, numElem)
    determ = VTD(undef, numElem)

    # Sum contributions to total stress tensor
    initStressTermsForElems(domain, sigxx, sigyy, sigzz)

    #   call elemlib stress integration loop to produce nodal forces from
    #   material stresses.
    integrateStressForElems(domain, sigxx, sigyy, sigzz, determ)

    # check for negative element volume and abort if found
    for i in 1:numElem
        if determ[i] <= 0.0
            error("Mid Volume Error")
        end
    end
    calcHourglassControlForElems(domain, determ, hgcoef)
end

function calcForceForNodes(domain::Domain)
    commRecv(domain, MSG_COMM_SBN, 3,
             domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
             true, false)

    domain.fx .= 0.0
    domain.fy .= 0.0
    domain.fz .= 0.0

    calcVolumeForceForElems(domain);
    fields = (domain.fx, domain.fy, domain.fz)
    commSend(domain, MSG_COMM_SBN, fields,
             domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
             true, false)
end

function calcAccelerationForNodes(domain::Domain)
    domain.xdd .= domain.fx ./ domain.nodalMass
    domain.ydd .= domain.fy ./ domain.nodalMass
    domain.zdd .= domain.fz ./ domain.nodalMass
end

function applyAccelerationBoundaryConditionsForNodes(domain::Domain)

    numNodeBC = (domain.sizeX+1)*(domain.sizeX+1)

    if length(domain.symmX) != 0
        for i in 1:numNodeBC
            domain.xdd[domain.symmX[i]+1] = 0.0
        end
    end
    if length(domain.symmY) != 0
        for i in 1:numNodeBC
            domain.ydd[domain.symmY[i]+1] = 0.0
        end
        end
    if length(domain.symmZ) != 0
        for i in 1:numNodeBC
            domain.zdd[domain.symmZ[i]+1] = 0.0
        end
    end
end

function calcVelocityForNodes(domain::Domain, dt, u_cut)


    numNode = domain.numNode

    for i in 1:numNode
        xdtmp = domain.xd[i] + domain.xdd[i] * dt
        if abs(xdtmp) < u_cut
            xdtmp = 0.0
        end
        domain.xd[i] = xdtmp

        ydtmp = domain.yd[i] + domain.ydd[i] * dt
        if abs(ydtmp) < u_cut
            ydtmp = 0.0
        end
        domain.yd[i] = ydtmp

        zdtmp = domain.zd[i] + domain.zdd[i] * dt
        if abs(zdtmp) < u_cut
            zdtmp = 0.0
        end
        domain.zd[i] = zdtmp
    end
end

function calcPositionForNodes(domain::Domain, dt)
    domain.x .= domain.x .+ domain.xd .* dt
    domain.y .= domain.y .+ domain.yd .* dt
    domain.z .= domain.z .+ domain.zd .* dt
end

function lagrangeNodal(domain::Domain)

    return nothing
end

@inline function areaFace( x0, x1, x2, x3,
                   y0, y1, y2, y3,
                   z0, z1, z2, z3  )

  fx = (x2 - x0) - (x3 - x1)
  fy = (y2 - y0) - (y3 - y1)
  fz = (z2 - z0) - (z3 - z1)
  gx = (x2 - x0) + (x3 - x1)
  gy = (y2 - y0) + (y3 - y1)
  gz = (z2 - z0) + (z3 - z1)

  area = (fx * fx + fy * fy + fz * fz) *
    (gx * gx + gy * gy + gz * gz) -
    (fx * gx + fy * gy + fz * gz) *
    (fx * gx + fy * gy + fz * gz)

  return area
end

@inline function calcElemCharacteristicLength( x, y, z, volume)

    charLength = 0.0

    a = areaFace(x[1],x[2],x[3],x[4],
                y[1],y[2],y[3],y[4],
                z[1],z[2],z[3],z[4])
    charLength = max(a,charLength)

    a = areaFace(x[5],x[6],x[7],x[8],
                y[5],y[6],y[7],y[8],
                z[5],z[6],z[7],z[8])
    charLength = max(a,charLength)

    a = areaFace(x[1],x[2],x[6],x[5],
                y[1],y[2],y[6],y[5],
                z[1],z[2],z[6],z[5])
    charLength = max(a,charLength)

    a = areaFace(x[2],x[3],x[7],x[6],
                y[2],y[3],y[7],y[6],
                z[2],z[3],z[7],z[6])
    charLength = max(a,charLength)

    a = areaFace(x[3],x[4],x[8],x[7],
                y[3],y[4],y[8],y[7],
                z[3],z[4],z[8],z[7])
    charLength = max(a,charLength)

    a = areaFace(x[4],x[1],x[5],x[8],
                y[4],y[1],y[5],y[8],
                z[4],z[1],z[5],z[8])
    charLength = max(a,charLength)

    charLength = 4.0 * volume / sqrt(charLength)

    return charLength
end

@inline function calcElemShapeFunctionDerivatives(x, y, z)
    @inbounds begin
    x0 = x[1]
    x1 = x[2]
    x2 = x[3]
    x3 = x[4]
    x4 = x[5]
    x5 = x[6]
    x6 = x[7]
    x7 = x[8]

    y0 = y[1]
    y1 = y[2]
    y2 = y[3]
    y3 = y[4]
    y4 = y[5]
    y5 = y[6]
    y6 = y[7]
    y7 = y[8]

    z0 = z[1]
    z1 = z[2]
    z2 = z[3]
    z3 = z[4]
    z4 = z[5]
    z5 = z[6]
    z6 = z[7]
    z7 = z[8]

    fjxxi = .125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) )
    fjxet = .125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) )
    fjxze = .125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) )

    fjyxi = .125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) )
    fjyet = .125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) )
    fjyze = .125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) )

    fjzxi = .125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) )
    fjzet = .125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) )
    fjzze = .125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) )

    # compute cofactors
    cjxxi =    (fjyet * fjzze) - (fjzet * fjyze)
    cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze)
    cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet)

    cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze)
    cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze)
    cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet)

    cjzxi =    (fjxet * fjyze) - (fjyet * fjxze)
    cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze)
    cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet)

    # calculate partials :
    #     this need only be done for l = 0,1,2,3   since , by symmetry ,
    #     (6,7,4,5) = - (0,1,2,3) .
    b = MMatrix{8, 3, Float64}(undef) # shape function derivatives
    b[1,1] =   -  cjxxi  -  cjxet  -  cjxze
    b[2,1] =      cjxxi  -  cjxet  -  cjxze
    b[3,1] =      cjxxi  +  cjxet  -  cjxze
    b[4,1] =   -  cjxxi  +  cjxet  -  cjxze
    b[5,1] = -b[3,1]
    b[6,1] = -b[4,1]
    b[7,1] = -b[1,1]
    b[8,1] = -b[2,1]

    b[1,2] =   -  cjyxi  -  cjyet  -  cjyze
    b[2,2] =      cjyxi  -  cjyet  -  cjyze
    b[3,2] =      cjyxi  +  cjyet  -  cjyze
    b[4,2] =   -  cjyxi  +  cjyet  -  cjyze
    b[5,2] = -b[3,2]
    b[6,2] = -b[4,2]
    b[7,2] = -b[1,2]
    b[8,2] = -b[2,2]

    b[1,3] =   -  cjzxi  -  cjzet  -  cjzze
    b[2,3] =      cjzxi  -  cjzet  -  cjzze
    b[3,3] =      cjzxi  +  cjzet  -  cjzze
    b[4,3] =   -  cjzxi  +  cjzet  -  cjzze
    b[5,3] = -b[3,3]
    b[6,3] = -b[4,3]
    b[7,3] = -b[1,3]
    b[8,3] = -b[2,3]
    end #inbounds

    # calculate jacobian determinant (volume)
    el_volume = 8.0 * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet)
    return SMatrix(b), el_volume
end

@inline function calcElemVelocityGradient( xvel, yvel, zvel, b, detJ)
    @inbounds begin
    inv_detJ = 1.0 / detJ
    pfx = b[:,1]
    pfy = b[:,2]
    pfz = b[:,3]

    d1 = inv_detJ * (   pfx[1] * (xvel[1]-xvel[7])
                        + pfx[2] * (xvel[2]-xvel[8])
                        + pfx[3] * (xvel[3]-xvel[5])
                        + pfx[4] * (xvel[4]-xvel[6]) )

    d2 = inv_detJ * (   pfy[1] * (yvel[1]-yvel[7])
                        + pfy[2] * (yvel[2]-yvel[8])
                        + pfy[3] * (yvel[3]-yvel[5])
                        + pfy[4] * (yvel[4]-yvel[6]) )

    d3 = inv_detJ * (   pfz[1] * (zvel[1]-zvel[7])
                        + pfz[2] * (zvel[2]-zvel[8])
                        + pfz[3] * (zvel[3]-zvel[5])
                        + pfz[4] * (zvel[4]-zvel[6]) )

    dyddx = inv_detJ * (  pfx[1] * (yvel[1]-yvel[7])
                        + pfx[2] * (yvel[2]-yvel[8])
                        + pfx[3] * (yvel[3]-yvel[5])
                        + pfx[4] * (yvel[4]-yvel[6]) )

    dxddy = inv_detJ * (  pfy[1] * (xvel[1]-xvel[7])
                        + pfy[2] * (xvel[2]-xvel[8])
                        + pfy[3] * (xvel[3]-xvel[5])
                        + pfy[4] * (xvel[4]-xvel[6]) )

    dzddx = inv_detJ * (  pfx[1] * (zvel[1]-zvel[7])
                        + pfx[2] * (zvel[2]-zvel[8])
                        + pfx[3] * (zvel[3]-zvel[5])
                        + pfx[4] * (zvel[4]-zvel[6]) )

    dxddz = inv_detJ * (  pfz[1] * (xvel[1]-xvel[7])
                        + pfz[2] * (xvel[2]-xvel[8])
                        + pfz[3] * (xvel[3]-xvel[5])
                        + pfz[4] * (xvel[4]-xvel[6]) )

    dzddy = inv_detJ * (  pfy[1] * (zvel[1]-zvel[7])
                        + pfy[2] * (zvel[2]-zvel[8])
                        + pfy[3] * (zvel[3]-zvel[5])
                        + pfy[4] * (zvel[4]-zvel[6]) )

    dyddz = inv_detJ * (  pfz[1] * (yvel[1]-yvel[7])
                        + pfz[2] * (yvel[2]-yvel[8])
                        + pfz[3] * (yvel[3]-yvel[5])
                        + pfz[4] * (yvel[4]-yvel[6]) )
    end #inbounds

    d6 = 0.5 * ( dxddy + dyddx )
    d5 = 0.5 * ( dxddz + dzddx )
    d4 = 0.5 * ( dzddy + dyddz )

    return SVector(d1, d2, d3, d4, d5, d6)
end

@inline function collectNodal(nodelist, src, i)
    @inbounds begin
        s1 = src[nodelist[i+1]+1]
        s2 = src[nodelist[i+2]+1]
        s3 = src[nodelist[i+3]+1]
        s4 = src[nodelist[i+4]+1]
        s5 = src[nodelist[i+5]+1]
        s6 = src[nodelist[i+6]+1]
        s7 = src[nodelist[i+7]+1]
        s8 = src[nodelist[i+8]+1]
    end

    return SVector(s1, s2, s3, s4, s5, s6, s7, s8)
end


function calcKinematicsForElems(domain::Domain, numElem, dt)

    nodelist = domain.nodelist
    # loop over all elements
    # printAllFields(domain, "$(@__FILE__):$(@__LINE__)")
    for k in 1:numElem
        # get nodal coordinates from global arrays and copy into local arrays
        @inbounds begin
            x_local = collectNodal(nodelist, domain.x, (k-1)*8)
            y_local = collectNodal(nodelist, domain.y, (k-1)*8)
            z_local = collectNodal(nodelist, domain.z, (k-1)*8)
        end

        # volume calculations
        volume = calcElemVolume(x_local, y_local, z_local)
        relativeVolume = volume / domain.volo[k]
        domain.vnew[k] = relativeVolume
        if domain.vnew[k] <= 0.0
            # @error "negative volume found in calcKinematicsForElems" rank=getMyRank(domain.comm) k volume volo=domain.volo[k] xyz=(x_local, y_local, z_local)
            error("Volume Error: calcKinematicsForElems")
        end
        domain.delv[k] = relativeVolume - domain.v[k]

        # set characteristic length
        domain.arealg[k] = calcElemCharacteristicLength(x_local, y_local,
                                                        z_local, volume)

        #  get nodal velocities from global array and copy into local arrays.
        @inbounds begin
            xd_local = collectNodal(nodelist, domain.xd, (k-1)*8)
            yd_local = collectNodal(nodelist, domain.yd, (k-1)*8)
            zd_local = collectNodal(nodelist, domain.zd, (k-1)*8)
        end

        dt2 = 0.5 * dt
        x_local = x_local .- dt2 .* xd_local
        y_local = y_local .- dt2 .* yd_local
        z_local = z_local .- dt2 .* zd_local

        B, detJ = calcElemShapeFunctionDerivatives(x_local, y_local, z_local)

        D = calcElemVelocityGradient(xd_local, yd_local, zd_local, B, detJ)

        # put velocity gradient quantities into their global arrays.
        domain.dxx[k] = D[1]
        domain.dyy[k] = D[2]
        domain.dzz[k] = D[3]
    end
end

function calcLagrangeElements(domain, delt)
    numElem = domain.numElem
    if numElem > 0
        calcKinematicsForElems(domain, numElem, delt)

        # element loop to do some stuff not included in the elemlib function.

        for k in 1:numElem
        # calc strain rate and apply as constraint (only done in FB element)
            vdov = domain.dxx[k] + domain.dyy[k] + domain.dzz[k]
            vdovthird = vdov/3.0

            # make the rate of deformation tensor deviatoric
            domain.vdov[k] = vdov
            domain.dxx[k] = domain.dxx[k] - vdovthird
            domain.dyy[k] = domain.dyy[k] - vdovthird
            domain.dzz[k] = domain.dzz[k] - vdovthird

            # See if any volumes are negative, and take appropriate action.
            if domain.vnew[k] <= 0.0
                error("Volume Error :2157")
            end
        end
    end
end

function calcMonotonicQGradientsForElems(domain::Domain)

  ptiny = 1.e-36

  numElem = domain.numElem

  for i in 1:numElem
        k = (i-1)*8
        n0 = domain.nodelist[k+1]
        n1 = domain.nodelist[k+2]
        n2 = domain.nodelist[k+3]
        n3 = domain.nodelist[k+4]
        n4 = domain.nodelist[k+5]
        n5 = domain.nodelist[k+6]
        n6 = domain.nodelist[k+7]
        n7 = domain.nodelist[k+8]

        x0 = domain.x[n0+1]
        x1 = domain.x[n1+1]
        x2 = domain.x[n2+1]
        x3 = domain.x[n3+1]
        x4 = domain.x[n4+1]
        x5 = domain.x[n5+1]
        x6 = domain.x[n6+1]
        x7 = domain.x[n7+1]

        y0 = domain.y[n0+1]
        y1 = domain.y[n1+1]
        y2 = domain.y[n2+1]
        y3 = domain.y[n3+1]
        y4 = domain.y[n4+1]
        y5 = domain.y[n5+1]
        y6 = domain.y[n6+1]
        y7 = domain.y[n7+1]

        z0 = domain.z[n0+1]
        z1 = domain.z[n1+1]
        z2 = domain.z[n2+1]
        z3 = domain.z[n3+1]
        z4 = domain.z[n4+1]
        z5 = domain.z[n5+1]
        z6 = domain.z[n6+1]
        z7 = domain.z[n7+1]

        xv0 = domain.xd[n0+1]
        xv1 = domain.xd[n1+1]
        xv2 = domain.xd[n2+1]
        xv3 = domain.xd[n3+1]
        xv4 = domain.xd[n4+1]
        xv5 = domain.xd[n5+1]
        xv6 = domain.xd[n6+1]
        xv7 = domain.xd[n7+1]

        yv0 = domain.yd[n0+1]
        yv1 = domain.yd[n1+1]
        yv2 = domain.yd[n2+1]
        yv3 = domain.yd[n3+1]
        yv4 = domain.yd[n4+1]
        yv5 = domain.yd[n5+1]
        yv6 = domain.yd[n6+1]
        yv7 = domain.yd[n7+1]

        zv0 = domain.zd[n0+1]
        zv1 = domain.zd[n1+1]
        zv2 = domain.zd[n2+1]
        zv3 = domain.zd[n3+1]
        zv4 = domain.zd[n4+1]
        zv5 = domain.zd[n5+1]
        zv6 = domain.zd[n6+1]
        zv7 = domain.zd[n7+1]

        vol = domain.volo[i]*domain.vnew[i]
        norm = 1.0 / ( vol + ptiny )

        function sum4(x1,x2,x3,x4)
            x1+x2+x3+x4
        end

        dxj = -0.25*(sum4(x0,x1,x5,x4) - sum4(x3,x2,x6,x7))
        dyj = -0.25*(sum4(y0,y1,y5,y4) - sum4(y3,y2,y6,y7))
        dzj = -0.25*(sum4(z0,z1,z5,z4) - sum4(z3,z2,z6,z7))

        dxi =  0.25*(sum4(x1,x2,x6,x5) - sum4(x0,x3,x7,x4))
        dyi =  0.25*(sum4(y1,y2,y6,y5) - sum4(y0,y3,y7,y4))
        dzi =  0.25*(sum4(z1,z2,z6,z5) - sum4(z0,z3,z7,z4))

        dxk =  0.25*(sum4(x4,x5,x6,x7) - sum4(x0,x1,x2,x3))
        dyk =  0.25*(sum4(y4,y5,y6,y7) - sum4(y0,y1,y2,y3))
        dzk =  0.25*(sum4(z4,z5,z6,z7) - sum4(z0,z1,z2,z3))

        # find delvk and delxk ( i cross j )

        ax = dyi*dzj - dzi*dyj
        ay = dzi*dxj - dxi*dzj
        az = dxi*dyj - dyi*dxj

        domain.delx_zeta[i] = vol / sqrt(ax*ax + ay*ay + az*az + ptiny)

        ax = ax * norm
        ay = ay * norm
        az = az * norm

        dxv = 0.25*(sum4(xv4,xv5,xv6,xv7) - sum4(xv0,xv1,xv2,xv3))
        dyv = 0.25*(sum4(yv4,yv5,yv6,yv7) - sum4(yv0,yv1,yv2,yv3))
        dzv = 0.25*(sum4(zv4,zv5,zv6,zv7) - sum4(zv0,zv1,zv2,zv3))

        domain.delv_zeta[i] = ax*dxv + ay*dyv + az*dzv

        # find delxi and delvi ( j cross k )

        ax = dyj*dzk - dzj*dyk
        ay = dzj*dxk - dxj*dzk
        az = dxj*dyk - dyj*dxk

        domain.delx_xi[i] = vol / sqrt(ax*ax + ay*ay + az*az + ptiny)

        ax = ax * norm
        ay = ay * norm
        az = az * norm

        dxv = 0.25*(sum4(xv1,xv2,xv6,xv5) - sum4(xv0,xv3,xv7,xv4))
        dyv = 0.25*(sum4(yv1,yv2,yv6,yv5) - sum4(yv0,yv3,yv7,yv4))
        dzv = 0.25*(sum4(zv1,zv2,zv6,zv5) - sum4(zv0,zv3,zv7,zv4))

        domain.delv_xi[i] = ax*dxv + ay*dyv + az*dzv ;

        # find delxj and delvj ( k cross i )

        ax = dyk*dzi - dzk*dyi ;
        ay = dzk*dxi - dxk*dzi ;
        az = dxk*dyi - dyk*dxi ;

        domain.delx_eta[i] = vol / sqrt(ax*ax + ay*ay + az*az + ptiny)

        ax = ax * norm
        ay = ay * norm
        az = az * norm

        dxv = -0.25*(sum4(xv0,xv1,xv5,xv4) - sum4(xv3,xv2,xv6,xv7))
        dyv = -0.25*(sum4(yv0,yv1,yv5,yv4) - sum4(yv3,yv2,yv6,yv7))
        dzv = -0.25*(sum4(zv0,zv1,zv5,zv4) - sum4(zv3,zv2,zv6,zv7))

        domain.delv_eta[i] = ax*dxv + ay*dyv + az*dzv ;
    end
    return nothing
end

function calcMonotonicQRegionForElems(domain::Domain, qlc_monoq, qqc_monoq,
                                        monoq_limiter_mult, monoq_max_slope,
                                        ptiny,
                                        elength
                                    )
    for ielem in 1:elength
        i = domain.matElemlist[ielem]
        bcMask = domain.elemBC[i]

    #   phixi
        norm = 1.0 / ( domain.delv_xi[i] + ptiny )

        case = bcMask & XI_M
        if case == 0 || case == XI_M_COMM
            # MAYBE
            delvm = domain.delv_xi[domain.lxim[i]+1]
        elseif case == XI_M_SYMM
            delvm = domain.delv_xi[i]
        elseif case == XI_M_FREE
            delvm = 0.0
        else
            error("Error")
            delvm = 0.0
        end

        case = bcMask & XI_P
        if case == 0 || case == XI_P_COMM
            delvp = domain.delv_xi[domain.lxip[i]]
        elseif case == XI_P_SYMM
            delvp = domain.delv_xi[i]
        elseif case == XI_P_FREE
            delvp = 0.0
        else
            error("Error")
            delvp = 0.0
        end

        delvm = delvm * norm
        delvp = delvp * norm

        phixi = 0.5 * ( delvm + delvp )

        delvm = delvm * monoq_limiter_mult
        delvp = delvp * monoq_limiter_mult

        if  delvm < phixi
            phixi = delvm
        end
        if  delvp < phixi
            phixi = delvp
        end
        if  phixi < 0.0
            phixi = 0.0
        end
        if  phixi > monoq_max_slope
            phixi = monoq_max_slope
        end


    #   phieta
        norm = 1.0 / ( domain.delv_eta[i] + ptiny )

        case = bcMask & ETA_M
        if case == 0 || case == ETA_M_COMM
            delvm = domain.delv_eta[domain.letam[i]+1]
        elseif case == ETA_M_SYMM
            delvm = domain.delv_eta[i]
        elseif case == ETA_M_FREE
            delvm = 0.0
        else
            delvm = 0.0
            error("Error")
        end

        case = bcMask & ETA_P
        if case == 0 || case == ETA_P_COMM
            delvp = domain.delv_eta[domain.letap[i]+1]
        elseif case == ETA_P_SYMM
            delvp = domain.delv_eta[i]
        elseif case == ETA_P_FREE
            delvp = 0.0
        else
            delvp = 0.0
            error("Error")
        end

        delvm = delvm * norm
        delvp = delvp * norm

        phieta = 0.5 * ( delvm + delvp )

        delvm = delvm * monoq_limiter_mult
        delvp = delvp * monoq_limiter_mult

        if delvm  < phieta
            phieta = delvm
        end
        if delvp  < phieta
            phieta = delvp
        end
        if phieta < 0.0
            phieta = 0.0
        end
        if phieta > monoq_max_slope
            phieta = monoq_max_slope
        end

    #   phizeta
        norm = 1.0 / ( domain.delv_zeta[i] + ptiny )

        case = bcMask & ZETA_M
        if case == 0 || case == ZETA_M_COMM
            delvm = domain.delv_zeta[domain.lzetam[i]+1]
        elseif case == ZETA_M_SYMM
            delvm = domain.delv_zeta[i]
        elseif case == ZETA_M_FREE
            delvm = 0.0
        else
            delvm = 0.0
            error("Error")
        end

        case = bcMask & ZETA_P
        if case == 0 || case == ZETA_P_COMM
            delvp = domain.delv_zeta[domain.lzetap[i]+1]
        elseif case == ZETA_P_SYMM
            delvp = domain.delv_zeta[i]
        elseif case == ZETA_P_FREE
            delvp = 0.0
        else
            delvp = 0.0
            error("Error")
        end

        delvm = delvm * norm
        delvp = delvp * norm

        phizeta = 0.5 * ( delvm + delvp )

        delvm = delvm * monoq_limiter_mult
        delvp = delvp * monoq_limiter_mult

        if delvm < phizeta
            phizeta = delvm
        end
        if delvp < phizeta
            phizeta = delvp
        end
        if phizeta < 0.0
            phizeta = 0.0
        end
        if phizeta > monoq_max_slope
            phizeta = monoq_max_slope
        end

    #   Remove length scale

        if domain.vdov[i] > 0.0
            qlin  = 0.0
            qquad = 0.0
        else
            delvxxi   = domain.delv_xi[i]   * domain.delx_xi[i]
            delvxeta  = domain.delv_eta[i]  * domain.delx_eta[i]
            delvxzeta = domain.delv_zeta[i] * domain.delx_zeta[i]

            if delvxxi   > 0.0
                delvxxi   = 0.0
            end
            if delvxeta  > 0.0
                delvxeta  = 0.0
            end
            if delvxzeta > 0.0
                delvxzeta = 0.0
            end

            rho = domain.elemMass[i] / (domain.volo[i] * domain.vnew[i])

            qlin = -qlc_monoq * rho *
                    (  delvxxi   * (1.0 - phixi)  +
                        delvxeta  * (1.0 - phieta) +
                        delvxzeta * (1.0 - phizeta)  )

            qquad = qqc_monoq * rho *
                    (  delvxxi*delvxxi     * (1.0 - phixi*phixi)   +
                        delvxeta*delvxeta   * (1.0 - phieta*phieta) +
                        delvxzeta*delvxzeta * (1.0 - phizeta*phizeta)  )
        end

        domain.qq[i] = qquad
        domain.ql[i] = qlin
    end
end


function calcMonotonicQForElems(domain::Domain)

    ptiny = 1e-36
    #
    # initialize parameters
    #
    monoq_max_slope    = domain.monoq_max_slope
    monoq_limiter_mult = domain.monoq_limiter_mult

    #
    # calculate the monotonic q for pure regions
    #
    elength = domain.numElem
    if elength > 0
        qlc_monoq = domain.qlc_monoq
        qqc_monoq = domain.qqc_monoq
        calcMonotonicQRegionForElems(domain, qlc_monoq, qqc_monoq,
                                        monoq_limiter_mult,
                                        monoq_max_slope,
                                        ptiny, elength )
    end
end



function calcQForElems(domain::Domain)


    qstop = domain.qstop
    numElem = domain.numElem

    # MONOTONIC Q option
    commRecv(domain, MSG_MONOQ, 3,
             domain.sizeX, domain.sizeY, domain.sizeZ,
             true, true)

    # Calculate velocity gradients
    calcMonotonicQGradientsForElems(domain)

    # Transfer veloctiy gradients in the first order elements
    # problem->commElements->Transfer(CommElements::monoQ)

    fields = (domain.delv_xi, domain.delv_eta, domain.delv_zeta)
    commSend(domain, MSG_MONOQ, fields,
             domain.sizeX, domain.sizeY, domain.sizeZ,
             true, true)
    commMonoQ(domain)

    calcMonotonicQForElems(domain)

    # Don't allow excessive artificial viscosity
    if numElem != 0
        idx = -1
        for i in 1:numElem
            if domain.q[i] > qstop
                idx = i
                break
            end
        end

        if idx >= 0
            if domain.comm !== nothing
                MPI.Abort(MPI.COMM_WORLD, 1)
            end
            error("QStopError")
        end
    end
end

function calcPressureForElems(domain::Domain, p_new, bvc,
                                 pbvc, e_old,
                                 compression, vnewc,
                                 pmin,
                                 p_cut,eosvmax,
                                 length              )

    c1s = 2.0/3.0

    for i in 1:length
        bvc[i] = c1s * (compression[i] + 1.0)
        pbvc[i] = c1s
    end

    for i in 1:length
        p_new[i] = bvc[i] * e_old[i]

        if abs(p_new[i]) < p_cut
            p_new[i] = 0.0
        end

        if vnewc[i] >= eosvmax # impossible condition here?
            p_new[i] = 0.0
        end

        if p_new[i] < pmin
            p_new[i] = pmin
        end
    end
end


function calcEnergyForElems(domain::Domain, p_new,  e_new,  q_new,
                                bvc,  pbvc,
                                p_old,  e_old,  q_old,
                                compression,  compHalfStep,
                                vnewc,  work,  delvc,  pmin,
                                p_cut,   e_cut,  q_cut,  emin,
                                qq,  ql,
                                rho0,
                                eosvmax,
                                length                          )

    TINY1 = 0.111111e-36
    TINY3 = 0.333333e-18
    SIXTH = 1.0 / 6.0


    pHalfStep = Vector{Float64}(undef, length)

    for i in 1:length
        e_new[i] = e_old[i] - 0.5 * delvc[i] * (p_old[i] + q_old[i]) + 0.5 * work[i]

        if e_new[i]  < emin
            e_new[i] = emin
        end
    end

    calcPressureForElems(domain, pHalfStep, bvc, pbvc, e_new, compHalfStep,
                            vnewc, pmin, p_cut, eosvmax, length)
    for i in 1:length
        vhalf = 1.0 / (1.0 + compHalfStep[i])

        if  delvc[i] > 0.0
        #      q_new(i) /* = qq(i) = ql(i) */ = Real_t(0.) ;
            q_new[i] = 0.0
        else
            ssc = (( pbvc[i] * e_new[i]
                + vhalf * vhalf * bvc[i] * pHalfStep[i] ) / rho0)
            if ssc <= TINY1
                ssc = TINY3
            else
                ssc = sqrt(ssc)
            end

            q_new[i] = (ssc*ql[i] + qq[i])
        end

        e_new[i] = (e_new[i] + 0.5 * delvc[i] * (  3.0*(p_old[i]     + q_old[i])
            - 4.0*(pHalfStep[i] + q_new[i])))
    end

    for i in 1:length
        e_new[i] = e_new[i] + 0.5 * work[i]
        if abs(e_new[i]) < e_cut
            e_new[i] = 0.0
        end
        if e_new[i]  < emin
            e_new[i] = emin
        end
    end

    calcPressureForElems(domain, p_new, bvc, pbvc, e_new, compression,
                            vnewc, pmin, p_cut, eosvmax, length)

    for i in 1:length
        if delvc[i] > 0.0
            q_tilde = 0.0
        else
            ssc = ( pbvc[i] * e_new[i]
                + vnewc[i] * vnewc[i] * bvc[i] * p_new[i] ) / rho0

            if ssc <= TINY1
                ssc = TINY3
            else
                ssc = sqrt(ssc)
            end

            q_tilde = (ssc*ql[i] + qq[i])
        end

        e_new[i] = (e_new[i] - (  7.0*(p_old[i]     + q_old[i])
                            -    8.0*(pHalfStep[i] + q_new[i])
                            + (p_new[i] + q_tilde)) * delvc[i]*SIXTH)

        if abs(e_new[i]) < e_cut
            e_new[i] = 0.0
        end
        if e_new[i]  < emin
            e_new[i] = emin
        end
    end

    calcPressureForElems(domain, p_new, bvc, pbvc, e_new, compression,
                            vnewc, pmin, p_cut, eosvmax, length)

    for i in 1:length

        if delvc[i] <= 0.0
            ssc = (( pbvc[i] * e_new[i]
                + vnewc[i] * vnewc[i] * bvc[i] * p_new[i] ) / rho0)

            if ssc <= TINY1
                ssc = TINY3
            else
                ssc = sqrt(ssc)
            end

            q_new[i] = (ssc*ql[i] + qq[i])

            if abs(q_new[i]) < q_cut
                q_new[i] = 0.0
            end
        end
    end
end

function calcSoundSpeedForElems(domain::Domain, vnewc,  rho0, enewc,
                                  pnewc, pbvc,
                                  bvc, ss4o3, nz       )
    TINY1 = 0.111111e-36
    TINY3 = 0.333333e-18

    for i in 1:nz
        iz = domain.matElemlist[i]
        ssTmp = ((pbvc[i] * enewc[i] + vnewc[i] * vnewc[i] *
                            bvc[i] * pnewc[i]) / rho0)
        if ssTmp <= TINY1
            ssTmp = TINY3
        else
            ssTmp = sqrt(ssTmp)
        end
        domain.ss[iz] = ssTmp
    end
end


function evalEOSForElems(domain::Domain, vnewc, length)

    e_cut = domain.e_cut
    p_cut = domain.p_cut
    ss4o3 = domain.ss4o3
    q_cut = domain.q_cut

    eosvmax = domain.eosvmax
    eosvmin = domain.eosvmin
    pmin    = domain.pmin
    emin    = domain.emin
    rho0    = domain.refdens

    e_old = Vector{Float64}(undef, length)
    delvc = Vector{Float64}(undef, length)
    p_old = Vector{Float64}(undef, length)
    q_old = Vector{Float64}(undef, length)
    compression = Vector{Float64}(undef, length)
    compHalfStep = Vector{Float64}(undef, length)
    qq = Vector{Float64}(undef, length)
    ql = Vector{Float64}(undef, length)
    work = Vector{Float64}(undef, length)
    p_new = Vector{Float64}(undef, length)
    e_new = Vector{Float64}(undef, length)
    q_new = Vector{Float64}(undef, length)
    bvc = Vector{Float64}(undef, length)
    pbvc = Vector{Float64}(undef, length)

    # compress data, minimal set
    for i in 1:length
        zidx = domain.matElemlist[i]
        e_old[i] = domain.e[zidx]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        delvc[i] = domain.delv[zidx]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        p_old[i] = domain.p[zidx]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        q_old[i] = domain.q[zidx]
    end

    for i in 1:length
        compression[i] = 1.0 / vnewc[i] - 1.0
        vchalf = vnewc[i] - delvc[i] * 0.5
        compHalfStep[i] = (1.0 / vchalf) - 1.0
    end

    # Check for v > eosvmax or v < eosvmin
    if eosvmin != 0.0
        for i in 1:length
            if vnewc[i] <= eosvmin  # impossible due to calling func?
                compHalfStep[i] = compression[i]
            end
        end
    end

    if eosvmax != 0.0
        for i in 1:length
            if vnewc[i] >= eosvmax  # impossible due to calling func?
                p_old[i]        = 0.0
                compression[i]  = 0.0
                compHalfStep[i] = 0.0
            end
        end
    end
    for i in 1:length
        zidx = domain.matElemlist[i]
        qq[i] = domain.qq[zidx]
        ql[i] = domain.ql[zidx]
        work[i] = 0.0
    end

    calcEnergyForElems(domain, p_new, e_new, q_new, bvc, pbvc,
                            p_old, e_old,  q_old, compression,
                            compHalfStep, vnewc, work,  delvc, pmin,
                            p_cut, e_cut, q_cut, emin,
                            qq, ql, rho0, eosvmax, length)


    for i in 1:length
        zidx = domain.matElemlist[i]
        domain.p[zidx] = p_new[i]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        domain.e[zidx] = e_new[i]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        domain.q[zidx] = q_new[i]
    end


    calcSoundSpeedForElems(domain, vnewc, rho0, e_new, p_new,
                                pbvc, bvc, ss4o3, length)
end


function applyMaterialPropertiesForElems(domain::Domain)

    length = domain.numElem

    if length != 0
        #  Expose all of the variables needed for material evaluation
        eosvmin = domain.eosvmin
        eosvmax = domain.eosvmax
        vnewc = Vector{Float64}(undef, length)

        for i in 1:length
            zn = domain.matElemlist[i]
            vnewc[i] = domain.vnew[zn]
        end

        if eosvmin != 0.0
            for i in 1:length
                if vnewc[i] < eosvmin
                    vnewc[i] = eosvmin
                end
            end
        end

        if eosvmax != 0.0
            for i in 1:length
                if vnewc[i] > eosvmax
                    vnewc[i] = eosvmax
                end
            end
        end

        for i in 1:length
            zn = domain.matElemlist[i]
            vc = domain.v[zn]
            if eosvmin != 0.0
                if vc < eosvmin
                    vc = eosvmin
                end
            end
            if eosvmax != 0.0
                if vc > eosvmax
                    vc = eosvmax
                end
            end
            if vc <= 0.0
                error("Volume Error :2887")
            end
        end
        evalEOSForElems(domain::Domain, vnewc, length)
    end
end

function updateVolumesForElems(domain::Domain)
    numElem = domain.numElem

    if numElem != 0
        v_cut = domain.v_cut

        for i in 1:numElem
            tmpV = domain.vnew[i]

            if abs(tmpV - 1.0) < v_cut
                tmpV = 1.0
            end
            domain.v[i] = tmpV
            if tmpV <= 0.0
                error("Volume Error :2908")
            end
        end
    end
end

function lagrangeElements(domain::Domain)

    delt = domain.deltatime
    domain.vnew = Vector{Float64}(undef, domain.numElem)
    domain.dxx = Vector{Float64}(undef, domain.numElem)
    domain.dyy = Vector{Float64}(undef, domain.numElem)
    domain.dzz = Vector{Float64}(undef, domain.numElem)

    domain.delx_xi = Vector{Float64}(undef, domain.numElem)
    domain.delx_eta = Vector{Float64}(undef, domain.numElem)
    domain.delx_zeta = Vector{Float64}(undef, domain.numElem)

    allElem = domain.numElem +  # local elem
            2*domain.sizeX*domain.sizeY + # plane ghosts
            2*domain.sizeX*domain.sizeZ + # row ghosts
            2*domain.sizeY*domain.sizeZ

    domain.delv_xi = Vector{Float64}(undef, allElem)
    domain.delv_eta = Vector{Float64}(undef, allElem)
    domain.delv_zeta = Vector{Float64}(undef, allElem)

    calcLagrangeElements(domain, delt)

    # Calculate Q.  (Monotonic q option requires communication)
    calcQForElems(domain)

    applyMaterialPropertiesForElems(domain)

    updateVolumesForElems(domain)

end

function calcCourantConstraintForElems(domain::Domain)

    dtcourant    = 1.0e+20
    courant_elem = -1

    qqc = domain.qqc
    length = domain.numElem

    qqc2 = 64.0 * qqc * qqc

    # Rewritten OpenMP code to sequential code
    courant_elem_per_thread = -1
    dtcourant_per_thread =  1.0e+20


    for i in 1:length
        indx = domain.matElemlist[i]

        dtf = domain.ss[indx] * domain.ss[indx]

        if domain.vdov[indx] < 0.0

        dtf = (dtf + qqc2 * domain.arealg[indx] * domain.arealg[indx]
                    * domain.vdov[indx]* domain.vdov[indx])
        end

        dtf = sqrt(dtf)

        dtf = domain.arealg[indx] / dtf

        #  determine minimum timestep with its corresponding elem
        if domain.vdov[indx] != 0.0
            if dtf < dtcourant_per_thread

                dtcourant_per_thread = dtf
                courant_elem_per_thread = indx
            end
        end
    end

    if dtcourant_per_thread < dtcourant
        dtcourant = dtcourant_per_thread
        courant_elem =  courant_elem_per_thread
    end


    # Don't try to register a time constraint if none of the elements
    # were active
    if courant_elem != -1
        domain.dtcourant = dtcourant
    end

    return nothing
end


function calcHydroConstraintForElems(domain::Domain)

    dthydro = 1.0e+20
    hydro_elem = -1
    dvovmax = domain.dvovmax
    length = domain.numElem

    # Rewritten OpenMP code to sequential code

    hydro_elem_per_thread = hydro_elem
    dthydro_per_thread = dthydro

    for i in 1:length
        indx = domain.matElemlist[i]

        if domain.vdov[indx] != 0.0
            dtdvov = dvovmax / (abs(domain.vdov[indx])+1.e-20)

            if dthydro_per_thread > dtdvov
                dthydro_per_thread = dtdvov
                hydro_elem_per_thread = indx
            end
        end
    end

    if dthydro_per_thread < dthydro
      dthydro = dthydro_per_thread
      hydro_elem =  hydro_elem_per_thread
    end

    if hydro_elem != -1
        domain.dthydro = dthydro
    end
    return nothing
end



function calcTimeConstraintsForElems(domain::Domain)
  # evaluate time constraint
  calcCourantConstraintForElems(domain::Domain)

  # check hydro constraint
  calcHydroConstraintForElems(domain::Domain)
end

function Isend(buf, dest::Integer, tag::Integer, comm::MPI.Comm)
    req = MPI.Request()
    ccall((:MPI_Isend, MPI.libmpi), Cint,
          (MPI.MPIPtr, Cint, MPI.MPI_Datatype, Cint, Cint, MPI.MPI_Comm, Ptr{MPI.MPI_Request}),
                  buf.data, buf.count, buf.datatype, dest, tag, comm, req)
    req.buffer = buf
    finalizer(MPI.free, req)
    return req
end

function lagrangeLeapFrog(domain, myRank)
        msgType = 0
	comm = MPI.COMM_WORLD
   
     data = MPI.Buffer(domain)
     if myRank == 0
      fromProc = 1
        MPI.Recv!(data, fromProc, msgType, comm)
     else
         otherRank = 0
         req = Isend(data, otherRank, msgType, comm)
      	 MPI.Wait!(req)
      end

   return nothing
end
