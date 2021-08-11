   # assume communication to 6 neighbors by default
using MPI
function  get_neighbors(domain::Domain)
   rowMin = rowMax = colMin = colMax = planeMin = planeMax = true

   domain.m_rowLoc == 0 ? rowMin = false : rowMin = true
   domain.m_rowLoc == domain.m_tp - 1 ? rowMax = false : rowMax = true
   domain.m_colLoc == 0 ? colMin = false : colMin = true
   domain.m_colLoc == domain.m_tp - 1 ? colMax = false : colMax = true
   domain.m_planeLoc == 0 ? planeMin = false : planeMin = true
   domain.m_planeLoc == domain.m_tp - 1 ? planeMax = false : planeMax = true

   return rowMin, rowMax, colMin, colMax, planeMin, planeMax
end

# Assume 128 byte coherence
# Assume Float64 is an "integral power of 2" bytes wide
const CACHE_COHERENCE_PAD_REAL = div(128, sizeof(Float64))

function commRecv(domain::Domain, msgType, xferFields, dx, dy, dz, doRecv, planeOnly)
    if domain.comm === nothing
        return
    end

    error("Not implemented")
end

function commSend(domain::Domain, msgType, fields,
                  dx, dy, dz, doSend, planeOnly)

    if domain.comm === nothing
        return
    end

    xferFields = length(fields)
    error("not implemented")
end

function commSBN(domain::Domain, fields)
    if domain.comm === nothing
        return
    end
    error("not implemented")
end

function commMonoQ(domain::Domain)
    if domain.comm === nothing
        return
    end
    error("not implemented")
end

function commSyncPosVel(domain::Domain)
    if domain.comm === nothing
        return
    end
    error("not implemented")
end
