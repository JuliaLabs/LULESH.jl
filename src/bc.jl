# Stuff needed for boundary conditions */
# BCs on each of 6 hexahedral faces (12 bits) */
const XI_M        = 0x00007
const XI_M_SYMM   = 0x00001
const XI_M_FREE   = 0x00002
const XI_M_COMM   = 0x00004

const XI_P        = 0x00038
const XI_P_SYMM   = 0x00008
const XI_P_FREE   = 0x00010
const XI_P_COMM   = 0x00020

const ETA_M       = 0x001c0
const ETA_M_SYMM  = 0x00040
const ETA_M_FREE  = 0x00080
const ETA_M_COMM  = 0x00100

const ETA_P       = 0x00e00
const ETA_P_SYMM  = 0x00200
const ETA_P_FREE  = 0x00400
const ETA_P_COMM  = 0x00800

const ZETA_M      = 0x07000
const ZETA_M_SYMM = 0x01000
const ZETA_M_FREE = 0x02000
const ZETA_M_COMM = 0x04000

const ZETA_P      = 0x38000
const ZETA_P_SYMM = 0x08000
const ZETA_P_FREE = 0x10000
const ZETA_P_COMM = 0x20000