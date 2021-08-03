function printUsage()
    @printf("Usage: \n")
    @printf("Unstructured grid:  %s -u <file.lmesh> \n", PROGRAM_FILE)
    @printf("Structured grid:    %s -s numEdgeElems \n", PROGRAM_FILE)
    @printf("\nExamples:\n")
    @printf("%s -s 45\n", PROGRAM_FILE)
    @printf("%s -u sedov15oct.lmesh\n", PROGRAM_FILE)
end

function getCacheCoherencePad(::T) where T
    return div(128, sizeof(T))
end

PAD_DIV(nbytes, align) = nbytes + div(align - 1, align)
PAD(nbytes, align) = PAD_DIV(nbytes,align) * align

function CalcElemVolume(
    x0, x1,x2, x3, x4, x5, x6, x7, y0, y1, y2, y3, y4, y5, y6, y7,
    z0, z1, z2, z3, z4, z5, z6, z7
)
    triple_product(x1, y1, z1, x2, y2, z2, x3, y3, z3) = ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))
    twelveth = 1.0/12.0

    dx61 = x6 - x1
    dy61 = y6 - y1
    dz61 = z6 - z1

    dx70 = x7 - x0
    dy70 = y7 - y0
    dz70 = z7 - z0

    dx63 = x6 - x3
    dy63 = y6 - y3
    dz63 = z6 - z3

    dx20 = x2 - x0
    dy20 = y2 - y0
    dz20 = z2 - z0

    dx50 = x5 - x0
    dy50 = y5 - y0
    dz50 = z5 - z0

    dx64 = x6 - x4
    dy64 = y6 - y4
    dz64 = z6 - z4

    dx31 = x3 - x1
    dy31 = y3 - y1
    dz31 = z3 - z1

    dx72 = x7 - x2
    dy72 = y7 - y2
    dz72 = z7 - z2

    dx43 = x4 - x3
    dy43 = y4 - y3
    dz43 = z4 - z3

    dx57 = x5 - x7
    dy57 = y5 - y7
    dz57 = z5 - z7

    dx14 = x1 - x4
    dy14 = y1 - y4
    dz14 = z1 - z4

    dx25 = x2 - x5
    dy25 = y2 - y5
    dz25 = z2 - z5

    # 11 + 3*14
    volume = (
        triple_product(
            dx31 + dx72, dx63, dx20,
            dy31 + dy72, dy63, dy20,
            dz31 + dz72, dz63, dz20
        ) +
        triple_product(
            dx43 + dx57, dx64, dx70,
            dy43 + dy57, dy64, dy70,
            dz43 + dz57, dz64, dz70
        ) +
        triple_product(
            dx14 + dx25, dx61, dx50,
            dy14 + dy25, dy61, dy50,
            dz14 + dz25, dz61, dz50
        )
    )
    volume *= twelveth;
    return volume
end

CalcElemVolume(x, y, z) = CalcElemVolume(
    x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8],
    y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8],
    z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]
)
