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
