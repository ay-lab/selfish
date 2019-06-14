params.f1 = ''
params.f2 = ''
params.o = ''
params.r = 0

params.str = 'Hello World'

process splitLetters {

    output:
    file 'chunk_*' into letters mode flatten

    script:
    if (params.f1 != '' && params.f2 != '' && params.o != '' && params.r != 0)
        """
        printf '${params.str}' | split -b 6 - chunk_
        """
    else
        error "Invalid arguments. The arguments --f1, --f2, --o, --r are required."
}


process convertToUpper {

    input:
    file x from letters

    output:
    stdout result

    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

result.println { it.trim() }