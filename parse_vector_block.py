# parse_vector_block.py

def parse_vector_block(lines,header,nstates):

    start=None

    for i,l in enumerate(lines):

        if header in l:
            start=i+1
            break

    if start is None:
        raise ValueError(header)

    data=[]

    for line in lines[start:]:

        toks=line.split()

        if len(toks)!=3:
            continue

        try:

            vec=list(map(float,toks))

            data.append(vec)

        except:

            continue

        if len(data)==nstates:
            break

    return np.array(data)