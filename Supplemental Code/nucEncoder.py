nucl = 'AAAAAAAA'

def nucEncoder(nucl):
    expo = len(nucl)
    valu = 0
    for n in nucl:
        tempsum = 4**(expo - 1)
        if n == 'A':
            valu += (tempsum * 1)
        elif n == 'T':
            valu += (tempsum * 2)
        elif n == 'C':
            valu += (tempsum * 3)
        elif n == 'G':
            valu += (tempsum * 4)
        expo = expo - 1
    valu = str(valu)
    return valu
