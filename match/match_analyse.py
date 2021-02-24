import math

def vector(a1,b1):
    return [a1[i]-b1[i] for i in range(3)]

def vector_len(a1,b1):
    vec = vector(a1,b1)
    return sum([x**2 for x in vec])**0.5

def angle(a1,b1,c1):
    b1a1 = vector(a1,b1)
    b1c1 = vector(c1,b1)
    DOT = sum([b1a1[i]*b1c1[i] for i in range(3)])
    cos = DOT / (vector_len(b1a1,[0,0,0])*vector_len(b1c1,[0,0,0]))
    if cos > 1.0:
        return 0.0
    elif cos < -1.0:
        return 180.0
    else:
        return round(math.acos(cos)*180.0/math.pi,1)

def normal(a1,a2,b1):
    a1a2 = vector(a2,a1)
    a1b1 = vector(b1,a1)
    return [(a1a2[1]*a1b1[2]-a1a2[2]*a1b1[1]), (a1a2[2]*a1b1[0]-a1a2[0]*a1b1[2]), (a1a2[0]*a1b1[1]-a1a2[1]*a1b1[0])]
    
def dihedral(a1,a2,b1,b2):
    norm1 = normal(a1,a2,b1)
    norm2 = normal(b2,b1,a2)
    direct = normal(norm1,(0,0,0),norm2)
    di_angle = round(180 - angle(norm1,[0,0,0],norm2),2)
    if angle(vector(a2,b1),(0,0,0),direct) > 0.0:
        return di_angle
    else:
        return -1*di_angle

match_result_file = open('result.txt','r+')

para_stat = []

for file in match_result_file:
    match_in = open(file.replace('\n',''),'r+')

    for line in match_in:
        lis_line = line.split(' ')
        lis = list(filter(None,lis_line))
        if lis[0] in {'ATOM','HETATM'}:
            xyz = [float(i) for i in lis[-7:-4]]
        if 'CA' in lis and 'CYZ' in lis:
            CA_CYZ = xyz
        if 'CB' in lis and 'CYZ' in lis:
            CB_CYZ = xyz
        if 'SG' in lis and 'CYZ' in lis:
            SG_CYZ = xyz
        if 'C3' in lis and 'XLK' in lis:
            C3_XLK = xyz
        if 'C2' in lis and 'XLK' in lis:
            C2_XLK = xyz
        if 'C1' in lis and 'XLK' in lis:
            C1_XLK = xyz
        if 'N1' in lis and 'XLK' in lis:
            N1_XLK = xyz
        if 'CB' in lis and 'ALX' in lis:
            CB_ALX = xyz
        if 'CA' in lis and 'ALX' in lis:
            CA_ALX = xyz
        if 'N' in lis and 'ALX' in lis:
            N_ALX  = xyz

    para = [[],[]]
    para[0].append(float(str(vector_len(SG_CYZ,C3_XLK))[0:4]))
    para[0].append(angle(C2_XLK,C3_XLK,SG_CYZ))
    para[0].append(angle(C3_XLK,SG_CYZ,CB_CYZ))
    para[0].append(dihedral(C1_XLK,C2_XLK,C3_XLK,SG_CYZ))
    para[0].append(dihedral(C2_XLK,C3_XLK,SG_CYZ,CB_CYZ))
    para[0].append(dihedral(C3_XLK,SG_CYZ,CB_CYZ,CA_CYZ))

    para[1].append(float(str(vector_len(N1_XLK,CB_ALX))[0:4]))
    para[1].append(angle(C1_XLK,N1_XLK,CB_ALX))
    para[1].append(angle(N1_XLK,CB_ALX,CA_ALX))
    para[1].append(dihedral(C2_XLK,C1_XLK,N1_XLK,CB_ALX))
    para[1].append(dihedral(C1_XLK,N1_XLK,CB_ALX,CA_ALX))
    para[1].append(dihedral(N1_XLK,CB_ALX,CA_ALX,N_ALX))

    print(para)
    #para_stat.append(para)

#print(para_stat)