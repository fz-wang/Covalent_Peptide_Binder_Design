#用于Match结果处理，可将Match的蛋白+配体+多肽的形式处理为蛋白+共价多肽的形式

import sys

#处理PDB中的行，加入处理第九个和第十个元素没有空格分离的特殊机制
#for example:
#ATOM      5  CB  VAL C  94       5.159 -57.155   5.053  1.00122.51           C  
#1.00和122.51之间没有空格，用split(' ')时会出现问题
def generate_lis(line):  
    lis = list(filter(None,line.split(' ')))
    if len(lis) < 9:
        return lis
    elif lis[9].count('.') == 2:
        lis.insert(10,lis[9][4:])
        lis[9] = lis[9][0:4]
        return lis
    else:
        return lis
        
#普通的行直接写入新文件
def process_regular_line(line):
    return line

#Match步骤中使用CYZ进行匹配，在后续操作中需要改为CYX
def process_CYZ_line(line):
    lis = generate_lis(line)
    lis[3] = 'CYX'
    return f'''{lis[0]:<6} {int(lis[1]):>4} {lis[2]:^4} {lis[3]} {lis[4]} {int(lis[5]):>3}     {lis[6]:>7} {lis[7]:>7} {lis[8]:>7} {lis[9]:>5}{lis[10]:>6}           {lis[11]:<3}{lis[12]}'''

#ALX在后续操作中需要改为UAA的名字，这里是23P
def process_ALX_line(line):
    lis = generate_lis(line)
    UAA_name = '23P'
    lis[3] = UAA_name
    return f'''{lis[0]:<6} {int(lis[1]):>4} {lis[2]:^4} {lis[3]} {lis[4]} {int(lis[5]):>3}     {lis[6]:>7} {lis[7]:>7} {lis[8]:>7} {lis[9]:>5}{lis[10]:>6}           {lis[11]:<3}{lis[12]}'''

#ALX在变为UAA后，原子数目变多了，故ALX后的所有原子的序号都需要增加相应的值
def process_after_ALXline(line):
    lig_atom_number = 10 #XLK配体有10个原子
    lis = generate_lis(line)
    return f'''{lis[0]:<6} {int(lis[1])+lig_atom_number:>4} {lis[2]:^4} {lis[3]} {lis[4]} {int(lis[5]):>3}     {lis[6]:>7} {lis[7]:>7} {lis[8]:>7} {lis[9]:>5}{lis[10]:>6}           {lis[11]:<3}{lis[12]}'''

#XLK行的处理，其中原子序号需要增加到ALX的对应值，原子类型按照字典进行修改，链编号对应ALX
def process_XLK_line(line):
    lis = generate_lis(line)
    lis[0], lis[3]= 'ATOM', '23P'
    return f'''{lis[0]:<6} {int(lis[1])+ALX_info['ALX_atom_maxnum']:>4} {XLK2NCAA_AtomType[lis[2]]:^4} {lis[3]} {ALX_info['ALX_chain']} {ALX_info['ALX_pos']:>3}     {lis[6]:>7} {lis[7]:>7} {lis[8]:>7} {lis[9]:>5}{lis[10]:>6}           {lis[11]:<3}{lis[12]}'''

#XLK和ALX在转化成NCAA时，原子坐标参照Match结果，但原子类型需要和Params文件中的保持一致
#使用字典
XLK2NCAA_AtomType ={'N1':'NG', 'C1':'CD', 'O1':'OE', 'C2':'CE', 'C3':'CZ', \
'H1':'1HG', 'H2':'1HE', 'H3':'2HE', 'H4':'1HZ', 'H5':'2HZ'}

match_result_file = open(sys.argv[1],'r+')

for file in match_result_file:
    match_in = open(file.replace('\n',""),'r+')
    match_out = open('Merge_'+file.replace('\n',""),'w')
    match_out_lis = []

    ALX_info = {}
    
    for line in match_in:   #读取一次，获取该Match结果中的ALX残基信息
        lis = generate_lis(line)
        if len(lis)>8 and lis[3] == "ALX":
            ALX_info['ALX_pos'] = int(lis[5])
            ALX_info['ALX_chain'] = lis[4]
            ALX_info['ALX_atom_maxnum'] = int(lis[1])
    match_in.seek(0)  #读取完一次后，光标位于文件末尾，需要将其移到文件开头，否则下一次读取将为空

    for line in match_in:   #再读取一次，写入新文件
        lis = generate_lis(line)

        if lis[0] not in {'ATOM','HETATM'} or len(lis) != 13:
            #match_out.write(process_regular_line(line))
            pass
        
        elif lis[3] == "CYZ":
            match_out_lis.append(process_CYZ_line(line))
       
        elif lis[3] == "ALX":
            match_out_lis.append(process_ALX_line(line))

        elif lis[4] == ALX_info['ALX_chain'] and int(lis[5])>ALX_info['ALX_pos']:
            match_out_lis.append(process_after_ALXline(line))
            

        elif lis[3] == "XLK":
            match_out_lis.append(process_XLK_line(line))
            
        else:
            match_out_lis.append(process_regular_line(line))
   
    match_out_lis.sort(key = lambda x: int(x[7:11])) #按照原子序号对文件排序，保证PDB中原子按序写入文件
    for x in match_out_lis:
        match_out.write(x)