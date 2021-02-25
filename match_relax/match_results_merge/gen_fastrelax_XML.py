#用于对拼接(Merge)后的PDB文件(Match结果)自动化生成进行FastRelax的XML文件

import os,sys,shutil

#处理PDB中的行
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

merge_in_list = open(sys.argv[1],'r+')

NCAA_lis = ['23P','CYX']  #待查找的NCAA的三字符

for file in merge_in_list:
    #读入模板文件
    template = open('fastrelax_template.xml','r+')

    match_result_filename = file.split('_')
    merge_in = open(file.replace('\n',''), 'r+')

    #为每个Match结果的每个构象建立文件夹，方便后续优化及分析
    new_dir = f'Merge_UM_{match_result_filename[2]}_{match_result_filename[4]}'
    os.mkdir(new_dir)
    shutil.copy(file.replace('\n',''),new_dir)
    os.chdir(new_dir) #进入各个子文件夹

    residue_info = []
    #将PDB中的残基按顺序排列起来，用于find共价残基的位置

    #接下来的两个循环可以找出PDB文件中NCAA_lis中共价残基的序号(Rosetta序号，从1开始数)
    for line in merge_in:
        lis = generate_lis(line)
        res = (lis[3],int(lis[5]))
        if res not in residue_info:
            residue_info.append(res)

    for residue in residue_info:
        for NCAA in NCAA_lis:
            if NCAA in residue:
                position = residue_info.index(residue) + 1
                exec(f"""res_{NCAA}_pos = {position}""")
    
    print(res_CYX_pos,res_23P_pos)
    
    sh_merge_out = open(f'fastrelax_Merge_UM_{match_result_filename[2]}_{match_result_filename[4]}.xml','w+')

    #将fastrelax_template.xml格式化，因为不同的match结果匹配不同的位点，bond Mover的参数需要修改以适配
    for line in template:
        if '{}' in line:
            sh_merge_out.write(line.format(res_CYX_pos,res_23P_pos))
        else:
            sh_merge_out.write(line)

    os.chdir('../')