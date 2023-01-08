import re

def SRD101():
    '''
    读取SRD101 database
    :return:a bond dict seems as {"H":{"H":"0.741","Li":"1.595",...},...}
    '''
    file=open("./bond_data/Standard Reference Database 101.txt")
    lines=file.readlines()
    for line in lines:
        line.replace("\n","")
    site1_list=lines[0].split("\t")
    bond_dict={x:{} for x in site1_list[1:]}
    #print(lines[0].split("\t"))
    for m in range(1,len(lines)):
        #site2_list = re.findall("\t(.+?)\t",lines[m])
        site2_list=lines[m].split("\t")
        for n in range(1,len(lines)):
            site1=site1_list[n]
            site2=site2_list[0]
            if not site2_list[n]==" ":
                bond_dict[site1][site2]=site2_list[n]
            else:
                pass
    #print(bond_dict)
    return bond_dict

def CSD_radius():
    '''
    读取经典成键半径
    :return: like {"H":"0.31"，"He"....}
    '''
    path="./bond_data/CSD_radius.txt"
    element_radius_dict={}
    with open(path,"r") as f:
        chart=f.readlines()
        f.close()
    for i in chart[1:]:
        block=i.split(" ")
        element_radius_dict[block[1]]=block[2]
    element_radius_dict["C"]="0.73"
    return element_radius_dict

if __name__=="__main__":
    pass