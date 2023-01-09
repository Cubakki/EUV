import copy
import re
import os
import os.path as path
from bond_generate import bond_generate
from pymatgen.core import Molecule


class seperator:
    def __init__(self,file_path,output_path="default"):
        '''
        :param file_path: xyz文件路径
        :param output_path: 输出路径主目录
        '''
        self.load_path=file_path
        self.out_path=output_path
        self.load_xyz()
        self.output_init()
        self.bond_init()
        self.element_init()

    def load_xyz(self):
        '''
        读取xyz文件转换为pymatgen的Molecule类实例
        :return: True
        '''
        molecule = Molecule.from_file(self.load_path)
        self.molecule=molecule
        return True

    def output_init(self):
        if self.out_path == "default":
            self.out_path="./"+path.basename(self.load_path).split(".")[0]
        if not path.exists(self.out_path):
            os.mkdir(self.out_path)
        return True

    def bond_init(self):
        self.bond_dict=bond_generate(self.molecule)
        self.neighbour_list=[]
        self.site_list=[x for x in list(self.bond_dict.keys())]
        for site in self.site_list:
            self.neighbour_list.append([self.site_list.index(x) for x in self.bond_dict[site][0]])
        return True

    def element_init(self):
        '''
        生成元素--原子表
        :return:
        '''
        self.element_atom_dict={x:[] for x in self.molecule.symbol_set}
        for atom in self.site_list:
            self.element_atom_dict[re.findall("([^0-9]+)",atom)[0]].append(atom)
        return True

    def DFS(self,neighbour_list,start_index,visited):
        '''
        搜索neighbour_list(graph)
        :return: list
        '''
        for inferior_atom in neighbour_list[start_index]:
            if inferior_atom in visited:
                pass
            else:
                visited.append(inferior_atom)
                visited=self.DFS(neighbour_list,inferior_atom,visited)
        return visited

    def cut(self,core_atom,ligand_atom):
        """
        :param core_atom: core_atom,形式类似"Ti2"
        :param ligand_atom:ligand_atom,形式类似"C1"
        切断一个bond，返回切断后形成的core和ligand原子list
        :return:list，list(core_molecule_list,ligand_molecule_list) *其中元素是int，表示site在site_list中的位置索引
        """
        core_atom_index=self.site_list.index(core_atom);ligand_atom_index=self.site_list.index(ligand_atom)
        fixed_bond_dict=copy.deepcopy(self.bond_dict)
        fixed_bond_dict[core_atom][0].remove(ligand_atom)
        fixed_bond_dict[ligand_atom][0].remove(core_atom)
        fixed_neighbour_list=copy.deepcopy(self.neighbour_list)
        fixed_neighbour_list[core_atom_index].remove(ligand_atom_index)
        fixed_neighbour_list[ligand_atom_index].remove(core_atom_index)
        ligand_molecule_list=self.DFS(fixed_neighbour_list,ligand_atom_index,[ligand_atom_index])
        #以下拟可作为检验
        #core_molecule_list=self.DFS(fixed_neighbour_list,core_atom_index,[core_atom])
        core_molecule_list=list(set([x for x in range(0,len(self.site_list))])-set(ligand_molecule_list))
        return [core_molecule_list,ligand_molecule_list]

    def xyz_generator(self,x_molecule_list):
        """
        situation_write的配套函数，生成.xyz文件格式的文本
        :param x_molecule_list:
        :return: str
        """
        atom_num=0
        element_atom_dict={}
        text=""
        for index in x_molecule_list:
            atom_num+=1
            site=self.molecule.sites[index]
            coordinate=site.coords
            element=site.specie.name
            line="{} {} {} {}\n".format(element,coordinate[0],coordinate[1],coordinate[2])
            text+=line
            try:
                element_atom_dict[element].append(1)
            except:
                element_atom_dict[element]=[1]
        second_line=""
        for element_ in list(element_atom_dict.keys()):
            second_line+="{}{} ".format(element_,len(element_atom_dict[element_]))
        second_line=second_line.strip()
        head="{}\n{}\n".format(atom_num,second_line)
        return head+text

    def cut_information_generator(self,core,ligand):
        core_index=self.site_list.index(core)
        ligand_index=self.site_list.index(ligand)
        core_coor=self.molecule.sites[core_index].coords
        ligand_coor=self.molecule.sites[ligand_index].coords
        text="{} {} {} {}\n{} {} {} {}\n".format(core,core_coor[0],core_coor[1],core_coor[2],ligand,ligand_coor[0],ligand_coor[1],ligand_coor[2])
        return text

    def situation_write(self,core_molecule_list,ligand_mocule_list,situation_name,cut_information_text):
        '''
        输出一种core-ligand对至output_path下的situation_name文件夹
        :param core_molecule_list:
        :param ligand_mocule_list:
        :param cut_information_text:断键信息文本，由cut_information_generator方法生成
        :return: True
        '''
        save_path=self.out_path+"/"+situation_name
        core_text=self.xyz_generator(core_molecule_list)
        ligand_text=self.xyz_generator(ligand_mocule_list)
        if not path.exists(save_path):
            os.mkdir(save_path)
        f=open(save_path+"/core.xyz","w")
        f.write(core_text)
        f.close()
        f=open(save_path+"/ligand.xyz","w")
        f.write(ligand_text)
        f.close()
        f=open(save_path+"/bond_cutoff.txt","w")
        f.write(cut_information_text)
        f.close()
        return True


    def run(self,core,ligand):
        '''
        进行配体分离，每一种可能性都会输出core和ligand的xyz文件至output_path下的一个文件夹中
        :param core: eg:"Ti"
        :param ligand: eg:"C"
        :return: int(分离的组数)
        '''
        all_situation_list=[];bond_cut_infor=[]
        if not (core in self.molecule.symbol_set) and (ligand in self.molecule.symbol_set):
            raise ValueError("指定的core或配体不存在于化合物中")
        #主循环，遍历core元素的原子的成键情况，分离配体；切断core-ligand键后，以ligand原子为起点做graph的遍历，分离出core和ligand两组粒子位置
        for core_atom in self.element_atom_dict[core]:
            for associated_atom in self.bond_dict[core_atom][0]:
                if re.findall("[^0-9]+",associated_atom)[0] == ligand:
                    all_situation_list.append(self.cut(core_atom,associated_atom))
                    bond_cut_infor.append(self.cut_information_generator(core_atom,associated_atom))
        for i in range(0,len(all_situation_list)):
            self.situation_write(all_situation_list[i][0],all_situation_list[i][1],str(i+1),bond_cut_infor[i])
        print("成功生成{}组core-ligand对于{}目录".format(i+1,self.out_path))
