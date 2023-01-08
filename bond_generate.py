from bond_judge import bond_judge
from pymatgen.core import Molecule

def bond_generate(molecule : Molecule,mode=2):
    '''

    :param molecule:
    :param mode: 1--经典键长,2--原子成键半径
    :return: bond_dict
    '''
    bond_dict={}
    distance_matrix=molecule.distance_matrix
    sites=molecule.sites
    used_ele={x:1 for x in molecule.symbol_set}
    for site in sites:
        bond_dict[site.specie.name+str(used_ele[site.specie.name])]=[[],{}]
        used_ele[site.specie.name]+=1
    judger=bond_judge(2)
    shape=distance_matrix.shape
    elements=list(bond_dict.keys())
    for m in range(0,shape[0]):
        for n in range(m,shape[0]):
            site1_element=sites[m].specie.name;site1=elements[m]
            site2_element=sites[n].specie.name;site2=elements[n]
            distance=distance_matrix[(m,n)]
            try:
                judge=judger.judge(site1_element,site2_element,distance)
                if judge:
                    bond_dict[site1][0].append(site2)
                    bond_dict[site1][1][site2] = distance
                    bond_dict[site2][0].append(site1)
                    bond_dict[site2][1][site1] = distance
            except NotImplementedError:
                print(f"{site1} and {site2} is not in bond_lenth database")
            except:
                print("warning")
    return bond_dict