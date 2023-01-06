from load_data import SRD101

class bond_judge():
    def __init__(self):
        self.bond_dict=SRD101()

    def judge(self,site1,site2,distance):
        '''
        单位为埃
        :param site1: atom1
        :param site2: atom2
        :param bond_lenth: 几何距离
        :return: bool
        '''
        try:
            typical_lenth=self.bond_dict[site1][site2]
            if float(typical_lenth)-0.5 <= distance <= float(typical_lenth)+0.5:
                return True
        except:
            raise NotImplementedError(f"the bond between {site1} and {site2} is not in database")
        return False