from load_data import SRD101,CSD_radius

class bond_judge():
    def __init__(self,mode):
        self.bond_dict=SRD101()
        self.radius_dict=CSD_radius()
        self.mode=mode

    def mode(self,mode):
        '''
        更改模式
        :param mode:1--经典键长,2--成键半径
        :return: int (self.mode)
        '''
        self.mode=mode
        return self.mode

    def judge(self,site1,site2,distance):
        '''
        单位为埃
        :param site1: atom1
        :param site2: atom2
        :param bond_lenth: 几何距离
        :return: bool
        '''
        try:
            if self.mode==1:
                typical_lenth=float(self.bond_dict[site1][site2])
                if typical_lenth-0.6 <= distance <= typical_lenth+0.6:
                    return True
            elif self.mode==2:
                radius_sum=float(self.radius_dict[site1])+float(self.radius_dict[site2])
                if radius_sum-0.6 <= distance <= radius_sum+0.6:
                    return True
            else:
                raise Exception
        except:
            raise NotImplementedError(f"the bond between {site1} and {site2} is not in database")
        return False