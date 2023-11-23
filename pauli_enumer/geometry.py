import ase.dft
class mole(object):
        
        def __init__(self, d):
                self.d=d


        def H_chain(self, hydrogen_num=2):
                'qubits = hydrogen number*2'
                geometry=[ ['H',[0, 0, i*self.d]] for i in range(hydrogen_num)]
                return geometry
   
        
        def H_net(self):
                'qubits = hydrogen number*2'
                geometry=[ ['H',[0, 0, 0]], ['H',[0, 0, self.d]],['H',[0, self.d, 0]],['H',[0, self.d, self.d]]]
                return geometry
   

        def LiH(self):
                print('bond length:',self.d)
                """qubit = 12"""
                LiH=[   ["Li",[0,0,0]], 
                        ["H",[0,0,self.d]]]
                return LiH


        def HF(self):
                """qubit = 12"""
                HF=[   ["F",[0,0,0]], 
                        ["H",[0,0,self.d]]]
                return HF


        def BeH2(self):
                """qubit = 14"""
                BeH2=[  ["H",[0,0,0]], 
                        ["Be",[0,0,self.d]], 
                        ["H",[0,0,2*self.d]]]
                return BeH2


        def H2O(self):
                """qubit = 14"""        
                H2O=  [ ["O", [+2.5369*self.d, -0.1550*self.d, 0.0000]],
                        ["H", [+3.0739*self.d, +0.1550*self.d, 0.0000]],
                        ["H", [+2.0000*self.d, +0.1550*self.d, 0.0000]]]
                return H2O


        def NH3(self):
                """qubit = 16"""
                NH3=  [ ["N",[4.877300262,4.771499932,0.004800000]],
                        ["H",[4.507300258,5.042399764,8.993999958]],
                        ["H",[4.507300258,5.511500239,0.744799972]],
                        ["H",[5.987299681,4.771499932,0.004800000]]]
                return NH3


        def CH4(self):
                """qubit =18 """
                CH4=  [ ["C", [1.000000000, 1.000000000, 1.000000000]],
                        ["H", [0.051054399, 0.051054399, 0.051054399]],
                        ["H", [1.948945601, 1.948945601, 0.051054399]],
                        ["H", [0.051054399, 1.948945601, 1.948945601]],
                        ["H", [1.948945601, 0.051054399, 1.948945601]]]
                return CH4


        def NaH(self):
                """qubit = 20"""
                NaH=  [ ["Na",[0,0,0]], 
                        ["H",[0,0,self.d]]]
                return NaH


        def N2(self):
                """qubit=20"""       
                N2=   [ ["N",[0,0,0]], 
                        ["N",[0,0,self.d]]] 
                return N2


        def C2H4(self):
                """qubit=28"""      
                C2H4= [ ["C", [+0.000000000, +0.000000000, -0.659227739]],
                        ["C", [+0.000000000, +0.000000000, +0.659227739]],
                        ["H", [+0.000000000, +0.929388479, -1.233149398]],
                        ["H", [-0.000000000, -0.929388479, -1.233149398]],
                        ["H", [-0.000000000, +0.929388479, +1.233149398]],
                        ["H", [+0.000000000, -0.929388479, +1.233149398]]]
                return C2H4

                
        def CO2(self):
                """qubit=28"""      
                geometry = [
                                ["O", [0., 0., 0. *self.d]],
                                ["C", [0., 0., 1. * self.d]],
                                ["O", [0., 0., 2. * self.d]]]

                return geometry

