#!/usr/bin/env python3

from dataclasses import dataclass
import numpy as np

def read_bxsf(filename):

    class Flag(object):
        def __init__(self,flag_char):
            self.flag=False
            self.flag_char=flag_char
        def flag_switch(self,fline):
            if fline.find('BEGIN_'+self.flag_char)!=-1:
                self.flag=True
            elif fline.find('END_'+self.flag_char)!=-1:
                self.flag=False
            else:
                pass
    try:
        f=open(filename,'r')
    except:
        print("Error: Cannot open file")
    else:
        f.close()
        count=0
        axis=[0.]*3
        info=Flag('INFO')
        block_bandgrid_3d=Flag('BLOCK_BANDGRID_3D')
        bandgrid_3d_fermi=Flag('BANDGRID_3D')
        nband=-1
        band_num=[]

        for f in open(filename,'r'):
            info.flag_switch(f)
            block_bandgrid_3d.flag_switch(f)
            bandgrid_3d_fermi.flag_switch(f)
            if info.flag:
                tmp=f.split('#')[0]
                if tmp.find('Fermi Energy')!=-1:
                    EF=float(tmp.split(':')[1])
            if block_bandgrid_3d.flag:
                if (not bandgrid_3d_fermi.flag) and f.find('BANDGRID_3D')==-1:
                    index=f.strip()
            if bandgrid_3d_fermi.flag:
                if count==0:
                    pass
                elif count==1:
                    num_bands=int(f)
                    E_list=[[] for i in range(num_bands)]
                elif count==2:
                    tmp=f.split()
                    a,b,c=(int(t) for t in tmp)
                    sumk=a*b*c
                    k_list=np.array([[l%c,l%(c*b)//c ,l//(c*b)]
                                     for l in range(sumk)])
                elif count==3:
                    tmp=f.split()
                    center=[float(t) for t in tmp]
                elif 3<count<7:
                    tmp=f.split()
                    axis[count-4]=[float(t) for t in tmp]
                else:
                    if f.find('BAND:')!=-1:
                        tmp=f.split(':')
                        nband=nband+1
                        band_num.append(int(tmp[1]))
                    else:
                        tmp=[float(ff) for ff in f.split()]
                        E_list[nband]=E_list[nband]+tmp
                count=count+1
        axis=np.array(axis)
        E_list=np.array(E_list)

        return axis,E_list,band_num,index,EF,center,k_list

# Define a function to read the file
def read_frmsf(filename):
    with open(filename, 'r') as f:
        
        #kpoint spacing
        nk1, nk2, nk3 = map(int, f.readline().split())
        
        # grid type
        ishift = int(f.readline())

        # number of bands
        nbnd = int(f.readline())

        # Read the lattice vectors
        bvec1 = list(map(float, f.readline().split()))
        bvec2 = list(map(float, f.readline().split()))
        bvec3 = list(map(float, f.readline().split()))


        # Initialize arrays for energy and matrix elements
        eig = np.zeros((nbnd, nk1, nk2, nk3))
        x = np.zeros((nbnd, nk1, nk2, nk3))

        # Read energy values
        for ibnd in range(nbnd):
            for ik1 in range(nk1):
                for ik2 in range(nk2):
                    for ik3 in range(nk3):
                        try:
                            eig[ibnd, ik1, ik2, ik3] = float(f.readline())
                        except ValueError:
                            line_num = f.tell()
                            raise ValueError(f"Failed to convert line {line_num} to float")

        # Read matrix element values
        for ibnd in range(nbnd):
            for ik1 in range(nk1):
                for ik2 in range(nk2):
                    for ik3 in range(nk3):
                        try:
                            x[ibnd, ik1, ik2, ik3] = float(f.readline())
                        except ValueError:
                            continue
                            # line_num = f.tell()
                            # raise ValueError(f"Failed to convert line {line_num} to float")

    return bvec1, bvec2, bvec3, nk1, nk2, nk3, ishift, nbnd, eig, x

class bxsf:

    __slots__=['axis','E_list','band_num','index','EF','center','k_list']
    def __init__(self,axis=[0.]*3,elist=[],bnum=[],index='',ef=0.0,center=[0.,0.,0.],klist=[]):
        self.axis=axis
        self.E_list=elist
        self.band_num=bnum
        self.index=index
        self.EF=ef
        self.center=center
        self.k_list=klist


    @classmethod
    def from_file(cls, filename: str):
        try:
            axis, E_list, band_num, index, EF, center, k_list = read_bxsf(filename)
            return cls(axis, E_list, band_num, index, EF, center, k_list)
        except TypeError:
            print('cannot read %s'%filename)
            
            return None

    def write(self,filename):
        try:
            f=open(filename,'w')
        except:
            print('Cannot open file')
        else:
            f.write('  BEGIN_INFO\n'
                    +'       Fermi Energy: %23.15E\n'%self.EF
                    +'  END_INFO\n\n')
            f.write('  BEGIN_BLOCK_BANDGRID_3D\n'
                    +'       %s\n'%self.index
                    +'  BEGIN_BANDGRID_3D_fermi\n'
                    +'         %d\n'%len(self.E_list)
                    +'         %d       %d       %d\n'%tuple( i+1 for i in self.k_list[-1])
                    +'    %16.14f  %16.14f  %16.14f\n'%tuple(self.center))
            for a in self.axis:
                f.write('    %16.14f  %16.14f  %16.14f\n'%tuple(a))
            for i,El in enumerate(self.E_list):
                f.write('    BAND: %9s%d\n'%('',i+1))
                for e in El:
                    f.write('      %15.8E\n'%e)
            f.write('  END_BANDGRID_3D\n'
                    +'  END_BLOCK_BANDGRID_3D\n')
            f.close()

    def obtain_EF_band(self):
        cross_band=[]
        for e1 in self.E_list:
            if(max(e1)>self.EF and min(e1)<self.EF):
                cross_band.append(e1)
        cross_band=np.array(cross_band)
        
        return cross_band
    
    def get_2D_Fermi_data(self,kz):
        cross_band=self.obtain_EF_band()
        E_kz=[[]]*len(cross_band)
        k_kz=[[]]*len(cross_band)
        for j,ce in enumerate(cross_band):
            for i,e in enumerate(ce):
                if self.k_list[i][0]==kz:
                    E_kz[j]=E_kz[j]+[e-self.EF]
                    k_kz[j]=k_kz[j]+list([self.k_list[i]])

        E_kz=np.array(E_kz)
        k_kz=np.array(k_kz)

        return E_kz,k_kz
    

@dataclass  
class frmsf:
    bvec1: list
    bvec2: list
    bvec3: list
    nk1: int
    nk2: int
    nk3: int
    ishift: int
    nbnd: int
    eig: np.ndarray
    x: np.ndarray


    @classmethod
    def from_file(cls, filename: str):
        '''Creates frmsf object from file'''
        bvec1, bvec2, bvec3, nk1, nk2, nk3, ishift, nbnd, eig, x = read_frmsf(filename)

        return cls(bvec1, bvec2, bvec3, nk1, nk2, nk3, ishift, nbnd, eig, x)
    
    def to_bxsf(self, filename: str):
        '''Writes frmsf object to bxsf file'''
        # Create bxsf object
        bxsf_obj = bxsf()

        # Set axis
        bxsf_obj.axis = np.array([self.bvec1, self.bvec2, self.bvec3])

        # Set band_num
        bxsf_obj.band_num = np.arange(1, self.nbnd+1)

        # Set index
        bxsf_obj.index = 'Generated by WikipediaToothBrushes Inc'

        # Set k_list
        bxsf_obj.k_list = np.array([[ik1, ik2, ik3] for ik1 in range(self.nk1) for ik2 in range(self.nk2) for ik3 in range(self.nk3)])

        # Set center
        bxsf_obj.center = np.array([0, 0, 0])

        # Set E_list
        bxsf_obj.E_list = self.eig.reshape(self.nbnd, -1).T


        # Write to file
        bxsf_obj.write(filename)