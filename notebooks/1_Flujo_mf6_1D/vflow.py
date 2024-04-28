import numpy as np

class MeshDis():
    def __init__(self, nrow = 1, ncol = 1, nlay = 1, 
                 column_length = 1.0, row_length = 1.0,
                 top = 0.0, bottom = 0.0):
        self.__nx = ncol  # Number of columns
        self.__ny = nrow  # Number of rows
        self.__nz = nlay  # Number of layers
        self.__lx = row_length 
        self.__ly = column_length
        self.__lz = np.abs(top - bottom)            
        self.__dx = self.__calc_dx()
        self.__dy = self.__calc_dy()
        self.__dz = 1.0
        self.__top = top       # Top of the model 
        self.__bottom = bottom # Layer bottom elevation 

    def print(self):
        print('NX : {}'.format(self.__nx))
        print('NY : {}'.format(self.__ny))
        print('NZ : {}'.format(self.__nz))
        print('LX : {}'.format(self.__lx))
        print('LY : {}'.format(self.__ly))
        print('LZ : {}'.format(self.__lz))    
            
    @property
    def ncol(self):
        return self.__nx
        
    @ncol.setter
    def ncol(self, ncol):
        if ncol > 0:
            self.__nx = ncol
            self.__calc_dx()
        else:
            print('El valor {} no es válido para ncol, debe ser mayor que cero'.format(ncol))
            
    @property
    def nrow(self):
        return self.__ny
        
    @nrow.setter
    def nrow(self, nrow):
        if nrow > 0:
            self.__ny = nrow
            self.__calc_dy()
        else:
            print('El valor {} no es válido para nrow, debe ser mayor que cero'.format(nrow))

    @property
    def nlay(self):
        return self.__nz
        
    @nlay.setter
    def nlay(self, nlay):
        if nlay > 0:
            self.__nz = nlay
        ## TODO: calcular un delta para la direccion z
        else:
            print('El valor {} no es válido para nlay, debe ser mayor que cero'.format(nlay))

    @property
    def row_length(self):
        return self.__lx
        
    @row_length.setter
    def row_length(self, r_l):
        if r_l > 0:
            self.__lx = r_l
        else:
            print('El valor {} no es válido para row_length, debe ser mayor que cero'.format(r_l))

    @property
    def col_length(self):
        return self.__ly
        
    @col_length.setter
    def col_length(self, c_l):
        if c_l > 0:
            self.__ly = c_l
        else:
            print('El valor {} no es válido para col_length, debe ser mayor que cero'.format(c_l))

    @property
    def lay_length(self):
        return self.__lz
        
    @lay_length.setter
    def lay_length(self, l_l):
        if l_l > 0:
            self.__lz = l_l
        else:
            print('El valor {} no es válido para lay_length, debe ser mayor que cero'.format(l_l))

    @property
    def delr(self):
        return self.__dx

    @property
    def delc(self):
        return self.__dy

    @property
    def delz(self):
        return self.__dz

    @property
    def top(self):
        return self.__top
        
    @property
    def bottom(self):
        return self.__bottom
        
    def __calc_dx(self):
        return self.__lx / self.__nx

    def __calc_dy(self):
        return self.__ly / self.__ny 

    def __calc_dz(self):
        return self.__lz / self.__nz 
        
    def get_coords(self):
            return (np.linspace(0.5 * self.__dx, self.__lx - 0.5 *self.__dx, self.__nx),
                    np.linspace(0.5 * self.__dy, self.__ly - 0.5 *self.__dy, self.__ny),
                    np.linspace(0.5 * self.__dz, self.__lz - 0.5 *self.__dz, self.__nz))
        
    def get_grid(self):
        xg, yg, zg = np.meshgrid(np.linspace(0, self.__lx, self.__nx),
                                 np.linspace(0, self.__ly, self.__ny),
                                 np.linspace(0, self.__lz, self.__nz))

        
if __name__ == '__main__':
    mesh = MeshDis(
        nrow = 1,    # Number of rows
        ncol = 120,  # Number of columns
        nlay = 1,    # Number of layers
        row_length = 12.0,   # Length of rows
        column_length = 0.1, # Length of columns
        top = 1.0,   # Top of the model
        bottom = 0,  # Layer bottom elevation 
    )
    print(mesh.delr, mesh.delc)
    

