'''
Reads MESA logs (both history and profile files) and stores data in dictionaries

MesaInfo is a class used to load a history file or a profile one. But not both at the same time.
It contains header information as well. In addition, it can clean MESA logs from backups
'''

import numpy as np
import gzip
import h5py
from pathlib import Path
import subprocess


__author__ = 'Adolfo Simaz Bunzel'
__credits__ = ['Adolfo Simaz Bunzel']
__license__ = 'GPL'
__version__ = '1.0'
__maintainer__ = 'Adolfo Simaz Bunzel'
__email__ = 'asimazbunzel@iar.unlp.edu.ar'



class AttributeMapper():
    '''simple mapper to access dictionary items as attributes'''

    def __init__(self, obj):
        self.__dict__['data'] = obj

    def __getattr__(self, attr):
        if attr in self.data:
            found_attr = self.data[attr]
            if isinstance(found_attr, dict):
                return AttributeMapper(found_attr)
            else:
                return found_attr
        else:
            raise AttributeError

    def __setattr__(self, attr, value):
        if attr in self.data:
            self.data[attr] = value
        else:
            raise NotImplementedError

    def __dir__(self):
        return self.data.keys()


class NoSingleValueFoundException(Exception):
    pass



class MesaInfo(object):
    '''
    Class that loads a MESA log file

    Parameters
    ----------
    filename (str): name of MESA file. It could be either a profile or a history output
    prune_data (bool): flag to cleanse output from backups (only for history files)

    Returns
    -------
    MesaInfo class. This class has the following attributes to access:
     - fname: same as the input file_name, i.e., it is a string
     - header: dictionary with header info. Methods for Python dictionaries are available.
     - data: dictionary with the column-type arrays
     - get() : method that given a name that corresponds to a valid column, it returns
               the corresponding array
    '''

    def __init__(self, filename, prune_data=True):

        if Path(filename).is_file():
            is_gzip = False
        else:
            if Path(filename + '.gz').is_file():
                is_gzip = True
            else:
                raise FileNotFoundError
       
        self.fname = filename
        if is_gzip: self.fname += '.gz'
        self.header = {}
        self.data_dict = {}
        is_history = 'history' in self.fname

        if is_gzip:
            file = gzip.open(self.fname, 'rb')
        else:
            file = open(self.fname, 'r')

        # First line is not used
        file.readline()

        # Header names
        if is_gzip:
            header_names = [name.decode('utf8') for name in file.readline().strip().split()]
        else:
            header_names = file.readline().strip().split()

        # After that are header names values
        if is_gzip:
            header_values = [val.decode('utf8') for val in file.readline().strip().split()]
        else:
            header_values = [val for val in file.readline().strip().split()]

        for i, name in enumerate(header_names): self.header[name] = header_values[i]

        # After header there is a blank line followed by an unused line.
        file.readline()
        file.readline()

        # Next are the column names
        if is_gzip:
            col_names = [name.decode('utf8') for name in file.readline().strip().split()]
        else:
            col_names = file.readline().strip().split()

        file.close()

        # Arrays are loaded using numpy an treated them as np.arrays
        file_data = np.loadtxt(self.fname, skiprows=6, unpack=True)

        for i, name in enumerate(col_names): self.data_dict[name] = np.array(file_data[i])

        if is_history:
            # Clean up log history
            # --------------------
            #
            #  The way it works is simple. It starts from the last model
            #  number and checks if that number is not repeated in the
            #  line before. This is combined with a mask so that it has
            #  a 0 (= True) for the lines that are not repeated and a 1
            #  (= False) for the repeated ones.
            #  After the mask is created, each array-type column is
            #  filtered with this mask using some numpy methods
            #  (numpy.ma.masked_array & compressed)
            models_number = self.data_dict['model_number']
            last_model = int(models_number[-1])
            mask = np.zeros(len(models_number))
            for i in range(len(models_number)-2, -1, -1):
                if int(models_number[i]) >= last_model:
                    mask[i] = 1
                else:
                    last_model = int(models_number[i])

            for name in col_names:
                self.data_dict[name] = np.ma.masked_array(self.data_dict[name], mask=mask).compressed()


    @property
    def data(self):
        return AttributeMapper(self.data_dict)

    
    def get(self, arg):
        '''
        Given a column name, it returns its values.

        Parameters
        -----------
        arg (str): column name

        Returns
        -------
        a: array of elements corresponding to the column name given
        '''
        return self.data[arg]


    def compress(self):
        '''
        Compress MESA output using gzip
        '''
        try:
            p = subprocess.Popen('gzip {}'.format(self.fname), stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, shell=True)
            stdout, stderr = p.communicate()
        except Exception:
            print('could not compress')



if __name__ == '__main__':
    pass
