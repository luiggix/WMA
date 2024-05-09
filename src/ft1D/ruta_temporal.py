import os, sys
c_path = os.getcwd()
l_path = c_path.split(sep="/")
i_wma = l_path.index('WMA')
a_path = '/'.join(l_path[:i_wma])
src_path = '/WMA/src/ft1D'
if not(src_path in sys.path[0]):
    sys.path.insert(0, os.path.abspath(a_path + src_path)) 