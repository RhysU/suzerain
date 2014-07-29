#!/usr/bin/env python
"""Usage: reacting_summary.py HDF5FILES
Obtain summary file from apps/reacting HDF5FILES,
or load an existing summary file.
Options:
  -h  --help      This help message.
"""
import ar
import gb
import getopt
import h5py
from itertools import izip
import numpy as np
import shutil
import sys


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class summary():

    def __init__(self, **kwargs):

        self.hdf5files = kwargs['hdf5files']
        self.nfiles    = len(self.hdf5files)

        self.load_metadata(self.hdf5files[0])

        # Get sequences to construct variable names
        self.get_vars_sequence()

        # Generate list of summary variables
        self.scalars = ['t']

        self.planes  = []
        self.planes.extend([('bar_rho'     ,              'bar_rho',            None)])
        self.planes.extend([('bar_rho_E'   ,              'bar_rho_E',          None)])
        self.planes.extend([('bar_rho_'+x,                'bar_rho_u',          i   )   for i,  x        in enumerate(self.seq_ui)])

        if self.Ns > 1:
            self.planes.extend([('bar_rho'+x                 ,'bar_rho_s'            , i     ) for i,  x     in enumerate(self.seq_ci)])

        self.planes.extend([('bar_Cv' ,                   'bar_Cv' ,              None)])
        self.planes.extend([('bar_Cp' ,                   'bar_Cp' ,              None)])
        self.planes.extend([('bar_gamma' ,                'bar_gamma' ,           None)])
        self.planes.extend([('bar_D0' ,                   'bar_D0' ,              None)])
        self.planes.extend([('bar_E'  ,                   'bar_E'  ,              None)])
        self.planes.extend([('bar_M'  ,                   'bar_M'  ,              None)])
        self.planes.extend([('bar_M_M',                   'bar_M_M',              None)])
        self.planes.extend([('bar_T'  ,                   'bar_T'  ,              None)])
        self.planes.extend([('bar_T_T',                   'bar_T_T',              None)])
        self.planes.extend([('bar_a'  ,                   'bar_a'  ,              None)])
        self.planes.extend([('bar_a_a',                   'bar_a_a',              None)])
        self.planes.extend([('bar_grad'+x+'_T',           'bar_grad_T',           i     ) for i,  x        in enumerate(self.seq_xi)])
        self.planes.extend([('bar_kappa',                 'bar_kappa',            None)])
        self.planes.extend([('bar_kappa_grad'+x+'_T',     'bar_kappa_grad_T',     i     ) for i,  x        in enumerate(self.seq_xi)])
        self.planes.extend([('bar_mu',                    'bar_mu',               None)])
        self.planes.extend([('bar_mu_S'+x+y,              'bar_mu_S',             i     ) for i, (x, y)    in enumerate(self.seq_xi_xj)])
        self.planes.extend([('bar_mu_div_u',              'bar_mu_div_u',         None)])
        self.planes.extend([('bar_mu_mu'   ,              'bar_mu_mu'   ,         None)])
        self.planes.extend([('bar_nu'      ,              'bar_nu'      ,         None)])
        self.planes.extend([('bar_p'       ,              'bar_p'       ,         None)])
        self.planes.extend([('bar_p_div_u' ,              'bar_p_div_u' ,         None)])
        self.planes.extend([('bar_p_p'     ,              'bar_p_p'     ,         None)])
        self.planes.extend([('bar_p_'+x    ,              'bar_p_u'     ,         i     ) for i,  x        in enumerate(self.seq_ui)])
        self.planes.extend([('bar_qb'      ,              'bar_qb',               None)])
        self.planes.extend([('bar_rho_E_'+x,              'bar_rho_E_u',          i     ) for i,  x        in enumerate(self.seq_ui)])
        self.planes.extend([('bar_rho_T'   ,              'bar_rho_T',            None)])
        self.planes.extend([('bar_rho_T_'+x,              'bar_rho_T_u',          i     ) for i,  x        in enumerate(self.seq_ui)])
        self.planes.extend([('bar_rho_grad'+x+'_T',       'bar_rho_grad_T',       i     ) for i,  x        in enumerate(self.seq_xi)])
        self.planes.extend([('bar_rho_mu'  ,              'bar_rho_mu'    ,       None)])
        self.planes.extend([('bar_rho_rho' ,              'bar_rho_rho'   ,       None)])
        self.planes.extend([('bar_'+x+'_'+y,              'bar_u_u'       ,       i     ) for i, (x, y)    in enumerate(self.seq_ui_uj)])
        self.planes.extend([('bar_rho_'+x+'_'+y,          'bar_rho_u_u'   ,       i     ) for i, (x, y)    in enumerate(self.seq_ui_uj)])
        self.planes.extend([('bar_rho_'+x+'_'+y+'_'+z,    'bar_rho_u_u_u' ,       i     ) for i, (x, y, z) in enumerate(self.seq_ui_uj_uk)])
        self.planes.extend([('bar_sym'+x+y+'_grad_u',     'bar_sym_grad_u',       i     ) for i, (x, y)    in enumerate(self.seq_xi_xj)])
        self.planes.extend([('bar_sym'+x+y+'_rho_grad_u', 'bar_sym_rho_grad_u',   i     ) for i, (x, y)    in enumerate(self.seq_xi_xj)])
        self.planes.extend([('bar_tau'+x+y ,              'bar_tau'       ,       i     ) for i, (x, y)    in enumerate(self.seq_xi_xj)])
        self.planes.extend([('bar_tau_colon_grad_u' ,     'bar_tau_colon_grad_u', None)])
        self.planes.extend([('bar_tauu'+x,                'bar_tau_u'     ,       i     ) for i,  x        in enumerate(self.seq_xi)])
        self.planes.extend([('bar_'+x,                    'bar_u'         ,       i     ) for i,  x        in enumerate(self.seq_ui)])

        if self.Ns > 1:
            self.planes.extend([('bar_om'+x                  ,'bar_om_s'             , i     ) for i,  x     in enumerate(self.seq_ci)])
            self.planes.extend([('bar_rho_D0'                ,'bar_rho_D0'           , None)])
            self.planes.extend([('bar_rho_Ds_grad'+x+'_cs'   ,'bar_rho_Ds_grad_cs'   , i     ) for i,  x     in enumerate(self.seq_xi)])
            self.planes.extend([('bar_rho_Ds_grad'+x+'_cs_hs','bar_rho_Ds_grad_cs_hs', i     ) for i,  x     in enumerate(self.seq_xi)])
            self.planes.extend([('bar_rho'+x+'_'+y           ,'bar_rho_s_u'          , i     ) for i, (x, y) in enumerate(self.seq_ci_uj)])


        if self.processfiles:
            # Check if each input variable is in files,
            # remove from list otherwise
            for index, file in enumerate(self.hdf5files):
                print 'Checking file', self.hdf5files[index]
                f = h5py.File(file, "r")

                # check for input data
                for var, in_var, varindex in list(self.planes):
                    if in_var not in f.keys():
                        self.planes.remove((var, in_var, varindex))
                        print "--- removing", var, " --- ", in_var, ", not in file"
                f.close()
        else:
            # Check if each processed variable is in files,
            # remove from list otherwise
            for index, file in enumerate(self.hdf5files):
                print 'Checking file', self.hdf5files[index]
                f = h5py.File(file, "r")

                # check for processed  data
                for var, in_var, varindex in list(self.planes):
                    if var not in f.keys():
                        self.planes.remove((var, in_var, varindex))
                        print "--- removing", var, ", not in file"
                f.close()

        # Set quantities to ignore if not loaded
        self.set_ignore()

        # Input variables
        self.inputscalars = ['t']
        self.inputdata    = list(set([planes[1] for planes in self.planes]))

        # init scalar storage
        for scalar in self.scalars:
            self.init_line(scalar, self.nfiles)

        # init plane storage
        for plane in [planes[0] for planes in self.planes]:
            #print "... init storage for", plane
            self.init_plane(plane, self.nfiles, self.Ny)

        # Initialize arsel storage for all variables
        for var in [planes[0] for planes in self.planes]:
            self.init_line('mean_'   +var, self.Ny)
            self.init_line('sigma_'  +var, self.Ny)
            self.init_line('eff_N_'  +var, self.Ny)
            self.init_line('eff_var_'+var, self.Ny)
            self.init_line('T0_'     +var, self.Ny)

        # Populate from samples and compute arsel stuff
        if self.processfiles:
            self.process_from_samples()
        else:
            self.load(self.hdf5files[0])

        # Get all mean profile derived quantities
        self.get_mean_derived()

        # Get all boundary layer derived quantities
        self.get_bl()

        # Get transformed quantities (eg, van Driest)
        self.get_transformed()

        # Get quantities scaled by wall units
        self.get_scaled_plus()

        # Get quantities scaled by semi-local units
        self.get_scaled_sl()


    def process_from_samples(self):
        # populate data from files
        for index, file in enumerate(self.hdf5files):
            print 'Processing data from file', self.hdf5files[index]
            f = h5py.File(file, "r")

            # load input scalars
            for scalar in self.inputscalars:
                self.__dict__['in_'+scalar] = np.squeeze(f[scalar][()])

            # load input data
            for data in list(self.inputdata):
                self.__dict__['in_'+data]   = np.squeeze(f[data][()])
            f.close()

            # populate scalars
            for scalar in self.scalars:
                #self.__dict__[scalar][index] = howtocompute[scalar]()
                self.__dict__[scalar][index] = self.get_t()

            # populate arrays
            for var, data, varindex in list(self.planes):
                self.__dict__[var][index,:] = self.get_scalar(data, varindex)

        # Compute mean with arsel
        # Compute arsel for all variables
        if self.nfiles == 1:
            for var in [planes[0] for planes in self.planes]:
                #print "... Computing for", var
                for index in xrange(self.Ny):
                    self.__dict__['mean_'   +var][index] = self.__dict__[var][0,index]
                    self.__dict__['sigma_'  +var][index] = 0
                    self.__dict__['eff_N_'  +var][index] = 0
                    self.__dict__['eff_var_'+var][index] = 0
                    self.__dict__['T0_'     +var][index] = 0
        else:
            for var in [planes[0] for planes in self.planes]:
                #print "... Computing for", var
                for index in xrange(self.Ny):
                    t = ar.arsel(self.__dict__[var][:,index])
                    self.__dict__['mean_'   +var][index] = t.mu[0]
                    self.__dict__['sigma_'  +var][index] = t.mu_sigma[0]
                    self.__dict__['eff_N_'  +var][index] = t.eff_N[0]
                    self.__dict__['eff_var_'+var][index] = t.eff_var[0]
                    self.__dict__['T0_'     +var][index] = t.T0[0]


    def init_line(self, var, a=1):
        self.__dict__[var] = np.empty([a])


    def init_plane(self, var, a=1, b=1):
        self.__dict__[var] = np.empty([a,b])


    def clone(self, outfile="summary.h5"):
        print "Cloning metadata"
        f = h5py.File(self.hdf5files[0], "r")
        g = h5py.File(outfile, "w")

        # Entries to exclude
        exclude_vars  = ['t']

        # Roots of entries to exclude from cloning
        exclude_roots = ['bar', 'rho', 'bulk', 'min', 'max',
                 'xmin', 'xmax', 'zmin', 'zmax', 'fneg',
                 'twopoint']

        # List of names to copy
        names = [x for x in f.keys()
                 if not(x in exclude_vars or any([x.startswith(r) for r in exclude_roots]))]

        # Entries to copy with a different name
        maps  = [('y', 'collocation_points_y')]

        # Straight copy with same name
        for key in names:
            g['/'].copy(f[key], key)

        # Copies entries with a new name in the cloned file
        for destination, origin in maps:
            g['/'].copy(f[origin], destination)

        g.close()
        f.close()


    def save(self, outfile="summary.h5"):

        # clone file with metadata
        self.clone()

        print "Saving data"
        f = h5py.File(outfile, "a")

        # Scalar data
        for scalar in self.scalars:
            f[scalar] = self.__dict__[scalar]

        # Data from field averages
        for plane in [planes[0] for planes in self.planes]:
            f[plane] = self.__dict__[plane]
            f[plane].attrs.create('mean'   , self.__dict__['mean_'   +plane])
            f[plane].attrs.create('sigma'  , self.__dict__['sigma_'  +plane])
            f[plane].attrs.create('eff_N'  , self.__dict__['eff_N_'  +plane])
            f[plane].attrs.create('eff_var', self.__dict__['eff_var_'+plane])
            f[plane].attrs.create('T0'     , self.__dict__['T0_'     +plane])

        f.close()


    def load_metadata(self, file):

        print "Loading metadata from", file
        f = h5py.File(file, "r")

        # Input metadata
        self.metascalars = ['Nx', 'Ny', 'Nz', 'k', 'Lx', 'Ly', 'Lz']
        self.metalines   = ['breakpoints_y',
                            'collocation_points_x',
                            'collocation_points_y',
                            'collocation_points_z',
                            'Dy0T', 'Dy1T', 'Dy2T',
                            'integration_weights']

        # Open first file to load metadata
        f = h5py.File(file, "r")

        # Load reference metadata
        for scalar in self.metascalars:
            self.__dict__[scalar] = np.squeeze(f[scalar][()])

        for line in self.metalines:
            self.__dict__[line] = np.squeeze(f[line][()])

        # Grab number of species
        self.Ns=f['antioch_constitutive_data'].attrs['Ns'][0]

        # Grab species names
        self.sname= []
        for s in xrange(self.Ns):
            self.sname.append(f['antioch_constitutive_data'].attrs['Species_'+str(s)])

        # Grab delta growth rate
        self.grD = f['largo'].attrs['grdelta'][0]

        # Grab largo formulation
        self.formulationi = f['largo'].attrs['formulation']

        # Grab mean amplitude growth rates
        if 'largo_gramp_mean' in f:
            self.largo_gramp = f['largo_gramp_mean'].value
        else:
            self.largo_gramp = None #np.zeros((6,1))

        # Grab baseflow coefficients
        if 'largo_baseflow' in f:
            self.coeff_baseflow = np.squeeze(f['largo_baseflow'][()])
        else:
            self.coeff_baseflow = None #np.zeros((6,2))

        if 'largo_baseflow_dx' in f:
            self.coeff_baseflow_dx = np.squeeze(f['largo_baseflow_dx'][()])
        else:
            self.coeff_baseflow_dx = None #np.zeros((6,2))

        # Process reference metadata
        # Check if the first file is a sample file
        # or a summary file (by looking for y)
        if 'y' not in f.keys():
            self.y   = np.squeeze(f[line][()])
            self.processfiles = True
        else:
            if self.nfiles > 1:
                #f.close()
                print >>sys.stderr, "The first file is a summary file, only one can be loaded"
            else:
                self.y   = self.collocation_points_y
                self.processfiles = False

        f.close()

        self.yb  = self.breakpoints_y
        self.iw  = self.integration_weights
        self.D0T = gb.gb2ge(self.Dy0T, self.Ny, self.k-2)
        self.D1T = gb.gb2ge(self.Dy1T, self.Ny, self.k-2)
        self.D2T = gb.gb2ge(self.Dy2T, self.Ny, self.k-2)

        # Matrices to compute first and second derivatives
        # from data on collocation points.
        # Eg, first derivative computed as df/dy = D1 * D0^-1 f(y)
        self.invD0T = np.linalg.solve(self.D0T, np.eye(self.Ny))
        self.invD0T_D1T = np.dot(self.invD0T, self.D1T)
        self.invD0T_D2T = np.dot(self.invD0T, self.D2T)
        return

    def load(self, file):

        print "Loading data from summary file", file

        # Open first file to load data
        f = h5py.File(file, "r")

        # Select variables to load
        for var, in_var, index in list(self.planes):
            if var not in f.keys():
                self.planes.remove((var, in_var, index))

        # Scalar data
        for scalar in self.scalars:
            self.__dict__[scalar] = np.squeeze(f[scalar][()])

        # Data from field averages
        for plane in [planes[0] for planes in self.planes]:
            self.__dict__[plane]            = np.squeeze(f[plane][()])
            self.__dict__['mean_'   +plane] = np.squeeze(f[plane].attrs.get('mean'   ))
            self.__dict__['sigma_'  +plane] = np.squeeze(f[plane].attrs.get('sigma'  ))
            self.__dict__['eff_N_'  +plane] = np.squeeze(f[plane].attrs.get('eff_N'  ))
            self.__dict__['eff_var_'+plane] = np.squeeze(f[plane].attrs.get('eff_var'))
            self.__dict__['T0_'     +plane] = np.squeeze(f[plane].attrs.get('T0'     ))

        f.close()
        return


    def get_t(self):
        return self.in_t


    def get_scalar(self, data, varindex):
        if varindex is None:
            return np.dot(self.__dict__['in_'+data], self.D0T)
        else:
            return np.dot(self.__dict__['in_'+data][varindex,:], self.D0T)


    def get_vars_sequence(self):
        # Coordinates
        self.seq_xi       = ['x', 'y', 'z']
        self.seq_xi_xi    = [(x,x) for    x in self.seq_xi]
        self.seq_xi_xj    = [(x,y) for i, x in enumerate(self.seq_xi)
                                   for    y in self.seq_xi[i:]]
        # Velocities
        self.seq_ui       = ['u', 'v', 'w']
        self.seq_ui_ui    = [(x,x) for    x in self.seq_ui]
        self.seq_ui_uj    = [(x,y) for i, x in enumerate(self.seq_ui)
                                   for    y in self.seq_ui[i:]]
        self.seq_ui_uj_uk = [(x,y,z) for i, x in enumerate(self.seq_ui)
                                     for j, y in enumerate(self.seq_ui[i:])
                                     for    z in           self.seq_ui[i+j:]]
        # Species
        self.seq_ci       = self.sname
        self.seq_ci_ci    = [(x,x) for    x in self.seq_ci]

        # Species-velocities
        self.seq_ci_uj    = [(x,y) for    x in self.seq_ui
                                   for    y in self.seq_ci]

        return

    def get_points_below(self, target, y):
        pb = self.Ny
        for j, pb in enumerate(y):
            if (y[j] > target):
                pb = j-1
                break
        return pb

    def set_ignore(self):
        # Set quantities to ignore 
        self.mean_bar_rho_T   = None
        self.mean_bar_p_u     = None
        self.mean_bar_p_v     = None
        self.mean_bar_p_w     = None
        self.mean_bar_u_u     = None
        self.mean_bar_u_v     = None
        self.mean_bar_u_w     = None
        self.mean_bar_v_v     = None
        self.mean_bar_v_w     = None
        self.mean_bar_w_w     = None
        self.mean_bar_T_T     = None 
        self.mean_bar_rho_rho = None
        self.mean_bar_p_p     = None
        self.mean_bar_a_a     = None
        self.mean_bar_M_M     = None
        self.mean_bar_mu_mu   = None
        self.mean_bar_M       = None
        self.mean_bar_gamma   = None
        self.mean_bar_Cp      = None
        self.mean_bar_Cv      = None

    def get_mean_derived(self):
        # Favre averages
        self.fav_u = self.get_fav(self.mean_bar_rho_u, self.mean_bar_rho)
        self.fav_v = self.get_fav(self.mean_bar_rho_v, self.mean_bar_rho)
        self.fav_w = self.get_fav(self.mean_bar_rho_w, self.mean_bar_rho)

        self.fav_T = self.get_fav(self.mean_bar_rho_T, self.mean_bar_rho)
        self.fav_E = self.get_fav(self.mean_bar_rho_E, self.mean_bar_rho)
        self.fav_H = self.fav_E + self.mean_bar_p / self.mean_bar_rho

        # Favre fluctuations
        self.upp   = self.mean_bar_u - self.fav_u
        self.vpp   = self.mean_bar_v - self.fav_v
        self.wpp   = self.mean_bar_w - self.fav_w
        if self.fav_T is not None:
            self.Tpp   = self.mean_bar_T - self.fav_T

        # Reyolds stresses
        self.rho_upp_upp = self.get_rho_fpp_gpp(self.mean_bar_rho_u_u, self.mean_bar_rho, self.fav_u, self.fav_u)
        self.rho_upp_vpp = self.get_rho_fpp_gpp(self.mean_bar_rho_u_v, self.mean_bar_rho, self.fav_u, self.fav_v)
        self.rho_upp_wpp = self.get_rho_fpp_gpp(self.mean_bar_rho_u_w, self.mean_bar_rho, self.fav_u, self.fav_w)
        self.rho_vpp_vpp = self.get_rho_fpp_gpp(self.mean_bar_rho_v_v, self.mean_bar_rho, self.fav_v, self.fav_v)
        self.rho_vpp_wpp = self.get_rho_fpp_gpp(self.mean_bar_rho_v_w, self.mean_bar_rho, self.fav_v, self.fav_w)
        self.rho_wpp_wpp = self.get_rho_fpp_gpp(self.mean_bar_rho_w_w, self.mean_bar_rho, self.fav_w, self.fav_w)

        # bar_rho_uipp_ujpp_ukpp
        for ui, uj, uk in self.seq_ui_uj_uk:
            self.__dict__                           ['rho_%spp_%spp_%spp'      % (ui, uj, uk)] =\
              self.get_rho_fpp_gpp_hpp(self.__dict__['mean_bar_rho_%s_%s_%s'   % (ui, uj, uk)]
                                      ,self.__dict__['mean_bar_rho_%s_%s'      % (ui, uj    )]
                                      ,self.__dict__['mean_bar_rho_%s_%s'      % (    uj, uk)]
                                      ,self.__dict__['mean_bar_rho_%s_%s'      % (ui,     uk)]
                                      ,self.mean_bar_rho
                                      ,self.__dict__['fav_%s'                  %  ui         ]
                                      ,self.__dict__['fav_%s'                  %      uj     ]
                                      ,self.__dict__['fav_%s'                  %          uk ])

        # Favre Reyolds stresses
        self.fav_upp_upp = self.get_fav(self.rho_upp_upp, self.mean_bar_rho)
        self.fav_upp_vpp = self.get_fav(self.rho_upp_vpp, self.mean_bar_rho)
        self.fav_upp_wpp = self.get_fav(self.rho_upp_wpp, self.mean_bar_rho)
        self.fav_vpp_vpp = self.get_fav(self.rho_vpp_vpp, self.mean_bar_rho)
        self.fav_vpp_wpp = self.get_fav(self.rho_vpp_wpp, self.mean_bar_rho)
        self.fav_wpp_wpp = self.get_fav(self.rho_wpp_wpp, self.mean_bar_rho)

        # Temperature-velocity Favre fluctuations
        self.rho_Tpp_upp = self.get_rho_fpp_gpp(self.mean_bar_rho_T_u, self.mean_bar_rho, self.fav_T, self.fav_u)
        self.rho_Tpp_vpp = self.get_rho_fpp_gpp(self.mean_bar_rho_T_v, self.mean_bar_rho, self.fav_T, self.fav_v)
        self.rho_Tpp_wpp = self.get_rho_fpp_gpp(self.mean_bar_rho_T_w, self.mean_bar_rho, self.fav_T, self.fav_w)

        # Favre Temperature velocity fluctuations
        self.fav_Tpp_upp = self.get_fav(self.rho_Tpp_upp, self.mean_bar_rho)
        self.fav_Tpp_vpp = self.get_fav(self.rho_Tpp_vpp, self.mean_bar_rho)
        self.fav_Tpp_wpp = self.get_fav(self.rho_Tpp_wpp, self.mean_bar_rho)

        # Pressure - Favre velocity fluctuations
        self.pp_upp = self.get_fp_gpp(self.mean_bar_p_u, self.mean_bar_p, self.mean_bar_u)
        self.pp_vpp = self.get_fp_gpp(self.mean_bar_p_v, self.mean_bar_p, self.mean_bar_v)
        self.pp_wpp = self.get_fp_gpp(self.mean_bar_p_w, self.mean_bar_p, self.mean_bar_w)

        # Tau fluctutions dot Favre velocity fluctuations
        # NOTE: Compute these ones explicitly
        self.taupuppx = self.mean_bar_tauux - ( self.mean_bar_tauxx * self.mean_bar_u
                                              + self.mean_bar_tauxy * self.mean_bar_v
                                              + self.mean_bar_tauxz * self.mean_bar_w)
        self.taupuppy = self.mean_bar_tauuy - ( self.mean_bar_tauxy * self.mean_bar_u
                                              + self.mean_bar_tauyy * self.mean_bar_v
                                              + self.mean_bar_tauyz * self.mean_bar_w)
        self.taupuppz = self.mean_bar_tauuz - ( self.mean_bar_tauxz * self.mean_bar_u
                                              + self.mean_bar_tauyz * self.mean_bar_v
                                              + self.mean_bar_tauzz * self.mean_bar_w)

        # Tau fluctuations dot derivative of Favre velocity fluctuations
        # NOTE: Compute these ones explicitly
        #       Assume duppdx = duppdz = 0
        self.taup_colon_grad_upp = (self.mean_bar_tau_colon_grad_u
                                   - ( self.mean_bar_tauxy * self.dy(self.mean_bar_u)
                                     + self.mean_bar_tauyy * self.dy(self.mean_bar_v)
                                     + self.mean_bar_tauyz * self.dy(self.mean_bar_w)))

        # Reynolds average of Reynolds fluctuation pairs
        self.up_up       = self.get_fp_gp(self.mean_bar_u_u, self.mean_bar_u, self.mean_bar_u)
        self.up_vp       = self.get_fp_gp(self.mean_bar_u_v, self.mean_bar_u, self.mean_bar_v)
        self.up_wp       = self.get_fp_gp(self.mean_bar_u_w, self.mean_bar_u, self.mean_bar_w)
        self.vp_vp       = self.get_fp_gp(self.mean_bar_v_v, self.mean_bar_v, self.mean_bar_v)
        self.vp_wp       = self.get_fp_gp(self.mean_bar_v_w, self.mean_bar_v, self.mean_bar_w)
        self.wp_wp       = self.get_fp_gp(self.mean_bar_w_w, self.mean_bar_w, self.mean_bar_w)

        self.Tp_Tp       = self.get_fp_fp(self.mean_bar_T_T    , self.mean_bar_T  )
        self.rhop_rhop   = self.get_fp_fp(self.mean_bar_rho_rho, self.mean_bar_rho)
        self.pp_pp       = self.get_fp_fp(self.mean_bar_p_p    , self.mean_bar_p  )
        self.ap_ap       = self.get_fp_fp(self.mean_bar_a_a    , self.mean_bar_a  )
        self.Mp_Mp       = self.get_fp_fp(self.mean_bar_M_M    , self.mean_bar_M  )
        self.mup_mup     = self.get_fp_fp(self.mean_bar_mu_mu  , self.mean_bar_mu )

        # y * dudy
        self.y_dudy      = self.y * self.dy(self.mean_bar_u)

        # tau
        self.tau         = self.mean_bar_tauxy - self.rho_upp_vpp

        # generalized law of the wall,
        # u = int_0^y \mu dy
        self.generalized_uwall  = np.zeros(self.Ny)
        for j in xrange(1,self.Ny):
            dy = self.y[j] - self.y[j-1]
            self.generalized_uwall [j] = (self.generalized_uwall[j-1]
                    + 0.5 * self.mean_bar_tauxy[0] * (1/self.mean_bar_mu[j] + 1/self.mean_bar_mu[j-1]) * dy)
        del dy

        # Turbulent Mach based on Reynolds and Favre average of velocity fluctuations
        if (self.up_up      is None or 
            self.vp_vp      is None or
            self.wp_wp      is None or
            self.mean_bar_a is None):
            self.Mt = None
        else:
            self.Mt = np.sqrt(self.up_up + self.vp_vp + self.wp_wp) / self.mean_bar_a
    
        self.Mt_Favre = np.sqrt((self.rho_upp_upp + self.rho_vpp_vpp + self.rho_wpp_wpp)  / self.mean_bar_rho) / self.mean_bar_a

        # Prandtl Mixing length
        self.mixing_length =  np.sqrt(np.abs(-self.fav_upp_vpp)) / self.dy(self.mean_bar_u)

        # Turbulent Prandtl
        if (self.rho_Tpp_vpp is None):
            self.Prt = None
        else:
            kappa_u  = self.rho_upp_vpp / self.dy(self.fav_u)
            kappa_T  = self.rho_Tpp_vpp / self.dy(self.fav_T)
            self.Prt = kappa_u / kappa_T
            del kappa_u, kappa_T

        # Total temperature, calorically perfect gas
        if (self.mean_bar_u_u      is None or 
            self.mean_bar_v_v      is None or
            self.mean_bar_w_w      is None or
            self.mean_bar_Cp       is None):
            self.CP_Ttotal = None
        else:
            self.CP_Ttotal = self.mean_bar_T + 0.5 * (self.mean_bar_u_u + self.mean_bar_v_v + self.mean_bar_w_w) / self.mean_bar_Cp

        # Strong Reynolds Analogy(ies)
        # ... evaluated through the convenience quantity G of
        # ... Morinishi etal, JFM, 2004
        def CP_G_function(g=0, h=1):
            if (self.Tp_Tp          is None or
                self.mean_bar_gamma is None or 
                self.mean_bar_M     is None or 
                self.up_up          is None or 
                self.CP_Ttotal      is None):
                return None
            else:
                num        = np.sqrt(self.Tp_Tp) / self.mean_bar_T
                den        = (self.mean_bar_gamma-1) * np.power(self.mean_bar_M,2) * np.sqrt(self.up_up) / self.mean_bar_u
                ghfactor   = h * np.abs(g * self.dy(self.CP_Ttotal) / self.dy(self.mean_bar_T) - 1)
                return num / den * ghfactor
        # -- Reynolds (original)
        self.CP_G_SRA  = CP_G_function()
        # -- Gaviglio (1987)
        self.CP_G_GSRA = CP_G_function(g=1, h=1)
        # -- Rubesin  (1990)
        self.CP_G_RSRA = CP_G_function(g=1, h=1.34)
        # -- Huang    (1990)
        self.CP_G_HSRA = CP_G_function(g=1, h=self.Prt)

        # Turbulent kinetic energy
        self.rhok = 0.5 * (self.rho_upp_upp + self.rho_vpp_vpp + self.rho_wpp_wpp)

        # Kinetic energy budget terms
        # -- mean convection
        self.kbudget_convection = self.fav_v * self.dy(self.rhok) + self.rhok * self.dy(self.fav_v)

        # -- production
        self.kbudget_production = - ( self.rho_upp_vpp * self.dy(self.fav_u)
                                    + self.rho_vpp_vpp * self.dy(self.fav_v)
                                    + self.rho_vpp_wpp * self.dy(self.fav_w) )
        # -- turbulent transport
        self.kbudget_transport = - 0.5 * (self.dy(self.rho_upp_upp_vpp) + self.dy(self.rho_vpp_vpp_vpp) + self.dy(self.rho_vpp_wpp_wpp))

        # -- pressure diffusion
        if self.pp_vpp is None:
            self.kbudget_pressure_diffusion  = None
        else:
            self.kbudget_pressure_diffusion  = - self.dy(self.pp_vpp)

        # -- pressure dilatation
        # -- NOTE: computed explicitly, assuming dudx=dwdz=0
        self.kbudget_pressure_dilatation =  self.mean_bar_p_div_u - self.mean_bar_p * self.dy(self.mean_bar_v)

        # -- viscous diffusion
        self.kbudget_viscous_diffusion   = self.dy(self.taupuppy)

        # -- viscous dissipation
        self.kbudget_viscous_dissipation = self.taup_colon_grad_upp

        # -- compressibility terms
        self.kbudget_compressibility_pressure = - self.vpp  * self.dy(self.mean_bar_p)
        self.kbudget_compressibility_tau      = ( self.upp  * self.dy(self.mean_bar_tauxy)
                                                + self.vpp  * self.dy(self.mean_bar_tauyy)
                                                + self.wpp  * self.dy(self.mean_bar_tauyz))
        self.kbudget_compressibility_tke      = - self.rhok * self.dy(self.fav_v)
        self.kbudget_compressibility          = ( self.kbudget_compressibility_pressure
                                                + self.kbudget_compressibility_tau
                                                + self.kbudget_compressibility_tke     )

        # -- sum (no slow growth source term)
        if self.kbudget_pressure_diffusion is None:
            self.kbudget_sum = None 
        else:
            self.kbudget_sum = ( self.kbudget_convection
                               + self.kbudget_production
                               + self.kbudget_transport
                               + self.kbudget_pressure_diffusion
                               + self.kbudget_pressure_dilatation
                               + self.kbudget_viscous_diffusion
                               - self.kbudget_viscous_dissipation
                               + self.kbudget_compressibility     )


    def get_transformed(self):
        # Van Driest velociself.kbudget_compressibility_tke     ty, integral in u
        self.bar_uVanDriest = np.zeros(self.Ny)
        for j in xrange(1,self.Ny):
            du = self.mean_bar_u[j] - self.mean_bar_u[j-1]
            self.bar_uVanDriest[j]  = (self.bar_uVanDriest[j-1]
                                      + 0.5 * ( np.sqrt(self.mean_bar_rho[j  ] / self.wall_rho)
                                              + np.sqrt(self.mean_bar_rho[j-1] / self.wall_rho)) * du)
        del du

        # Van Driest velocity, integral in u, consider wall injection
        self.bar_uVanDriest_inj = np.zeros(self.Ny)
        for j in xrange(1,self.Ny):
            du = self.mean_bar_u[j] - self.mean_bar_u[j-1]
            self.bar_uVanDriest_inj[j]  = (self.bar_uVanDriest_inj[j-1]
                        + 0.5 * self.u_tau * ( np.sqrt(self.mean_bar_rho[j  ] / (self.wall_V * self.wall_rho * self.fav_u[j  ] + self.wall_rho * np.power(self.u_tau,2)))
                                             + np.sqrt(self.mean_bar_rho[j-1] / (self.wall_V * self.wall_rho * self.fav_u[j-1] + self.wall_rho * np.power(self.u_tau,2)))) * du)

        # Van Driest transformed version of y * dudy
        self.y_duVanDriestdy = self.y * np.sqrt(self.mean_bar_rho / self.wall_rho) * self.dy(self.mean_bar_u)

        # generalized law of the wall, Van Driest transformed
        self.generalized_uwallVanDriest  = np.zeros(self.Ny)
        for j in xrange(1,self.Ny):
            du = self.generalized_uwall[j] - self.generalized_uwall[j-1]
            self.generalized_uwallVanDriest[j] = self.generalized_uwallVanDriest[j-1] + 0.5 * (np.sqrt(self.mean_bar_rho[j] / self.wall_rho) + np.sqrt(self.mean_bar_rho[j-1] / self.wall_rho)) * du
        del du

        return


    def get_scaled_plus(self):

        # Favre averages
        self.fav_u_plus = self.fav_u / self.u_tau
        self.fav_v_plus = self.fav_v / self.u_tau
        self.fav_w_plus = self.fav_w / self.u_tau

        # Favre Reyolds stresses
        self.fav_upp_upp_plus = self.fav_upp_upp / self.u_tau / self.u_tau
        self.fav_upp_vpp_plus = self.fav_upp_vpp / self.u_tau / self.u_tau
        self.fav_upp_wpp_plus = self.fav_upp_wpp / self.u_tau / self.u_tau
        self.fav_vpp_vpp_plus = self.fav_vpp_vpp / self.u_tau / self.u_tau
        self.fav_vpp_wpp_plus = self.fav_vpp_wpp / self.u_tau / self.u_tau
        self.fav_wpp_wpp_plus = self.fav_wpp_wpp / self.u_tau / self.u_tau

        # Incompressible Reyolds stresses
        if self.up_up      is not None:
            self.up_up_plus = self.up_up / self.u_tau / self.u_tau
        if self.up_vp      is not None:
            self.up_vp_plus = self.up_vp / self.u_tau / self.u_tau
        if self.up_wp      is not None:
            self.up_wp_plus = self.up_wp / self.u_tau / self.u_tau
        if self.vp_vp      is not None:
            self.vp_vp_plus = self.vp_vp / self.u_tau / self.u_tau
        if self.vp_wp      is not None:
            self.vp_wp_plus = self.vp_wp / self.u_tau / self.u_tau
        if self.wp_wp      is not None:
            self.wp_wp_plus = self.wp_wp / self.u_tau / self.u_tau

        # y * dudy
        self.y_dudy_plus = self.y_dudy / self.u_tau

        # Van Driest
        self.bar_uVanDriest_plus = self.bar_uVanDriest / self.u_tau

        # Van Driest, wall-injection corrected
        self.bar_uVanDriest_inj_plus  = self.bar_uVanDriest_inj / self.u_tau

        # Van Driest transformed version of y * dudy
        self.y_duVanDriestdy_plus = self.y_duVanDriestdy / self.u_tau

        # Generalized law of the wall, Van Driest transformed
        self.generalized_uwallVanDriest_plus = self.generalized_uwallVanDriest / self.u_tau

        # Mixing length
        self.mixing_length_plus = self.mixing_length / self.delta_nu

        # Van Driest velocity, injection corrected, composite layer
        def integrand_composite_plus(ku=0.41, A_plus=25.51, lmix='model'):
            irho       = self.mean_bar_rho/self.wall_rho
            imu        = self.mean_bar_mu/self.wall_mu
            f_wallV    = (self.wall_V_plus * self.fav_u/self.u_tau + 1)
            if lmix == "dns":
                l_plus = self.mixing_length_plus
            elif lmix == "model":
                l_plus = self.mixing_length_model_plus(rho=irho, mu=imu)
            l_plus_inc = self.mixing_length_model_plus()
            num    = 1 + np.sqrt(4 * np.power(l_plus/imu, 2) * irho * f_wallV + 1)
            den    = f_wallV / imu * (1 + np.sqrt(4 * np.power(l_plus_inc, 2) + 1))
            ret    = num/den
            ret[0] = 0
            return ret

        self.bar_uVanDriest_inj_composite_plus = np.zeros(self.Ny)
        integrand = integrand_composite_plus()
        for j in xrange(1,self.Ny):
            du = (self.mean_bar_u[j] - self.mean_bar_u[j-1]) / self.u_tau
            self.bar_uVanDriest_inj_composite_plus[j]  = (self.bar_uVanDriest_inj_composite_plus[j-1]
                        + 0.5 * ( integrand[j-1] + integrand[j  ]) * du)
        del integrand
        del du

        self.bar_uVanDriest_inj_dns_plus = np.zeros(self.Ny)
        integrand = integrand_composite_plus(lmix="dns")
        for j in xrange(1,self.Ny):
            du = (self.mean_bar_u[j] - self.mean_bar_u[j-1]) / self.u_tau
            self.bar_uVanDriest_inj_dns_plus[j]  = (self.bar_uVanDriest_inj_dns_plus[j-1]
                        + 0.5 * ( integrand[j-1] + integrand[j  ]) * du)
        del integrand
        del du

        return


    # Semi-local scaling
    def get_scaled_sl(self):

        # Favre averages
        self.fav_u_sl = self.fav_u / self.u_tau_sl
        self.fav_v_sl = self.fav_v / self.u_tau_sl
        self.fav_w_sl = self.fav_w / self.u_tau_sl

        # Favre Reyolds stresses
        self.fav_upp_upp_sl = self.fav_upp_upp / self.u_tau_sl / self.u_tau_sl
        self.fav_upp_vpp_sl = self.fav_upp_vpp / self.u_tau_sl / self.u_tau_sl
        self.fav_upp_wpp_sl = self.fav_upp_wpp / self.u_tau_sl / self.u_tau_sl
        self.fav_vpp_vpp_sl = self.fav_vpp_vpp / self.u_tau_sl / self.u_tau_sl
        self.fav_vpp_wpp_sl = self.fav_vpp_wpp / self.u_tau_sl / self.u_tau_sl
        self.fav_wpp_wpp_sl = self.fav_wpp_wpp / self.u_tau_sl / self.u_tau_sl

        return

    def get_bl(self):
        # Input parameters
        self.wall_T =  self.mean_bar_T[0]
        self.wall_U =  self.mean_bar_u[0]
        self.wall_V =  self.mean_bar_v[0]
        self.wall_W =  self.mean_bar_w[0]

        self.inf_T  =  self.mean_bar_T[-1]
        self.inf_U  =  self.mean_bar_u[-1]
        self.inf_V  =  self.mean_bar_v[-1]
        self.inf_W  =  self.mean_bar_w[-1]

        # Grid parameters
        self.y1     =  self.y[1]
        self.y1b    =  self.yb[1]

        # Wall values
        self.wall_p     =  self.get_wall( self.mean_bar_p            )
        self.wall_a     =  self.get_wall( self.mean_bar_a            )
        self.wall_M     =  self.get_wall( self.mean_bar_M            )
        self.wall_Cp    =  self.get_wall( self.mean_bar_Cp           )
        self.wall_Cv    =  self.get_wall( self.mean_bar_Cv           ) 
        self.wall_rho   =  self.get_wall( self.mean_bar_rho          ) 
        self.wall_mu    =  self.get_wall( self.mean_bar_mu           ) 
        self.wall_nu    =  self.get_wall( self.mean_bar_nu           ) 
        self.wall_D0    =  self.get_wall( self.mean_bar_D0           ) 
        self.wall_Cp    =  self.get_wall( self.mean_bar_Cp           ) 
        self.wall_q     =  self.get_wall(-self.mean_bar_kappa_grady_T)
        self.wall_H     =  self.get_wall( self.fav_H                 )
        self.wall_dudy  =  self.get_wall( self.dy(self.mean_bar_u)   )

        # Edge values
        self.inf_p      =  self.get_inf( self.mean_bar_p            )
        self.inf_a      =  self.get_inf( self.mean_bar_a            )
        self.inf_M      =  self.get_inf( self.mean_bar_M            )
        self.inf_Cp     =  self.get_inf( self.mean_bar_Cp           )
        self.inf_Cv     =  self.get_inf( self.mean_bar_Cv           )
        self.inf_rho    =  self.get_inf( self.mean_bar_rho          )
        self.inf_mu     =  self.get_inf( self.mean_bar_mu           )
        self.inf_nu     =  self.get_inf( self.mean_bar_nu           )
        self.inf_D0     =  self.get_inf( self.mean_bar_D0           )
        self.inf_Cp     =  self.get_inf( self.mean_bar_Cp           )
        self.inf_q      =  self.get_inf(-self.mean_bar_kappa_grady_T)
        self.inf_H      =  self.get_inf( self.fav_H                 )
        self.inf_dudy   =  self.get_inf( self.dy(self.mean_bar_u)   )

        # Compute delta (BL thickness)
        jdelta = self.Ny-1
        for j, u in enumerate(self.mean_bar_u):
            if (u > 0.99*self.inf_U):
                jdelta = j-1
                break
        # ... get delta by interpolation
        frac  = (0.99*self.inf_U - self.mean_bar_u[jdelta]) / (self.mean_bar_u[jdelta+1] - self.mean_bar_u[jdelta])
        self.delta = (self.y[jdelta+1] - self.y[jdelta]) * frac + self.y[jdelta]
        del jdelta
        del frac

        # TODO: Here the value at the edge is defined as the value at "infinity"
        #       This needs to be generalized for baseflow computations
        self.edge_rho   = self.inf_rho
        self.edge_rhoU  = self.inf_rho * self.inf_U
        self.edge_U     = self.inf_U
        self.edge_mu    = self.inf_mu

        # Compute delta_star (displacement thickness)
        self.delta_star = 0
        for j, ru in enumerate(self.mean_bar_rho_u):
            self.delta_star += (1 - ru / self.edge_rhoU) * self.iw[j]

        # Compute theta (momentum thickness)
        self.theta = 0
        for j, (ru, u) in enumerate(zip(self.mean_bar_rho_u, self.mean_bar_u)):
            self.theta += ru / self.edge_rhoU * (1 - u / self.edge_U) * self.iw[j]

        # Computed parameters
        self.Re_delta_star    = self.edge_rho * self.edge_U * self.delta_star / self.edge_mu
        self.Re_theta         = self.edge_rho * self.edge_U * self.theta      / self.edge_mu
        self.Re_delta         = self.edge_rho * self.edge_U * self.delta      / self.edge_mu
        self.H1               = self.delta_star / self.theta
        self.H2               = self.delta      / self.theta
        self.wall_tau         = self.wall_mu  * self.wall_dudy
        self.u_tau            = np.sqrt(self.wall_tau / self.wall_rho)
        self.delta_nu         = self.wall_mu    / self.wall_rho / self.u_tau
        self.Re_tau           = self.wall_rho * self.u_tau * self.delta / self.wall_mu
        self.turnover_time    = self.delta      / self.u_tau
        self.flowthrough_time = self.Lx         / self.edge_U
        self.skin_friction    = 2.0 * np.power((self.u_tau/self.edge_U),2) * self.wall_rho/self.edge_rho
        if self.inf_Cp is None:
            self.T_total          = None
        else:
            self.T_total          = self.inf_H      / self.inf_Cp

        self.M_tau            = self.u_tau      / self.wall_a
        if self.wall_Cp is None:
            self.T_tau            = None
        else:
            self.T_tau            = self.wall_q     / (self.wall_rho * self.wall_Cp * self.u_tau)

        if self.wall_Cp is None:
            self.B_q              = None
        else:
            self.B_q              = self.wall_q     / (self.wall_rho * self.wall_Cp * self.u_tau * self.wall_T)

        self.F_injection      = self.wall_rho * self.wall_V / (self.inf_rho * self.inf_U)

        # Box size and resolution parameters
        self.Lx_over_delta    = self.Lx  / self.delta
        self.Ly_over_delta    = self.Ly  / self.delta
        self.Lz_over_delta    = self.Lz  / self.delta
        self.y1_plus          = self.y1  / self.delta_nu
        self.y1b_plus         = self.y1b / self.delta_nu
        self.Dx               = self.Lx  / self.Nx
        self.Dz               = self.Lz  / self.Nz
        self.Dx_plus          = self.Dx  / self.delta_nu
        self.Dz_plus          = self.Dz  / self.delta_nu

        # Coordinate scalings
        # delta
        self.y_delta  = self.y  / self.delta

        # ... wall units
        self.y_plus   = self.y  / self.delta_nu
        self.yb_plus  = self.yb / self.delta_nu

        # ... semilocal viscous friction velocity, length scale, wall units, Re_tau
        self.u_tau_sl    = np.sqrt(self.wall_tau / self.mean_bar_rho)
        self.delta_nu_sl = self.mean_bar_mu / self.mean_bar_rho / self.u_tau_sl
        self.y_sl        = self.y      / self.delta_nu_sl
        self.Re_tau_sl   = self.delta  / self.delta_nu_sl

        # wall unit quantities (plus)
        self.wall_V_plus = self.wall_V / self.u_tau

        return


    def print_bl(self):
        print "rho_wall           = ", self.wall_rho
        print "p_wall             = ", self.wall_p
        print "a_wall             = ", self.wall_a
        print "mu_wall            = ", self.wall_mu
        print "mu_inf             = ", self.inf_mu
        print "nu_wall            = ", self.wall_nu
        print "nu_inf             = ", self.inf_nu
        print "Cp_wall            = ", self.wall_Cp
        print "Cp_inf             = ", self.inf_Cp
        print "V_wall             = ", self.wall_V
        print "V_wall_plus        = ", self.wall_V_plus
        print "V_wall/U_edge      = ", self.wall_V / self.edge_U
        print "U_edge             = ", self.edge_U
        print "U_inf              = ", self.inf_U
        print "V_inf              = ", self.inf_V
        print "W_inf              = ", self.inf_W
        print "T_inf              = ", self.inf_T
        print "T_total            = ", self.T_total

        if self.T_total is not None:
            print "T_wall/T_total     = ", self.wall_T / self.T_total

        print "M_inf              = ", self.inf_M
        print "rho_inf            = ", self.inf_rho
        print "p_inf              = ", self.inf_p
        print "a_inf              = ", self.inf_a
        print "q_wall             = ", self.wall_q

        print "delta_star         = ", self.delta_star
        print "theta              = ", self.theta
        print "delta              = ", self.delta
        print "H1                 = ", self.H1
        print "H2                 = ", self.H2

        print "dudy_wall          = ", self.wall_dudy
        print "tau_wall           = ", self.wall_tau
        print "u_tau              = ", self.u_tau
        print "delta_nu           = ", self.delta_nu
        print "Cf (skin friction) = ", self.skin_friction
        print "M_tau              = ", self.M_tau
        print "T_tau              = ", self.T_tau
        print "B_q                = ", self.B_q
        print "F_injection        = ", self.F_injection

        print "Re_tau             = ", self.Re_tau
        print "Re_delta_star      = ", self.Re_delta_star
        print "Re_theta           = ", self.Re_theta
        print "Re_delta           = ", self.Re_delta

        print "y1b_plus           = ", self.y1b_plus
        print "Dx_plus            = ", self.Dx_plus
        print "Dz_plus            = ", self.Dz_plus
        print "Lx_over_delta      = ", self.Lx_over_delta
        print "Ly_over_delta      = ", self.Ly_over_delta
        print "Lz_over_delta      = ", self.Lz_over_delta
        print "nyb_below_5plus    = ", self.get_points_below( 5, self.yb_plus)
        print "nyb_below_10plus   = ", self.get_points_below(10, self.yb_plus)
        print "nyb_below_delta    = ", self.get_points_below(self.delta, self.yb)
        print "turnover_time      = ", self.turnover_time
        print "flowthrough_time   = ", self.flowthrough_time


    def get_baseflow(self):
        # TODO, compute baseflow info
        return

    def get_j_for_y_plus(self, y_plus):
        for j, y in enumerate(self.y_plus):
            if y > y_plus:
                break
        return j-1

    def dy(self, var):
        return np.ravel(np.dot(var, self.invD0T_D1T))

    def ddy(self, var):
        return np.ravel(np.dot(var, self.invD0T_D2T))

    def get_u_loglaw_plus(self, kappa_vonKarman=0.4, B=4.7):
        return (1.0/kappa_vonKarman) * np.log(self.y_plus) + B

    def mixing_length_model_plus(self
            , rho    = 1
            , mu     = 1
            , wall_V = 0
            , ku=0.4, A_plus=26):
        fav_u   = self.fav_u_plus
        f_wallV = (wall_V * fav_u + 1)
        return ku * self.y_plus * (1-np.exp(-self.y_plus*np.sqrt(rho*f_wallV)/mu/A_plus))

    # Some convenience functions:
    def get_wall(self, f):
        if f is None:
            return None
        else:
            return f[0]

    def get_inf(self, f):
        if f is None:
            return None
        else:
            return f[-1]

    def get_fav(self, rho_f, rho):
        if rho_f is None:
            return None
        else:
            return rho_f / rho

    def get_rho_fpp_gpp(self, rho_f_g, rho, fav_f, fav_g):
        if (rho_f_g is None or
            fav_f   is None or
            fav_g   is None):
            return None
        else:
            return rho_f_g - rho * fav_f * fav_g

    def get_rho_fpp_fpp_gpp(self, rho_f_f_g, rho_f_f, rho_f_g, rho, fav_f, fav_g):
        return rho_f_f_g - 2 * rho_f_g * fav_f - rho_f_f* fav_g + 2 * rho * fav_g * fav_f * fav_f

    def get_rho_fpp_gpp_hpp(self, rho_f_g_h, rho_f_g, rho_g_h, rho_f_h, rho, fav_f, fav_g, fav_h):
        return (rho_f_g_h - (rho_f_g - rho * fav_f * fav_g) * fav_h
                          - (rho_g_h - rho * fav_g * fav_h) * fav_f
                          - (rho_f_h - rho * fav_f * fav_h) * fav_g
                          -            rho * fav_f * fav_g  * fav_h)

    def get_fp_gp(self, f_g, bar_f, bar_g):
        if (f_g   is None or
            bar_f is None or
            bar_g is None):
            return None
        else:
            return f_g - bar_f * bar_g

    def get_fp_fp(self, f_f, bar_f):
        return self.get_fp_gp(f_f, bar_f, bar_f)

    def get_fp_gpp(self, f_g, bar_f, bar_g):
        return self.get_fp_gp(f_g, bar_f, bar_g)

    def get_fpp_gpp(self, f_g, bar_f, bar_g, fav_f, fav_g):
        return f_g + fav_f*fav_g - fav_f*bar_g - bar_f*fav_g


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Parse and check incoming command line arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
        if len(args) < 1:
            print >>sys.stderr, "Incorrect number of arguments.  See --help."
            return 2
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

    # Input dict defaults
    input = {}
    input.update(dict(
        hdf5files = args,
    ))

    # Process all files
    sum = summary(**input)

    # save summarized data
    # and delete the object
    sum.save()
    del sum

    # Now try to load the data
    print "-"*80
    input = {}
    input.update(dict(
        hdf5files = ['summary.h5'],
        #conserved_only = conserved_only
    ))

    # Initialize summary by loading summary file
    sum = summary(**input)
    sum.print_bl()
    del sum


if __name__ == "__main__":
    sys.exit(main())
