import numpy as np
import h5py
import completeness
import dust_attenuation


class snapshot:
    '''
    This classes computes key quantities for
    a SP snapshot from the halo model.
    '''

    def __init__(self, SP_file_name, redshift):
        '''
        Load SP file.

        Parameters
        ----------
        SP_file_name : str
          The filename of SP file (hdf5 file).

        redshift : float
          The redshift of the snapshot considered.

        '''
        self.data = h5py.File(SP_file_name, 'r')
        self.redshift = redshift
        #self.stellar_mass = np.array([])
        print('snapshot has been successfully loaded')

    def get_halo_mass(self, exclude_contam_halos=True):
        '''
        Returns halo mass.

        Parameters
        ----------
        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        halo_mass: float
          Halo mass in Msun.

        '''
        halo_mass = self.data['DM/DM_M'][:]
        if exclude_contam_halos:
            idx = ~self.get_contaminated_halos()
        else:
            idx = (np.ones(len(halo_mass)) == 1.0)
        return(halo_mass[idx])

    def get_Z_hist(self, exclude_contam_halos=True):
        '''
        Returns metallicity.

        Parameters
        ----------
        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        Z_hist: float
          Z_hist = M_met/M_gas for all times.

        '''
        Z_hist = self.data['SFH/SFH_Z'][:]
        if exclude_contam_halos:
            idx = ~self.get_contaminated_halos()
        else:
            idx = (np.ones(Z_hist.shape[0]) == 1.0)
        return(Z_hist[idx])

    def get_Z(self, exclude_contam_halos=True):
        '''
        Returns metallicity.

        Parameters
        ----------
        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        Z: float
          Z = M_met/M_gas.

        '''
        Z = self.data['SFH/SFH_Z'][:, -1]
        if exclude_contam_halos:
            idx = ~self.get_contaminated_halos()
        else:
            idx = (np.ones(len(Z)) == 1.0)
        return(Z[idx])

    def get_logOH(self, exclude_contam_halos=True):
        '''
        Returns metallicty log(O/H), assuming
        Z_sun = 0.0143
        logOH (+ 12) = 8.69
        from Asplund et al. 2009, 2012.

        Parameters
        ----------
        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        logOH: float
          Metallicity in units of log O/H + 12

        '''
        Z_sun = 0.0143
        conversion = 8.69-np.log10(Z_sun)
        logOH = np.log10(self.get_Z(exclude_contam_halos=exclude_contam_halos))+conversion
        return(logOH)

    def get_spectrum(self, idx=None, SP_param_nr='4', exclude_contam_halos=True):
        '''
        Returns emission line luminosity.

        Parameters
        ----------
        idx : int
          Specifices the index of the galaxy.

        SP_param_nr : str
          Stellar population parameter number that
          defines choice of initial mass function,
          metallicity.
          fiducial: '4'

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        wavelength : np.array
          Wavelength in Angstrom.

        spectrum :  np.array
          Spectral energy density in Lsun/Hz.

        '''
        if exclude_contam_halos:
            idx_new = ~self.get_contaminated_halos() & idx
        else:
            idx_new = (np.ones(len(passband_lum)) == 1.0) & idx
        wavelength = self.data['SP/spec/wavelength'][:]
        spectrum = self.data['SP/spec/luminosity_' + SP_param_nr][idx_new]
        return(wavelength, spectrum)

    def get_EmL_lum(self, emission_line='L_Ha', SP_param_nr='4', exclude_contam_halos=True):
        '''
        Returns emission line luminosity.

        Parameters
        ----------
        emission_line : str
          Emission line; possible choices include:
          L_Lya, L_HeII, L_OIII_L1, L_OIII_L2, L_CIII_1,
          L_CIII_2, L_CIV, L_OII, L_Hb, L_OIII, L_Ha,
          L_NII, L_SII_1, L_SII_2

        SP_param_nr : str
          Stellar population parameter number that
          defines choice of initial mass function,
          metallicity.
          fiducial: '0'

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        line_luminosity: np.array
          Emission line luminosity in erg/s.

        '''
        idx_EL = (self.data['SP/EmL'].attrs['EL_info'] == emission_line)
        line_luminosity = self.data['SP/EmL/luminosity_' + SP_param_nr][:, idx_EL]
        if exclude_contam_halos:
            idx = ~self.get_contaminated_halos()
        else:
            idx = (np.ones(len(line_luminosity)) == 1.0)
        return(line_luminosity[idx])

    def get_band_lum(self, passband='i1500', SP_param_nr='4', add_dust=False, exclude_contam_halos=True):
        '''
        Returns luminosity in a certain passband.

        Parameters
        ----------
        passband : str
          Passband; possible choices include:
          i1500, i2300, i2800, v, u, 2mass_j

        SP_param_nr : str
          Stellar population parameter number that
          defines choice of initial mass function,
          metallicity.
          fiducial: '4'

        add_dust : bool
          Add dust to UV (ONLY UV 1500)
          fiducial: False

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        lum_passband: np.array
          Luminosity in a certain passband in erg/s.

        '''
        wave_dict = [('i1500', 1500.0), ('u', 3650.0), ('v', 5500.0), ('2mass_j', 12200.0)]
        wave_dict = dict(wave_dict)
        idx_passband = (self.data['SP/FilL'].attrs['FL_info'] == passband)
        passband_lum = self.data['SP/FilL/luminosity_' + SP_param_nr][:, idx_passband].flatten()
        if (add_dust):
            mag_list_UV = -48.6-2.5*np.log10(self.data['SP/FilL/luminosity_' + SP_param_nr][:, (self.data['SP/FilL'].attrs['FL_info'] == 'i1500')].flatten()/(4*np.pi*(3.086e+19)**2))
            AUV_list = dust_attenuation.add_dust_attenuation(mag_list_UV, self.redshift)-(mag_list_UV)
            mag_list = -48.6-2.5*np.log10(passband_lum/(4*np.pi*(3.086e+19)**2))
            mag_list_d = mag_list+AUV_list/dust_attenuation.calzetti(1500.0, tau_v=1, R_v=4.05)*dust_attenuation.calzetti(wave_dict[passband], tau_v=1, R_v=4.05)
            passband_lum = (4*np.pi*(3.086e+19)**2)*10.0**(-1.0/2.5*(mag_list_d+48.6))
        if exclude_contam_halos:
            idx = ~self.get_contaminated_halos()
        else:
            idx = (np.ones(len(passband_lum)) == 1.0)
        return(passband_lum[idx])

    def get_stellar_mass(self, integral_of_SFR=True, SP_param_nr='4', Mseed=0.0, R=0.0, exclude_contam_halos=True):
        '''
        Computes total stellar mass; after returning
        fraction R to the interstellar medium.

        Parameters
        ----------
        integral_of_SFR : bool
          Specifies whether we define the stellar mass to be
          the integral of the past SFR or to be the actual
          mass in stars and remnants (i.e. subtracting return).

        Mseed : float
          The seed mass at BB.
          fiducial: 0

        R : float
          Fraction of the mass that is converted into stars,
          as measured by the SFR, is promptly (we will assume
          instantaneously) returned to the interstellar medium.
          The remaining fraction (1-R) stays in form of
          long-lived stars.
          see: http://adsabs.harvard.edu/abs/2016MNRAS.455.4183V
          fiducial: 0.29 (Salpeter)

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        stellar_mass: float
          Stellar mass in Msun (remaining fraction 1-R).

        '''
        if integral_of_SFR:
            # set parameters
            print 'calculate integral of SFR with following parameters...'
            print 'seed mass Mseed = ', Mseed
            print 'and return R = ', R
            # load SFH for given time interval
            SFR_list = self.data['SFH/SFH_SFR'][:]
            time_list = self.data['SFH/SFH_time'][:]
            #SFH = np.repeat(SFR_list, 2, axis=1)[:, 2:]
            #time = np.append(0.0, np.repeat(time_list[1:]+0.5*np.diff(time_list), 2))[:-1]
            # compute mass
            stellar_mass = Mseed+(1.0-R)*np.trapz(SFR_list, time_list*10**6)
        else:
            print 'get stellar mass in stars and remnants...'
            idx_M = (self.data['SP/FilL'].attrs['FL_info'] == 'stellar_mass')
            stellar_mass = self.data['SP/FilL/luminosity_' + SP_param_nr][:, idx_M].flatten()
        if exclude_contam_halos:
            idx = ~self.get_contaminated_halos()
        else:
            idx = (np.ones(len(stellar_mass)) == 1.0)
        return(stellar_mass[idx])

    def get_SFR(self, time_interval, exclude_contam_halos=True):
        '''
        Computes average SFR over certain time interval.

        Parameters
        ----------
        time_interval : float
          Time in which SFR is comupted, in Myr.

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        SFR: float
          SFR in Msun/yr.

        '''
        print 'compute average SFR over time interval (Myr): ', time_interval
        # load SFH for given time interval
        idx_time_interval = (self.data['SFH/SFH_time'][:] >= (self.data['SFH/SFH_time'][:][-1]-time_interval))
        SFR_list = self.data['SFH/SFH_SFR'][:, idx_time_interval]
        time_list = self.data['SFH/SFH_time'][idx_time_interval]
        time_interval = (time_list[-1]-time_list[0])*10**6  # in yrs
        # compute mass
        mass_tot_list = np.trapz(SFR_list, time_list*10**6)
        # compute average SFR
        SFR = mass_tot_list/time_interval
        if exclude_contam_halos:
            idx = ~self.get_contaminated_halos()
        else:
            idx = (np.ones(len(SFR)) == 1.0)
        return(SFR[idx])

    def get_SFH(self, exclude_contam_halos=True):
        '''
        Obtain SFHs of all galaxies.

        Parameters
        ----------
        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        SFH_time: np.array
          Time since Big Bang in Myr.

        SFH_SFR: np.array x np.array
          SFR in Msun/yr as a function of time for each galaxy.

        '''
        SFH_time = self.data['SFH/SFH_time'][:]
        SFH_SFR = self.data['SFH/SFH_SFR'][:]
        if exclude_contam_halos:
            idx = ~self.get_contaminated_halos()
        else:
            idx = (np.ones(SFH_SFR.shape[0]) == 1.0)
        return(SFH_time, SFH_SFR[idx])


    def get_SFH_int(self, exclude_contam_halos=True):
        '''
        Obtain integrated SFHs (i.e. stellar masses)
        of all galaxies.

        Parameters
        ----------
        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        SFH_time: np.array
          Time since Big Bang in Myr.

        SFH_SFR: np.array x np.array
          SFR in Msun/yr as a function of time for each galaxy.

        '''
        SFH_time, SFH_SFR_matrix = self.get_SFH(exclude_contam_halos=exclude_contam_halos)
        for ii in range(SFH_SFR_matrix.shape[1]):
            if (ii == 0):
                SFH_int = np.trapz(SFH_SFR_matrix[:, :ii],  SFH_time[:ii]*10**6)
            else:
                SFH_int = np.vstack([SFH_int, np.trapz(SFH_SFR_matrix[:, :ii],  SFH_time[:ii]*10**6)])
        return(SFH_time, SFH_int.T)


    def get_DM_accretion(self, exclude_contam_halos=True):
        '''
        Obtain dark matter accretion history of all galaxies.

        Parameters
        ----------
        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        Returns
        -------
        DM_time: np.array
          Time since Big Bang in Myr.

        DM_Mt: np.array x np.array
          Halo mass as a function of time for each galaxy.

        '''
        DM_time = self.data['DM/DM_time'][:]
        DM_Mt = self.data['DM/DM_Mt'][:]
        if exclude_contam_halos:
            idx = ~self.get_contaminated_halos()
        else:
            idx = (np.ones(DM_Mt.shape[0]) == 1.0)
        return(DM_time, DM_Mt[idx])

    def get_contaminated_halos(self):
        '''
        Obtain dark matter accretion history of all galaxies.

        Parameters
        ----------
        None

        Returns
        -------
        contam_flag: np.array
          Returns array with True for contaminated halos,
          and False for halos that are not contaminated.

        '''
        contam_flag = (self.data['DM/DM_cont'][:] == 1.0)
        return(contam_flag)

    def compute_SMF(self, integral_of_SFR=True, SP_param_nr='4', Mseed=0.0, R=0.0, volume_box=100.0**3, bin_size=0.1, cumulative=False, exclude_contam_halos=True, completeness_correction=False, completeness_correction_type='parameterized'):
        '''
        Computes stellar mass function (SMF).

        Parameters
        ----------
        integral_of_SFR : bool
          Specifies whether we define the stellar mass to be
          the integral of the past SFR or to be the actual
          mass in stars and remnants (i.e. subtracting return).

        Mseed : float
          The seed mass at BB.
          fiducial: 0

        R : float
          Fraction of the mass that is converted into stars,
          as measured by the SFR, is promptly (we will assume
          instantaneously) returned to the interstellar medium.
          The remaining fraction (1-R) stays in form of
          long-lived stars.
          see: http://adsabs.harvard.edu/abs/2016MNRAS.455.4183V
          fiducial: 0.29 (Salpeter)

        volume_box : float
          Volume of the box in Mpc^3.
          fiducial: 100.0^3

        bin_size : float
          Size of stellar mass bins (in log units).
          fiducial: 0.1

        cumulative : bool
          Compute cumulative, n(>M), or non-cumulative,
          dn/dlogM, stellar mass function.
          fiducial: False

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        completeness_correction : bool
          Apply completeness correction.
          fiducial: False

        completeness_correction_type : str
          Type of completeness correction:
            'parameterized' / 'numerical'
          fiducial: 'parameterized'

        Returns
        -------
        stellar_mass_bins : np.array
          Stellar mass in log units Msun.

        phi : np.array
          Number density of objects, units depend if
          cumulative or not.

        '''
        M_list = np.log10(self.get_stellar_mass(integral_of_SFR=integral_of_SFR, SP_param_nr=SP_param_nr, Mseed=Mseed, R=R, exclude_contam_halos=exclude_contam_halos))
        idx_good = np.isfinite(M_list)
        M_list = M_list[idx_good]
        M_bins = np.arange(np.max([6.0, np.min(M_list)]), np.min([11.5, np.max(M_list)]), bin_size)
        M_bins_center = M_bins[:-1] + 0.5*np.diff(M_bins)
        if completeness_correction:
            if (completeness_correction_type == 'numerical'):
                weights = 10**completeness.get_completeness_correction_numerical(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx_good]), self.redshift)
            elif (completeness_correction_type == 'parameterized'):
                weights = 10**completeness.get_completeness_correction_parametrized(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx_good]), self.redshift)
        else:
            weights = None
        hist, bin_edges = np.histogram(M_list, bins=M_bins, weights=weights)
        if cumulative:
            M_Fctcum = np.cumsum(hist[::-1])[::-1]/volume_box
            return(M_bins_center, M_Fctcum)
        else:
            M_Fct = hist/(np.diff(M_bins)*volume_box)
            return(M_bins_center, M_Fct)

    def compute_SFR_Fct(self, time_interval=200.0, volume_box=100.0**3, bin_size=0.1, cumulative=False, exclude_contam_halos=True, completeness_correction=False, completeness_correction_type='parameterized'):
        '''
        Computes stellar mass function (SMF).

        Parameters
        ----------
        time_interval : float
          Time in which SFR is comupted, in Myr.
          fiducial: 200

        volume_box : float
          Volume of the box in Mpc^3.
          fiducial: 100.0^3

        bin_size : float
          Size of stellar mass bins (in log units).
          fiducial: 0.1

        cumulative : bool
          Compute cumulative, n(>M), or non-cumulative,
          dn/dlogM, stellar mass function.
          fiducial: False

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        completeness_correction : bool
          Apply completeness correction.
          fiducial: False

        completeness_correction_type : str
          Type of completeness correction:
            'parameterized' / 'numerical'
          fiducial: 'parameterized'

        Returns
        -------
        SFR_bins : np.array
          SFR in log units Msun/yr.

        phi : np.array
          Number density of objects, units depend if
          cumulative or not.

        '''
        SFR_list = np.log10(self.get_SFR(time_interval, exclude_contam_halos=exclude_contam_halos))
        idx_good = np.isfinite(SFR_list)
        SFR_list = SFR_list[idx_good]
        SFR_bins = np.arange(np.max([-2.0, np.min(SFR_list)]), np.min([3.0, np.max(SFR_list)]), bin_size)
        SFR_bins_center = SFR_bins[:-1] + 0.5*np.diff(SFR_bins)
        if completeness_correction:
            if (completeness_correction_type == 'numerical'):
                weights = 10**completeness.get_completeness_correction_numerical(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx_good]), self.redshift)
            elif (completeness_correction_type == 'parameterized'):
                weights = 10**completeness.get_completeness_correction_parametrized(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx_good]), self.redshift)
        else:
            weights = None
        hist, bin_edges = np.histogram(SFR_list, bins=SFR_bins, weights=weights)
        if cumulative:
            SFR_Fctcum = np.cumsum(hist[::-1])[::-1]/volume_box
            return(SFR_bins_center, SFR_Fctcum)
        else:
            SFR_Fct = hist/(np.diff(SFR_bins)*volume_box)
            return(SFR_bins_center, SFR_Fct)

    def compute_UVLF(self, SP_param_nr='4', volume_box=100.0**3, bin_size=0.25, cumulative=False, add_dust=False, exclude_contam_halos=True, completeness_correction=False, completeness_correction_type='parameterized'):
        '''
        Computes UV luminosity function (UVLF).

        Parameters
        ----------
        SP_param_nr : str
          Stellar population parameter number that
          defines choice of initial mass function,
          metallicity.
          fiducial: '0'

        volume_box : float
          Volume of the box in Mpc^3.
          fiducial: 100.0^3

        bin_size : float
          Size of UV mag bins.
          fiducial: 0.25

        cumulative : bool
          Compute cumulative, n(>Muv), or non-cumulative,
          dn/dMuv, UV luminosity function.
          fiducial: False

        add_dust : bool
          Add dust attenuation to UV.
          fiducial: False

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        completeness_correction : bool
          Apply completeness correction.
          fiducial: False

        completeness_correction_type : str
          Type of completeness correction:
            'parameterized' / 'numerical'
          fiducial: 'parameterized'

        Returns
        -------
        mag_bins : np.array
          Absolute UV magnitude (1500 A).

        phi : np.array
          Number density of objects, units depend if
          cumulative or not.

        '''
        mag_list = -48.6-2.5*np.log10(self.get_band_lum(passband='i1500', SP_param_nr=SP_param_nr, exclude_contam_halos=exclude_contam_halos)/(4*np.pi*(3.086e+19)**2))
        idx_good = np.isfinite(mag_list)
        mag_list = mag_list[idx_good]
        mag_bins = np.arange(np.max([-25.0, np.min(mag_list)]), np.min([-12.5, np.max(mag_list)]), bin_size)
        mag_bins_center = mag_bins[:-1] + 0.5*np.diff(mag_bins)
        if completeness_correction:
            if (completeness_correction_type == 'numerical'):
                weights = 10**completeness.get_completeness_correction_numerical(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx_good]), self.redshift)
            elif (completeness_correction_type == 'parameterized'):
                weights = 10**completeness.get_completeness_correction_parametrized(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx_good]), self.redshift)
        else:
            weights = None
        hist, bin_edges = np.histogram(mag_list, bins=mag_bins, weights=weights)
        if add_dust:
            mag_bins_center_d = dust_attenuation.add_dust_attenuation(mag_bins_center, self.redshift)
        else:
            mag_bins_center_d = mag_bins_center
        if cumulative:
            LFcum = np.cumsum(hist)/volume_box
            return(mag_bins_center_d, LFcum)
        else:
            LF = hist/(np.diff(mag_bins)*volume_box)
            return(mag_bins_center_d, LF)

    def compute_EmLLF(self, emission_line='L_Ha', SP_param_nr='4', volume_box=100.0**3, bin_size=0.25, cumulative=False, add_dust=False, exclude_contam_halos=True, completeness_correction=False, completeness_correction_type='parameterized'):
        '''
        Computes emission line luminosity function (EmLLF).

        Parameters
        ----------
        emission_line : str
          Emission line; possible choices include:
          L_Lya, L_HeII, L_OIII_L1, L_OIII_L2, L_CIII_1,
          L_CIII_2, L_CIV, L_OII, L_Hb, L_OIII, L_Ha,
          L_NII, L_SII_1, L_SII_2

        SP_param_nr : str
          Stellar population parameter number that
          defines choice of initial mass function,
          metallicity.
          fiducial: '0'

        volume_box : float
          Volume of the box in Mpc^3.
          fiducial: 100.0^3

        bin_size : float
          Size of UV mag bins.
          fiducial: 0.25

        cumulative : bool
          Compute cumulative, n(>Muv), or non-cumulative,
          dn/dMuv, UV luminosity function.
          fiducial: False

        add_dust : bool
          Add dust attenuation to UV.
          fiducial: False

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        completeness_correction : bool
          Apply completeness correction.
          fiducial: False

        completeness_correction_type : str
          Type of completeness correction:
            'parameterized' / 'numerical'
          fiducial: 'parameterized'

        Returns
        -------
        lum_bins : np.array
          Halpha luminosity (in erg/s).

        phi : np.array
          Number density of objects, units depend if
          cumulative or not.

        '''
        lum_list = np.log10(self.get_EmL_lum(emission_line=emission_line, SP_param_nr=SP_param_nr, exclude_contam_halos=exclude_contam_halos))
        idx_good = np.isfinite(lum_list)
        lum_list = lum_list[idx_good]
        lum_bins = np.arange(np.max([35.0, np.min(lum_list)]), np.min([45.0, np.max(lum_list)]), bin_size)
        lum_bins_center = lum_bins[:-1] + 0.5*np.diff(lum_bins)
        if completeness_correction:
            if (completeness_correction_type == 'numerical'):
                weights = 10**completeness.get_completeness_correction_numerical(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx_good]), self.redshift)
            elif (completeness_correction_type == 'parameterized'):
                weights = 10**completeness.get_completeness_correction_parametrized(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx_good]), self.redshift)
        else:
            weights = None
        hist, bin_edges = np.histogram(lum_list, bins=lum_bins, weights=weights)
        if add_dust:
            lum_bins_center_d = dust_attenuation.add_dust_attenuation_Ha(lum_bins_center, self.redshift)
        else:
            lum_bins_center_d = lum_bins_center
        if cumulative:
            LFcum = np.cumsum(hist)/volume_box
            return(lum_bins_center_d, LFcum)
        else:
            LF = hist/(np.diff(lum_bins)*volume_box)
            return(lum_bins_center_d, LF)

    def compute_cSFRD(self, SFR_limit_in=0.3, time_interval=200.0, volume_box=100.0**3, exclude_contam_halos=True, completeness_correction=False, completeness_correction_type='parameterized'):
        '''
        Computes cosmic star-formation rate density.

        Parameters
        ----------
        SFR_limit_in : float
          SFR limit in Msun/yr to which galaxies contribute
          to the cosmic star-formation rate density.
          fiducial: 0.3

        time_interval : float
          Time in which SFR is comupted, in Myr.
          fiducial: 200

        volume_box : float
          Volume of the box in Mpc^3.
          fiducial: 100.0^3

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        completeness_correction : bool
          Apply completeness correction.
          fiducial: False

        completeness_correction_type : str
          Type of completeness correction:
            'parameterized' / 'numerical'
          fiducial: 'parameterized'

        Returns
        -------
        total_SFRD : float
          cSFRD in log units Msun/yr/Mpc^3.

        '''
        SFR_list = self.get_SFR(time_interval, exclude_contam_halos=exclude_contam_halos)
        idx = (SFR_list >= SFR_limit_in) & np.isfinite(SFR_list)
        if completeness_correction:
            if (completeness_correction_type == 'numerical'):
                weights = 10**completeness.get_completeness_correction_numerical(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx]), self.redshift)
            elif (completeness_correction_type == 'parameterized'):
                weights = 10**completeness.get_completeness_correction_parametrized(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx]), self.redshift)
        else:
            weights = 1.0
        total_SFR = np.sum(weights*SFR_list[idx])
        total_SFRD = total_SFR/volume_box
        return(total_SFRD)

    def compute_cSMD(self, integral_of_SFR=True, SP_param_nr='4', Mseed=0.0, R=0.0, M_limit_in=10**8, volume_box=100.0**3, exclude_contam_halos=True, completeness_correction=False, completeness_correction_type='parameterized'):
        '''
        Computes cosmic stellar mass density.

        Parameters
        ----------
        integral_of_SFR : bool
          Specifies whether we define the stellar mass to be
          the integral of the past SFR or to be the actual
          mass in stars and remnants (i.e. subtracting return).

        Mseed : float
          The seed mass at BB.
          fiducial: 0

        R : float
          Fraction of the mass that is converted into stars,
          as measured by the SFR, is promptly (we will assume
          instantaneously) returned to the interstellar medium.
          The remaining fraction (1-R) stays in form of
          long-lived stars.
          see: http://adsabs.harvard.edu/abs/2016MNRAS.455.4183V
          fiducial: 0.29 (Salpeter)

        M_limit_in : float
          Mass limit in Msun to which galaxies contribute
          to the cosmic stellar mass density.
          fiducial: 10^8

        volume_box : float
          Volume of the box in Mpc^3.
          fiducial: 100.0^3

        exclude_contam_halos : bool
          Exclude contaminated halos.
          fiducial: True

        completeness_correction : bool
          Apply completeness correction.
          fiducial: False

        completeness_correction_type : str
          Type of completeness correction:
            'parameterized' / 'numerical'
          fiducial: 'parameterized'

        Returns
        -------
        total_SMD : float
          SMD in log units Msun/Mpc^3.

        '''
        M_list = self.get_stellar_mass(integral_of_SFR=integral_of_SFR, SP_param_nr=SP_param_nr, Mseed=Mseed, R=R, exclude_contam_halos=exclude_contam_halos)
        idx = (M_list >= M_limit_in) & np.isfinite(M_list)
        if completeness_correction:
            if (completeness_correction_type == 'numerical'):
                weights = 10**completeness.get_completeness_correction_numerical(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx]), self.redshift)
            elif (completeness_correction_type == 'parameterized'):
                weights = 10**completeness.get_completeness_correction_parametrized(np.log10(self.get_halo_mass(exclude_contam_halos=exclude_contam_halos)[idx]), self.redshift)
        else:
            weights = 1.0
        total_M = np.sum(weights*M_list[idx])
        total_SMD = total_M/volume_box
        return(total_SMD)



      