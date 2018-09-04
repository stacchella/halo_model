import numpy as np
from scipy.optimize import fsolve


def remove_dust_attenuation(magUV_in, dust_type, redshift_in, with_scatter=False):
    '''
    This function removes the dust attenuation / performs
    the dust attenuation correction by computing the
    dust attenuation at 1600 A (A1600) from the
    Bouwens et al 2015 relation.

    dust_type (str):
        Meurer+99
        Calzetti+00
        Reddy+15
        Gordon+03 (SMC)
    '''
    # ensure inputs are arrays
    if type(magUV_in) is not np.ndarray:
        magUV_in = np.array([magUV_in])
    if type(redshift_in) is not np.ndarray:
        redshift_in = np.array([redshift_in])
    # # values for beta-MUV relation from Bouwens et al 2015 + extrapolation to z=10
    B15_z_list = np.array([2.5, 3.8, 5.0, 5.9, 7.0, 8.0, 10.0, 14.0])
    B15_beta_list = np.array([-1.70, -1.85, -1.91, -2.00, -2.05, -2.13, -2.20, -2.20])
    B15_dbdM_list = np.array([-0.20, -0.11, -0.14, -0.20, -0.20, -0.20, -0.20, -0.20])
    # # interpolate beta-MUV relation
    beta_value = np.interp(redshift_in, B15_z_list, B15_beta_list, left=B15_beta_list[0], right=np.nan)
    dbdM_value = np.interp(redshift_in, B15_z_list, B15_dbdM_list, left=B15_dbdM_list[0], right=np.nan)
    #beta_value = -0.09*redshift_in-1.49
    #dbdM_value = -0.007*redshift_in-0.09
    # get average beta
    avg_b = dbdM_value*(magUV_in+19.5) + beta_value
    # convert beta to A1600
    # see Table 3 of Reddy+18
    if dust_type not in ['Meurer+99', 'Calzetti+00', 'Reddy+15', 'Gordon+03']:
        print 'NOT KNOWN DUST TYPE (IRX-beta)!'
    if (dust_type == 'Meurer+99'):
        a_dust = 1.99
        b_dust = 4.43
    elif (dust_type == 'Calzetti+00'):
        a_dust = 2.13
        b_dust = 5.57
    elif (dust_type == 'Reddy+15'):
        a_dust = 1.82
        b_dust = 4.47
    elif (dust_type == 'Gordon+03'):
        a_dust = 1.07
        b_dust = 2.79
    if with_scatter:
        A1600 = b_dust + a_dust*(avg_b + np.random.normal(loc=0.0, scale=0.34, size=len(avg_b)))
    else:
        A1600 = b_dust + 0.2*np.log(10)*a_dust**2*0.34**2+a_dust*avg_b
    A1600[np.isnan(A1600)] = np.zeros(np.sum(np.isnan(A1600)))
    A1600[A1600 < 0.0] = np.zeros(np.sum(A1600 < 0.0))
    magUV_without_dust = magUV_in - A1600
    return(magUV_without_dust)


# def remove_dust_attenuation(magUV_in, redshift_in, with_scatter=False):
#     '''
#     This function removes the dust attenuation / performs
#     the dust attenuation correction by computing the
#     dust attenuation at 1600 A (A1600) from the
#     Bouwens et al 2015 relation.
#     '''
#     # ensure inputs are arrays
#     if type(magUV_in) is not np.ndarray:
#         magUV_in = np.array([magUV_in])
#     if type(redshift_in) is not np.ndarray:
#         redshift_in = np.array([redshift_in])
#     # values for beta-MUV relation from Bouwens et al 2015 + extrapolation to z=10
#     Mag_corr = np.array([-22.979, -21.986, -20.993, -20.028, -18.994, -18.000, -16.280])
#     AV_corr = 1.0/(redshift_in-4.0+1.0)**(1./4.)*np.array([2.5, 1.7, 1.3, 0.9, 0.7, 0.4, 0.25])
#     A1600 = np.interp(magUV_in, Mag_corr, AV_corr, left=AV_corr[0], right=0.1)
#     A1600[np.isnan(A1600)] = np.zeros(np.sum(np.isnan(A1600)))
#     A1600[A1600 < 0.0] = np.zeros(np.sum(A1600 < 0.0))
#     magUV_without_dust = magUV_in - A1600
#     return(magUV_without_dust)


def fct_solve(mag_with_dust, mag_without_dust, dust_type, redshift):
    '''
    Define function to minimize.
    '''
    return(mag_without_dust-remove_dust_attenuation(mag_with_dust, dust_type, redshift))


def add_dust_attenuation(magUV_in, dust_type, redshift_in, with_scatter=False):
    '''
    Adds dust according to the Bouwens et al 2015 relation.
    '''
    mag_list = np.linspace(-30.0, -10.0)
    mag_list_with_dust = fsolve(fct_solve, mag_list, args=(mag_list, dust_type, redshift_in))
    mag_with_dust = np.interp(magUV_in, mag_list, mag_list_with_dust)
    return(mag_with_dust)


def calzetti(wave, tau_v=1, R_v=4.05, **kwargs):
    """Calzetti et al. 2000 starburst attenuation curve, with
    extrapolations to the FUV and NIR.

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the attenuation curve.

    :param R_v: (default: 4.05)
        The ratio of total selective extinction, parameterizing the slope of
        the attenuation curve.  A_v = R_v * E(B-V)

    :returns tau:
        The optical depth at each wavelength.
    """
    # optical/NIR
    k1 = lambda x: 2.659 * (-1.857 + 1.040 * x)
    # UV
    k2 = lambda x: 2.659 * (-2.156 + 1.509 * x - 0.198 * x**2. + 0.011 * x**3.)

    # get slopes at edges and k(5500)
    uv = np.array([0.12, 0.13]) * 1e4
    kuv = k2(1e4 / uv) + R_v
    uv_slope = np.diff(kuv) / np.diff(uv)
    ir = np.array([2.19, 2.20]) * 1e4
    kir = k1(1e4 / ir) + R_v
    ir_slope = np.diff(kir) / np.diff(ir)
    k_v = k2(1e4 / 5500.) + R_v

    # define segments
    uinds = (wave >= 1200.) & (wave < 6300)  # uv
    oinds = (wave >= 6300.) & (wave <= 22000)  # optical
    xinds = (wave < 1200.)  # xuv
    iinds = (wave > 22000.)  # ir

    # do it
    x = 1e4 / wave
    ktot = oinds * (k1(x) + R_v)
    ktot += uinds * (k2(x) + R_v)
    ktot += xinds * (kuv[0] + (wave - uv[0]) * uv_slope)
    ktot += iinds * (kir[1] + (wave - ir[1]) * ir_slope)

    ktot[ktot < 0] = 0
    tau_lambda = tau_v * (ktot / k_v)
    return tau_lambda


