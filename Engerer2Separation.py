# The Engerer Separation Model was first described and validated in the
# paper: Engerer, N. A. 2015. Minute resolution estimates of the diffuse
# fraction of global irradiance for southeastern Australia. Solar Energy.
# 116, 215-237. Section 4. Development of new models.
# The paper presents 3 variations of the Engerer separation model: 1, 2 & 3
# The Engerer1 is fit to Non-Cloud Enhancement data, the Engerer2 to Cloud
# Enhancement data and Engerer3 to clear-sky data only.

# The predictor variables to fit B0-5 are based on:
# kt = clearness index
# zen = zenith angle
# ast = apparent solar time
# dktc = ktc - kt (where ktc is the clear-sky clearness index = ghi,cs/e_exth.
# k_de = Proportion of kd attributable to cloud enhancement

# The original formulation was found to be not general enough to other parts of the world and at different temporal
# resolutions, and so in a new paper, the parameters have been updated.
# Bright, Jamie M. & Engerer, Nicholas A. 2019. Engerer2: Global re-parameterisation, update and validation of an
# irradiance separation model at different temporal resolutions. Journal of Renewable & Sustainable Energy. 11(2), xxx.

# The new parameters are defined in the paper as:
# ===================================================================================================================
# Parameter, Predictor,   Orig 1-min, New 1-min,  5-min       10-min      15-min      30-min,     1-hr,       1-day
# C,          - &         0.042336,   0.10562,    0.093936,   0.079965,   0.065972,   0.032675,   -0.0097539, 0.32726
# B0,         - &         -3.7912,    -4.1332,    -4.5771,    -4.8539,    -4.7211,    -4.8681,    -5.3169,    -9.4391
# B1,         K_t,         7.5479,    8.2578,     8.4641,     8.4764,     8.3294,     8.1867,     8.5084,     17.113
# B2,         AST,        -0.010036,  0.010087,   0.010012,   0.018849,   0.0095444,  0.015829,   0.013241,   0.13752
# B3,         zen,        0.003148,   0.00088801, 0.003975,   0.0051497,  0.0053493,  0.0059922,  0.0074356,  -0.024099
# B4,         dKtc,       -5.3146,    -4.9302,    -4.3921,    -4.1457,    -4.169,     -4.0304,    -3.0329,    6.6257
# B5,         Kde,        1.7073,     0.44378,    0.39331,    0.37466,    0.39526,    0.47371,    0.56403,    0.31419
# ===================================================================================================================

# INPUT PARAMETERS
# ghi               = global horizontal irradiance [W/m^2]. Numpy array.
# time              = in the format [yyyy, mm, dd, HH, MM, SS] [UTC]. Numpy array
# latitudes         = either a latitude per measurement or one for all. Numpy array.
# longitude         = either a longitude per measurement or one for all. Numpy array.
# averaging_period  = x in minutes. Must be 1, 5, 10, 15, 30, 60 or 1440. Numeric.

# OUTPUT PARAMETERS
# kd    = diffuse fraction = Edh/ghi [fraction]. Numpy array.
# dif   = diffuse horizontal irradiance [W/m^2]. Numpy array.
# dni   = direct/beam normal irradiance [W/m^2]. Numpy array.

# EXAMPLE USAGE
# time = np.array([[2018, 5, 1, 8, 1, 0], [2016, 3, 5, 10, 15, 0], [2015, 9, 10, 13, 17, 0], [2014,1,20,20,18,0]])
# ghi = np.array([[300], [400], [1000], [600]])
# lat = np.array([[-10], [10], [40], [-40]])
# lon = np.array([[-60], [-5], [20], [170]])
# averaging_period = 1
# kd, dif, dni = engerer2separation(ghi, lat, lon, time, averaging_period)
# print(kd, dif, dni)

import datetime
import numpy as np

# Expected performance of this model will compare values to nans, e.g. nan < 1. The answer we want from this is still
# nan, however, numpy will throw a warning, hence, we ignore it.
np.warnings.filterwarnings('ignore')

# generalised logistic function formulaton of the Engerer2 separation model.
def engerer2(cc, bb0, bb1, bb2, bb3, bb4, bb5, zzen, aast, ddktc, kk_de, kkt):
    diff_frac = cc + (1 - cc) / (1 + np.exp(bb0 + bb1 * kkt + bb2 * aast + bb3 * zzen / (np.pi / 180) + bb4 * ddktc)) + bb5 * kk_de
    return diff_frac

def engerer2separation(ghi, latitudes, longitudes, time, averaging_period):

    # Extract the time data
    year = time[:, 0]
    month = time[:, 1]
    day = time[:, 2]
    hour = time[:, 3]
    minute = time[:, 4]
    second = time[:, 5]

    # decimal time (e.g. 10:30am = 10.5, and 10:30pm = 22.5).
    dt = (hour * 60 + minute + second / 3600) / 60
    # NOTE: this code produces a 1-x-4 array, which causes errors later. The usual .transpose() or T do not seem to work
    # with dt, as such decimal_time is redefined within the loop following...

    # Day of year
    # NOTE: I could not figure a method of vectorising this as datetime.datetime() requires a single int to produce the
    # unix timestamps. I would recommend passing in time as a unix POSIX time so that day of year is simply derived
    # without this loop. If you do figure out how to remove this loop, i would appreciatae if you pushed the changes to
    # the github reopository detailed in the paper.
    day_of_year = np.zeros((len(ghi), 1))
    decimal_time = np.zeros((len(ghi), 1))

    for i in range(len(ghi)):
        t = datetime.datetime(year=year[i], month=month[i], day=day[i], hour=hour[i], minute=minute[i], second=second[i], tzinfo=datetime.timezone.utc).timestamp()
        t0 = datetime.datetime(year=year[i], month=1, day=1, tzinfo=datetime.timezone.utc).timestamp()
        day_of_year[i, 0] = np.ceil((t - t0)/(24*60*60))
        decimal_time[i] = dt[i]

    # Equation of time
    beta_eot = (360 / 365.242) * (day_of_year - 1)
    eot = (0.258 * np.cos(np.pi / 180 * beta_eot) - 7.416 * np.sin(np.pi / 180 * beta_eot) - 3.648 * np.cos(
        np.pi / 180 * 2 * beta_eot) - 9.228 * np.sin(np.pi / 180 * 2 * beta_eot))

    # Local solar noon
    lsn = 12 - longitudes / 15 - 1 * (eot / 60)

    # Hour angle
    hour_angle = (decimal_time - lsn) * 15
    hour_angle[hour_angle >= 180] = hour_angle[hour_angle >= 180] - 360
    hour_angle[hour_angle <= -180] = 360 + hour_angle[hour_angle <= -180]

    # Declination angle
    phi_r = (2 * np.pi / 365) * (day_of_year + (decimal_time / 24) - 1)
    delta_r = 0.006918 - 0.399912 * np.cos(phi_r) + 0.070257 * np.sin(phi_r) - 0.006758 * np.cos(2 * phi_r) \
        + 0.000907 * np.sin(2 * phi_r) - 0.002697 * np.cos(3 * phi_r) + 0.001480 * np.sin(3 * phi_r)

    # Zenith angle
    zen = np.arccos(np.sin(np.pi / 180 * latitudes) * np.sin(delta_r) + np.cos(np.pi / 180 * latitudes) *
                    np.cos(delta_r) * np.cos(np.pi / 180 * hour_angle))

    # Extraterrestrial irradiance
    solar_constant = 1366.1
    beta = (2 * np.pi * day_of_year) / 365
    e_extn = solar_constant * (1.00011 + 0.034221 * np.cos(beta) + 0.00128 * np.sin(beta) + 0.000719 * np.cos(2 * beta)
                               + 0.000077 * np.sin(2 * beta))

    # Extraterrestrial horizontal irradiance, with night time correction.
    e_exth = e_extn * np.cos(zen)
    e_exth[zen > np.pi / 180 * 90] = 0

    # The TJ clear-sky model
    a_tj = 1160 + 75 * np.sin((360 * (day_of_year - 275)) / 365)
    k_tj = 0.174 + 0.035 * np.sin((360 * (day_of_year - 100)) / 365)
    c_tj = 0.095 + 0.04 * np.sin((360 * (day_of_year - 100)) / 365)
    dnics = a_tj * np.exp(-k_tj / np.cos(zen))
    dnics[zen > np.pi / 180 * 90] = 0
    ghics = dnics * np.cos(zen) + c_tj * dnics
    ghics[zen > np.pi / 180 * 90] = 0

    # Calculate predictors.
    # clearness index
    kt = ghi / e_exth

    # clearness index of clear sky
    ktc = ghics / e_exth

    # deviation of clearness from clear-sky-clearness
    dktc = ktc - kt

    # apparent solar time
    ast = hour_angle / 15 + 12
    # quality check asserted in the original Engerer2 R code.
    ast[ast < 0] = abs(ast[ast < 0])

    # Proportion of Kd attributable to cloud enhancement from Engerer2 R code.
    cloud_enhancement_estimate = ghi - ghics
    cloud_enhancement_estimate[cloud_enhancement_estimate < 0.015] = 0
    k_de = cloud_enhancement_estimate / ghi

    # QC the predictors
    # it is likely/possible for there to be inf values and NA values.
    kt[kt == np.inf] = np.nan
    ktc[ktc == np.inf] = np.nan
    dktc[dktc == np.inf] = np.nan
    ast[ast == np.inf] = np.nan
    k_de[k_de == np.inf] = np.nan

    kt[zen > np.pi / 180 * 90] = np.nan
    ktc[zen > np.pi / 180 * 90] = np.nan
    dktc[zen > np.pi / 180 * 90] = np.nan
    ast[zen > np.pi / 180 * 90] = np.nan
    k_de[zen > np.pi / 180 * 90] = np.nan

    # Engerer version parametrisation

    if averaging_period == 1:
        parameters = (0.105620, -4.13320, 8.25780, 0.0100870, 0.000888010, -4.93020, 0.443780)
    elif averaging_period == 5:
        parameters = (0.0939360, -4.57710, 8.46410, 0.0100120, 0.00397500, -4.39210, 0.393310)
    elif averaging_period == 10:
        parameters = (0.0799650, -4.85390, 8.47640, 0.0188490, 0.00514970, -4.14570, 0.374660)
    elif averaging_period == 15:
        parameters = (0.0659720, -4.72110, 8.32940, 0.00954440, 0.00534930, -4.16900, 0.395260)
    elif averaging_period == 30:
        parameters = (0.0326750, -4.86810, 8.18670, 0.0158290, 0.00599220, -4.03040, 0.473710)
    elif averaging_period == 60:
        parameters = (-0.00975390, -5.31690, 8.50840, 0.0132410, 0.00743560, -3.03290, 0.564030)
    elif averaging_period == 1440:
        parameters = (0.327260, -9.43910, 17.1130, 0.137520, -0.0240990, 6.62570, 0.314190)
    else:
        raise Exception("not a valid averaging period, must be 1, 5, 10, 15, 30, 60, or 1440")

    c = parameters[0]
    b0 = parameters[1]
    b1 = parameters[2]
    b2 = parameters[3]
    b3 = parameters[4]
    b4 = parameters[5]
    b5 = parameters[6]

    # original Engerer 2 variables.
    # c = 4.2336E-2
    # b0 = -3.7912
    # b1 = 7.547948473
    # b2 = -1.0035892E-2
    # b3 = 3.147968E-3
    # b4 = -5.31460674
    # b5 = 1.707321551

    # Calculate the outputs
    # calculated the diffuse fraction
    kd = engerer2(c, b0, b1, b2, b3, b4, b5, zen, ast, dktc, k_de, kt)

    # quality control.
    kd[kd > 1] = 1
    kd[kd < 0] = 0

    # calculate the diffuse irradiance
    dif = ghi * kd

    # calculate the direct normal irradiance
    dni = (ghi - dif) / np.cos(zen)

    # QC
    dni[zen > np.pi / 180 * 90] = np.nan
    dif[zen > np.pi / 180 * 90] = np.nan

    return kd, dif, dni


time = np.array([[2018, 5, 1, 8, 1, 0], [2016, 3, 5, 10, 15, 0], [2015, 9, 10, 13, 17, 0], [2014, 1, 20, 20, 18, 0]])
ghi = np.array([[300], [400], [1000], [600]])
lat = np.array([[-10], [10], [40], [-40]])
lon = np.array([[-60], [-5], [20], [170]])
averaging_period = 1
kd, dif, dni = engerer2separation(ghi, lat, lon, time, averaging_period)
print("Diffuse fraction: ")
print(kd)
print("Diffuse Horizontal Irradiance: ")
print(dif)
print("Direct Normal Irradiance: ")
print(dni)
