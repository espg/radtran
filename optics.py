# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 17:04:03 2012

@author: - Shane Grigsby <shane@geog.ucsb.edu>
"""
from pylab import *
from scipy.constants import *
from pytz import timezone
from pytz import utc
from datetime import datetime

import time

def excel_fix(array):
    for i in array:
        print i

ang = 1E-10
nm = 1E-9
um = 1E-6
Cm = 1E-2
hH = 1.0
kH = 1E3
mH = 1E6
gH = 1E9
tH = 1E12
wn =0

def m(length, unit=1):
    return length*unit
    
def freq(wavelength, unit=1):
    frequency =(c / (wavelength*unit))
    return frequency

def wave(frequency, unit=1):
    if unit == wn:
        return (1.0/(frequency*100))
    else:
        wavelength = (c / (frequency*unit))
        return wavelength

def waven(wavelength, unit=1):
    if unit == (hH or kH or mH or gH or tH):
        wavenumber = wave(wavelength, unit)
    wavenumber = (1.0 / ((wavelength*unit)))/100.0
    return wavenumber
    
def Q(wavelength, unit=1):
    w = m(wavelength, unit)
    q = ((h*c)/w)
    return q
    

def Xwave(wavelength, temp, unit=1):
    X_wave = (h*c)/(k*(wavelength*unit)*temp)
    return X_wave
    
def Lwave(wavelength, temp, unit=1):
    """Calculates L given wavelength and Temp
    To get M, multiply by pi
    Note that units are:
        W * m**-2 * sr**-1 * m**-1
    I.e, units are given in meter of spectrum
    multiply by nm to get:
        W * m**-2 *sr**-1 nm**-1"""
    X_funct= Xwave(wavelength, temp, unit)
    L=(2*h*(c**2))/(((wavelength*unit)**5)*(exp(X_funct)-1))
    return L

## WTF is this ???    
def kineticT(wavelength, brightnessTemp, emissivity):
    """Converts from a brightness temp to kinetic temp given a known emissivity"""
    TempK = (h*c)/(k*wavelength*log(1 + emissivity*exp((h*c)/(k*wavelength*brightnessTemp)) - emissivity)) 
    return TempK

def kinetic_bright_T(wavelength, L):
    """Defines brightness temperature, (not adjusted for emissivity); use kineticT to convert to kinetic temperature
    Units are in meters, L is expected in Watts per meter^2 per meter (of spectrum) per steridan"""
    y = (2*h*(c**2))/(L*((wavelength)**5))
    TempK = (h*c)/(k*wavelength*log(1 + y)) 
    return TempK

def brightnessT(wavelength, kTemp, emissivity):
    """Returns the brightness temp that would be observed for a known kinetic temperature and emissivity""" 
    TempB = ((h*c)/(k*wavelength*log((exp((h*c)/(k*wavelength*kTemp)) + 
                    emissivity - 1.0)/emissivity))) 
    return TempB

def PlanetTemp(albedo,emissivity,LongM,ShortM):
    '''Caluclates temperature in kelvin.
    Expects top of atmosphere M values in W*m**-2
    Accounts for spherical shape of planet
    
    Note: Albedo and emissivity both expect decimal values''' 
    incoming = .25*(ShortM*(1.0-albedo)+LongM*emissivity)
    tempPlanet = (incoming/(emissivity*sigma))**.25
    return tempPlanet

def kelvin(celcius):
    return (celcius + 273.15)

def a2k(absorb_coeff,wave_length,alpha_unit=Cm):
    k_value = ((absorb_coeff*(alpha_unit**-1.0)*wave_length) / (4.0*pi))
    return k_value

def k2a(extinct_coeff, wavelength, output_unit_alpha=Cm):
    a_value = (4.0*pi*(extinct_coeff / wavelength)*output_unit_alpha)
    return a_value

def Rll(angleIncidence,n2,n1=1,k2=0,k1=0):
    thetaTrans = refractAngle(angleIncidence,n2,n1)
    top = ((n2*cos(angleIncidence) - n1*cos(thetaTrans))**2 +
          (k2*cos(angleIncidence)  - k1*cos(thetaTrans))**2)
    bot = ((n2*cos(angleIncidence) + n1*cos(thetaTrans))**2 +
          (k2*cos(angleIncidence)  + k1*cos(thetaTrans))**2)
    return (top/bot)
    
def Rppd(angleIncidence,n2,n1=1,k2=0,k1=0):
    thetaTrans = refractAngle(angleIncidence,n2,n1)
    top = ((n2*cos(thetaTrans) - n1*cos(angleIncidence))**2 +
          (k2*cos(thetaTrans)  - k1*cos(angleIncidence))**2)
    bot = ((n2*cos(thetaTrans) + n1*cos(angleIncidence))**2 +
          (k2*cos(thetaTrans)  + k1*cos(angleIncidence))**2)
    return (top/bot)
    
    
def refractAngle(angleIncidence,n2,n1=1):
    """Simple application of Snell's law.
    Given two refractive indexes and an incoing angle of
    incidence, will calculate transmission angle"""
    thetaR = arcsin((n1*sin(angleIncidence))/n2)
    return thetaR
    
###Time Equations###

##def ET(jn):
    
    
    
def derivative(f):
    """Needed for the newton function
    Don't be lazy----> errors propagate when taking numerical derivatives
    Don't use this for derivatives greater then the 2nd order"""

    def df(x, h=0.1e-5):
        return (f(x+h) - f(x))/h
    return df

def newton(g, x_inital, n):
    """Functions should be of one variable, and defined similarly
    to the following example:

        def g(x): return 1075 +2*x - 6*x**2 + 0.4*x**3
        
    Note that '2x' is not equivlent to '2*x'
    Also note that the function is 'unstable' with very, truely, bad guesses
    Also note that taking beyond the 2nd derivative will mess with the floating 
    point registers --> do the calculus symbolically to avoid this;
    i.e., only use this function to produce zeros!
    """
    dg = derivative(g)
    for i in range(n):
        g_prime = (((x_inital*dg(x_inital)) - g(x_inital)) / (dg(x_inital)))
        print(g_prime)
        x_inital = g_prime
    return g_prime
    
    
def bowen_ratio(top_temp, bot_temp,VP_top,VP_bot, surf_P=101325.,Cp=1005.,L_vap=2.468E6):
    bowen = ((surf_P*Cp*(top_temp-bot_temp)) / (0.622*L_vap*(VP_top - VP_bot)))
    return bowen
###Julian Date and UTC Stuff###
##Horribly ugly###

def o_thick(E_toa, E_local, zenith,z=0,H=8434.0):
    return ((-log((E_local/(E_toa*cos(zenith))))/(m_air(zenith)*exp((-z/H)))))

def m_air(zenith):
    return (1.0/(cos(zenith)+0.15*(((90-degrees(zenith)+3.885)**-1.253))))

def E_z(E_toa,zenith,o_thickness, z=0,H=8434.0):
    return ((E_toa*cos(zenith))*(exp((-o_thickness*m_air(zenith)*exp(-z/H)))))    

#String formatting
fmt = '%Y-%m-%d %H:%M:%S %Z%z'

def eq_time(jd):
    return (0.17*(sin(4*pi*(jd-80)/373))-(0.129*sin(2*pi*(jd-8)/355)))    
    
def solar_noon(lon,year,month,day,time_zone):
    utc_noon = 12.0 - eq_time(JD(year,month,day))
    local_noon_utc = utc_noon + (-lon/15.0)  ##local time of solarnoon in utc
    h_m_s = decimal_to_base60(local_noon_utc) ##h_m_s is hours_minutes_seconds
    utc_dt = datetime(year,month,day,h_m_s[1],h_m_s[2],int(round(h_m_s[3])),tzinfo=utc)
    local_noon_local_time = utc_dt.astimezone(timezone(time_zone))
    string_local_time = local_noon_local_time.strftime(fmt)
    print string_local_time    
    noon_time = (str(local_noon_local_time.hour) + ' ' + str(local_noon_local_time.minute) 
        + ' ' + str(local_noon_local_time.second))
    offset = float(string_local_time[23] + string_local_time[24] + string_local_time[25])
    return base60_to_decimal(noon_time), utc_noon, offset

def local_time(time_zone,year,month,day,hour,minute,second=0):
    utc_dt = datetime(year,month,day,hour,minute,second,tzinfo=utc)
    local_time_dt = utc_dt.astimezone(timezone(time_zone))
    print(local_time_dt)
    return local_time_dt
    
def h_theta(lon, solar_noon_utc, offset, hours=12, minutes=0, seconds=0):
    """Note that the hour angle is give in degrees
    Offset can be taken from solar_noon function; note the sign
    NOTE WELL: hours are in military / 24 hour clock
    
    This function will also calculate w, the solar azimuth"""
    time1 = (str(hours) + ' ' + str(minutes) 
        + ' ' + str(seconds))
    local_decimal_time = base60_to_decimal(time1)
    utc_time = local_decimal_time - offset
    w = (solar_noon_utc - utc_time)*15.0
    hour_angle = w - lon
    return hour_angle, w

def day_stats(lat,lon,year,month,day,time_zone,refraction=True):
    noon = solar_noon(lon,year,month,day,time_zone)
    dec = solar_declination(year,month,day)
    if refraction == False:
        h_angle = (arccos(-1*tan(dec)*tan(radians(lat))))/radians(15.0)
    else:
        h_angle = (arccos(cos(radians(90.8333))-tan(dec)*tan(radians(lat))) / 
                  radians(15.0))
    sunrise = noon[0] - h_angle
    sunset = noon[0] + h_angle
    daylength = sunset - sunrise
#    print decimal_to_base60(sunrise)
#    print decimal_to_base60(sunset)
#    print decimal_to_base60(daylength)
    return sunrise,sunset,daylength

def sun_z_and_a(lat,lon,year,month,day,hours,minutes,seconds,time_zone,exposure=0, slope_angle=0):
    """Calculates local solar zenith and azimuth given a slope, or normal zenith without"""
    """North is +180 degrees, South is 0 degrees; degrees decrease clockwise"""
    s_noon = solar_noon(lon,year,month,day,time_zone)
    s_dec = solar_declination(year,month,day)
    h_angle = h_theta(lon, s_noon[1], s_noon[2], hours, minutes, seconds)
    sun_z = arccos((sin(s_dec)*sin(radians(lat)) + 
            cos(s_dec)*cos(radians(lat))*cos(radians(h_angle[0]))))
    sun_a = arctan2((cos(s_dec)*sin(radians(h_angle[0]))), 
            ((cos(s_dec)*sin(radians(lat))*cos(radians(h_angle[0]))) - 
            sin(s_dec)*cos(radians(lat))))
    sun_z_local = arccos((cos(sun_z)*cos(radians(slope_angle)) + 
                  sin(sun_z)*sin(radians(slope_angle))*cos(sun_a-radians(exposure))))
    sun_a_local = arctan2((sin(sun_z)*sin(sun_a-radians(exposure))), 
                  (sin(sun_z)*cos(radians(slope_angle))*cos(sun_a - exposure) -
                  cos(sun_z)*sin(radians(slope_angle))))
 #   print h_angle
 #   print degrees(s_dec)
    print degrees(sun_a)     
    print("Azimuth (degrees): " + str(degrees(sun_a)))
    print degrees(sun_z)
    print("Zenith (degrees): " + str(degrees(sun_z)))
    print("Local Azimuth (degrees): " + str(degrees(sun_a_local)))
    print("Local Zenith (degrees): " + str(degrees(sun_z_local)))
    return sun_z  
                 
def solar_declination(year,month,day):
    """Note, returns radians"""
    mod_jd =  2*pi*JD(year,month,day)/365.0
    declination = (0.006918-0.399912*cos(mod_jd)+0.070257*sin(mod_jd) - 
                   0.006758*cos(2*mod_jd) +0.000907*sin(2*mod_jd) - 
                   0.002697*cos(3*mod_jd)+0.00148*sin(3*mod_jd))
    return declination

def GAMMA(year,month,day):
    result_g = (1.00011+ 0.034221 * cos((JD(year,month,day)*2*pi)/365.0) + 
             0.00128 * sin((JD(year,month,day)*2*pi)/365.0)+ 0.000719 * 
             cos(2* (JD(year,month,day)*2*pi)/365.0)+ 0.000077 * 
             sin(2* (JD(year,month,day)*2*pi)/365.0))
    return result_g

#def sun_a(s_noon,time,lat,j_day):

def JD(year,month,day):
    t = time.mktime((year,month,day,0,0,0,0,0,0))
    return float(time.gmtime(t)[7])
    
    
###Thanks to the authors of the below code for decimal degree conversions

def base60_to_decimal(xyz,delimiter=None):
    """Decimal value from numbers in sexagesimal system. 

    The input value can be either a floating point number or a string
    such as "hh mm ss.ss" or "dd mm ss.ss". Delimiters other than " "
    can be specified using the keyword ``delimiter``.
    """
    divisors = [1,60.0,3600.0]
    xyzlist = str(xyz).split(delimiter)
    sign = -1 if xyzlist[0].find("-") != -1 else 1
    xyzlist = [abs(float(x)) for x in xyzlist]
    decimal_value = 0 

    for i,j in zip(xyzlist,divisors): # if xyzlist has <3 values then
                                      # divisors gets clipped.
        decimal_value += i/j

    decimal_value = -decimal_value if sign == -1 else decimal_value
    return decimal_value
    
def decimal_to_base60(deci,precision=1e-8):
    """Converts decimal number into sexagesimal number parts. 

    ``deci`` is the decimal number to be converted. ``precision`` is how
    close the multiple of 60 and 3600, for example minutes and seconds,
    are to 60.0 before they are rounded to the higher quantity, for
    example hours and minutes.
    """
    sign = "+" # simple putting sign back at end gives errors for small
               # deg. This is because -00 is 00 and hence ``format``,
               # that constructs the delimited string will not add '-'
               # sign. So, carry it as a character.
    if deci < 0:
        deci = abs(deci)
        sign = "-" 

    frac1, num = math.modf(deci)
    num = int(num) # hours/degrees is integer valued but type is float
    frac2, frac1 = math.modf(frac1*60.0)
    frac1 = int(frac1) # minutes is integer valued but type is float
    frac2 *= 60.0 # number of seconds between 0 and 60 

    # Keep seconds and minutes in [0 - 60.0000)
    if abs(frac2 - 60.0) < precision:
        frac2 = 0.0
        frac1 += 1
    if abs(frac1 - 60.0) < precision:
        frac1 = 0.0
        num += 1 

    return (sign,num,frac1,frac2)
    
#############
#############
#Other used functions that I found with the decimal code


import math
MJD0 = 2400000.5 # 1858 November 17, 00:00:00 hours 



def julian_date(year,month,day,hour,minute,second):
    """Given year, month, day, hour, minute and second return JD.

    ``year``, ``month``, ``day``, ``hour`` and ``minute`` are integers,
    truncates fractional part; ``second`` is a floating point number.
    For BC year: use -(year-1). Example: 1 BC = 0, 1000 BC = -999.
    """
    MJD0 = 2400000.5 # 1858 November 17, 00:00:00 hours 

    year, month, day, hour, minute =\
    int(year),int(month),int(day),int(hour),int(minute)

    if month <= 2:
        month +=12
        year -= 1 

    modf = math.modf
    # Julian calendar on or before 1582 October 4 and Gregorian calendar
    # afterwards.
    if ((10000L*year+100L*month+day) <= 15821004L):
        b = -2 + int(modf((year+4716)/4)[1]) - 1179
    else:
        b = int(modf(year/400)[1])-int(modf(year/100)[1])+\
            int(modf(year/4)[1]) 

    mjdmidnight = 365L*year - 679004L + b + int(30.6001*(month+1)) + day

    fracofday = base60_to_decimal(\
      " ".join([str(hour),str(minute),str(second)])) / 24.0 

    return MJD0 + mjdmidnight + fracofday

def caldate(mjd):
    """Given mjd return calendar date. 

    Retrns a tuple (year,month,day,hour,minute,second). The last is a
    floating point number and others are integers. The precision in
    seconds is about 1e-4. 

    To convert jd to mjd use jd - 2400000.5. In this module 2400000.5 is
    stored in MJD0.
    """
    MJD0 = 2400000.5 # 1858 November 17, 00:00:00 hours 

    modf = math.modf
    a = long(mjd+MJD0+0.5)
    # Julian calendar on or before 1582 October 4 and Gregorian calendar
    # afterwards.
    if a < 2299161:
        b = 0
        c = a + 1524
    else:
        b = long((a-1867216.25)/36524.25)
        c = a+ b - long(modf(b/4)[1]) + 1525 

    d = long((c-122.1)/365.25)
    e = 365*d + long(modf(d/4)[1])
    f = long((c-e)/30.6001)

    day = c - e - int(30.6001*f)
    month = f - 1 - 12*int(modf(f/14)[1])
    year = d - 4715 - int(modf((7+month)/10)[1])
    fracofday = mjd - math.floor(mjd)
    hours = fracofday * 24.0 

    sign,hour,minute,second = decimal_to_base60(hours)

    return (year,month,day,int(sign+str(hour)),minute,second)

if __name__ == '__main__':
    print "Julian date for 2010/1/1 13:20:12.3456 : ",
    j = julian_date(2010,1,1,13,20,12.3456)
    print j
    print "Calendar date for MJD "+ str(j-MJD0) + " (jd = " + str(j)+" )"
    print "Year: {0}, Month: {1}, Day: {2}, Hour: {3}, Minute: {4},\
    Second: {5:8.5f}".format(*caldate(j-MJD0))
