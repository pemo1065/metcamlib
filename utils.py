from astropy import units as u
from astropy.coordinates import AltAz
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

def ra_dec_to_alt_az(ra, dec, lat, lon, elev, timestamp):
    observing_location = EarthLocation(lat=str(lat), lon=str(lon), height=elev*u.m)  
    observing_time = Time(timestamp, scale="utc")  
    aa = AltAz(location=observing_location, obstime=observing_time)

    skycoord = SkyCoord(ra * u.degree, dec * u.degree)
    altaz = skycoord.transform_to(aa)
    return altaz.alt.value, altaz.az.value
    
def alt_az_to_ra_dec(alt, az, lat, lon, elev, timestamp):
    observing_location = EarthLocation(lat=str(lat), lon=str(lon), height=elev*u.m)  
    observing_time = Time(timestamp, scale="utc")  
    altaz = SkyCoord(alt = alt*u.degree, az = az*u.deg, obstime = observing_time, frame = 'altaz', location = observing_location)
    icrs = altaz.icrs
    return icrs.ra.value, icrs.dec.value