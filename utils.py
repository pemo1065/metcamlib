from math import degrees
import ephem

def ra_dec_to_alt_az(ra, dec, lat, lon, elev, timestamp):
    home = ephem.Observer()

    home.lon = str(lon)
    home.lat = str(lat)
    home.elevation = elev
    home.date = timestamp
    body = ephem.FixedBody()
    body._epoch = ephem.J2000
    body._ra = ephem.degrees(str(ra))
    body._dec = ephem.degrees(str(dec))
    body.compute(home)
    return degrees(body.alt), degrees(body.az)


def alt_az_to_ra_dec(alt, az, lat, lon, elev, timestamp):
    home = ephem.Observer()

    home.lon = str(lon)
    home.lat = str(lat)
    home.elevation = elev
    home.date = timestamp
    ra, dec = home.radec_of(az=str(az), alt=str(alt))
    return degrees(ra), degrees(dec)