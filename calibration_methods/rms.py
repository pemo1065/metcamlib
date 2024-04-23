
from RMS.Astrometry.ApplyAstrometry import xyToRaDecPP, raDecToXYPP
from RMS.Formats import Platepar

class Calibration:
    def __init__(self, platepar_file):
        self.platepar: Platepar = Platepar.Platepar()
        self.platepar.read(platepar_file, "json")

    def calibrate(self, img_stars, cat_stars):
        self.platepar.fitAstrometry(self.platepar.JD, img_stars, cat_stars, first_platepar_fit=True, fit_only_pointing=False, fixed_scale=False)
        self.platepar.fitAstrometry(self.platepar.JD, img_stars, cat_stars, first_platepar_fit=False, fit_only_pointing=False, fixed_scale=False)
        return
    
    def get_star_list(self):
        return [(s[1], s[2], s[4], s[5]) for s in self.platepar.star_list]

    def ra_dec_to_xy(self, ra, dec):
        return raDecToXYPP(ra, dec, self.platepar.JD, self.platepar)
    
    def xy_to_ra_dec(self, x, y):
        return xyToRaDecPP(len(x)*[self.platepar.JD], x, y, len(x)*[1], self.platepar, extinction_correction=False, jd_time=True)
