# storing the data from the 6_25_18 mpfd irradiation

class Foil(object):
    def __init__(self, peak_area, peak_error, activity, activity_error, 
                 counting_time):
        self.peak = peak_area
        self.peak_error = peak_error
        self.activity = activity
        self.t_c = counting_time
        self.activity_error = activity_error


data = {}
data['al1'] = Foil(115.3, 14.62, 3.2033E-002, 12.68, 3600)
data['al2'] = Foil(150.9, 16.86, 4.1906E-002, 11.17, 3600)
data['au1'] = Foil(137326.1, 382.43, 1.318e+000, 7.382e-003, 60)
data['au2'] = Foil(175561.4, 433.84, 1.685E+000 , 9.186E-003, 60)
data['au3'] = Foil(185317.2, 447.60, 1.518E+001, 5.480E-001, 60)
data['au4'] = Foil(113048.9, 505.68, 2.938E-001, 3.235E-002, 60)
data['au_cd1'] = Foil(847634.3, 955.40, 1.388E+001, 5.004E-001, 300)
data['au_cd2'] = Foil(1202112.1, 1101.06, 2.308E+000, 1.142E-002, 300)
data['au_cd3'] = Foil(905406.3, 952.34, 1.739E+000, 8.645E-003, 300)
data['au_cd4'] = Foil(836461.0, 4653.72, 2.1532E+002, 6.48, 300)
data['in1'] = Foil(16385.8, 137.96, 9.897E-001, 2.669e-002, 60)
data['in2'] = Foil(18046.7, 144.63, 1.090E+000, 2.926E-002, 60)
data['in3'] = Foil(14003.7, 126.65, 8.458E-001, 2.298E-002, 60)
data['in4'] = Foil(13549.7, 125.62, 8.184E-001, 2.230E-002, 60)
data['in_cd1'] = Foil(40613.9, 219.35, -1, -1, 300)
data['in_cd3'] = Foil(3993.0, 67.70, 4.823E-002, 1.482E-003, 300)