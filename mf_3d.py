import mf_3d

if __name__ == '__main__':
    mf = mf_3d.MagneticField('Basic_PMS.txt')
    print(mf.get_field(0.1, 0.25, 0.333))