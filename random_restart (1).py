import start
import random
def random_restart(myqthreshold,mypidthreshold, qr_align):
    step1 = random.uniform(-myqthreshold, 1.0 -myqthreshold)
    step2 = random.uniform(-mypidthreshold, 1.0 - mypidthreshold)
    step3 = random.uniform(-qr_align, 1.0 - qr_align)
    return (step1+myqthreshold, step2+mypidthreshold, step3+qr_align)