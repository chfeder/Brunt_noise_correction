
'''
notes:
'''

def Produce_B(sigmasqrho,sigmasqvel):
    # compute the driving parameter
    # the arguments passed in squared values
    drive_param = np.sqrt(sigmasqrho)/np.sqrt(sigmasqvel)

    if 1:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\ndrive %f"%drive_param)
        alphaFile.close()

    return
