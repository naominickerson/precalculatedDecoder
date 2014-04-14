## Reads in the data from its stored form


def syndrome(filename):

    f=open(filename,"r")      
    data= [[[int(z) for z in y.split(" ") if z!=''] for y in x.split(",") if y!='\n'] for x in  f.readlines()]
    f.close()
    return data



    
def results(filename):
    f=open(filename,"r")
    data = [[float(y) for y in (x.rstrip("\n")).split(" ")] for x in f.readlines()]
    f.close()
    return data
    

#results2D("3results_y.txt")



