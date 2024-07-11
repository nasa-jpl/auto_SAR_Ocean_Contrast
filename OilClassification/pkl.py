import _pickle as pickle

def load(filename):
    with open(filename,'rb') as fh:
        data = pickle.load(fh)
    return data

def save(filename,data):
    # put in some file protection stuff
    with open(filename,'wb') as fh:
        pickle.dump(data, fh)
    return
