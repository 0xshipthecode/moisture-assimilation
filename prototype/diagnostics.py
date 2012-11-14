

from datetime import datetime
import cPickle


class Diagnostics:
    """
    Class that stores diagnostics pushed in by computational parts of the system,
    is able to print these out right away and/or store them.
    """
    
    def __init__(self, out_file):
        """
        Initialize the output channel, if any.
        """
        if out_file is not None:
            self.f = open(out_file, 'w')
            self.f.write('This file is [%s] opened on [%s].\n' % (out_file, datetime.now()))
        else:
            self.f = None
            
        self.cfg = {}
        self.store = {}
        
        
    def configure_tag(self, tag, stdout = True, fout = True, store = False):
        """
        Configure the diagnostics behavior of a given tag.
        """
        self.cfg[tag] = (stdout, fout, store)
        
        
    def get_tag_config(self, tag):
        """
        Retrieve the tag configuration from the internal dictionary.
        If none found, return default: stdout = True, fout = True, store = False.
        """
        return self.cfg[tag] if tag in self.cfg else (True, True, False)
        

    def push(self, tag, info):
        """
        Push the information marked by tag into the diagnostics system.
        Pushing performs the required actions according to the behavior
        """ 
        tcfg = self.get_tag_config(tag)
        
        if tcfg[0]:
            print("Tag [%s]\n%s" % (tag, str(info)))
            
        if tcfg[1] and self.f is not None:
            self.f.write("Tag [%s]\n%s\n" % (tag, str(info)))
            self.f.flush()
            
        if tcfg[2]:
            lst = self.store[tag] if tag in self.store else []
            lst.append(info)
            self.store[tag] = lst
            
            
    def pull(self, tag):
        """
        Pull the list of stored values for tag (if any).
        If nothing is stored for the tag, method returns None.
        """
        return self.store[tag] if tag in self.store else None
    
    
    def dump_store(self, fname):
        """
        Dump the store of this diagnostics object to file given as filename.
        """
        with open(fname, 'w') as f:
            cPickle.dump(self.store, f)  



_diagnostic_singleton = None


def diagnostics():
    return _diagnostic_singleton

def init_diagnostics(fname = None):
    global _diagnostic_singleton
    _diagnostic_singleton = Diagnostics(fname)

