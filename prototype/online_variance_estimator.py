

class OnlineVarianceEstimator:
    """
    This class keeps an estimate of the running mean and variance of a field.
    Online algorithm taken from wikipedia [attributed to D. Knuth]
    http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    """
    
    def __init__(self, imean, ivar, iN):
        """
        Initialize with prior information.  
        """
        self.mean = imean
        self.M2 = ivar
        self.N = iN
        

    def update_with(self, ndata):
        """
        Acquire new sample and update the field statistics.
        """
        self.N += 1
        delta = ndata - self.mean
        self.mean += delta / self.N
        self.M2 += delta * (ndata - self.mean)
         
         
    def get_variance(self):
        """
        Return the current variance estimate.
        """
        return self.M2 / max(self.N - 1, 1)
    

    def get_mean(self):
        """
        Returns the estimate of the current mean.
        """
        return self.mean
        
