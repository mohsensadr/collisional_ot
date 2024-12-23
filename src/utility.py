import numpy as np


def fit_xp(x, y):
    log_x = np.log10(x)
    log_y = np.log10(y)

    slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)
    slope_str = f'{slope:.2f}'

    fitted_y = 10**(intercept + slope * log_x)
    return slope, slope_str, fitted_y


def swissroll(n=1000, noise=0.5):
    """
    Generate samples of a Swiss roll dataset.

    Parameters:
    n (int)     : Number of points to generate (default: 1000)
    noise (float): Amount of Gaussian noise to add (default: 0.05)

    Returns:
    X (ndarray)  : An n x 2 matrix, where each row is a point (x, y) in 2D space
    color (ndarray): A vector of values corresponding to the angle (theta), useful for coloring
    """
    # Generate random values for the angle (theta) and height (z)
    theta = (3 * np.pi / 2) * (1 + 2 * np.random.rand(n, 1))  # Angle (theta) values for the spiral

    # Generate the (x, y) coordinates
    x = theta * np.cos(theta)
    y = theta * np.sin(theta)

    # Add Gaussian noise
    x += noise * np.random.randn(n, 1)
    y += noise * np.random.randn(n, 1)


    # Combine the coordinates into a matrix X
    X = np.hstack((x / 5., y / 5.))  # Dividing by 5 to match the scale

    return X

def sample_banana(Np):
    x1 = np.random.normal(0,1,(Np,1))
    x2 = x1**2 + np.random.normal(0,1,(Np,1))
    return np.concatenate([x1,x2], axis=1)

class Funnel:
    """
    Defines a multivariate distribution with 2 dimensions that
    have a funnel shape and d-2 dimensional Gaussian components.

    Author:
    Date:   September 2019
    """

    def __init__(self, d, sigma=2.0, limit_min=0.0, limit_max=10.0):
        """
        Initializes the Funnel object.

        Parameters:
        d (int): Dimension of distribution (must be at least 2)
        sigma (float): Standard deviation of the density (default: 2.0)
        limit_min (float): Minimum limit for threshold (default: 0.0)
        limit_max (float): Maximum limit for threshold (default: 10.0)
        """

        if d < 2:
            raise ValueError('Funnel: dimension must be at least 2')

        self.d = d
        self.sigma = sigma
        self.limit_min = limit_min
        self.limit_max = limit_max
        self.name = 'funnel'

    def sample(self, N):
        """
        Samples from the funnel distribution.

        Parameters:
        N (int): Number of samples.

        Returns:
        numpy.ndarray: Samples from the distribution.
        """
        X = np.random.randn(N, self.d)
        X[:, 0] = self.sigma * X[:, 0]
        for i in range(1, self.d):
            X[:, i] = X[:, i] * np.sqrt(self.threshold(X[:, 0]))
        return X

    def threshold(self, x):
        """
        Applies the threshold function to a given value.

        Parameters:
        x (numpy.ndarray): Input array.

        Returns:
        numpy.ndarray: Thresholded values.
        """
        v = np.exp(-x)
        v = np.clip(v, self.limit_min, self.limit_max)
        return v

    def log_pdf(self, X):
        """
        Computes the log of the probability density function.

        Parameters:
        X (numpy.ndarray): Input samples.

        Returns:
        numpy.ndarray: Logarithm of the PDF evaluated at X.
        """
        if X.shape[1] != self.d:
            raise ValueError('Funnel: dimension mismatch for input samples')

        # Log PDF for the first dimension
        log_pi = self.norm_log_pdf(X[:, 0], 0, self.sigma)

        # Log PDF for the remaining dimensions
        for i in range(1, self.d):
            v = self.threshold(X[:, 0])
            log_pi += self.norm_log_pdf(X[:, i], 0, np.sqrt(v))

        return log_pi

    @staticmethod
    def norm_log_pdf(X, mean, sigma):
        """
        Computes the log of the normal probability density function.

        Parameters:
        X (numpy.ndarray): Input samples.
        mean (float): Mean of the distribution.
        sigma (float): Standard deviation of the distribution.

        Returns:
        numpy.ndarray: Logarithm of the normal PDF evaluated at X.
        """
        return -0.5 * np.log(2 * np.pi * sigma**2) - 0.5 * ((X - mean)**2 / sigma**2)

class Ring:
    def __init__(self, d, sigma=0.2, radia=5.0):
        """
        Defines a d-dimensional distribution with two components lying along a ring
        and (d-2) Gaussian components.
        
        :param d: Dimension of distribution
        :param sigma: Standard deviation of density (default: 0.2)
        :param radia: Radial parameter of distribution (default: 5.0)
        """
        if d < 2:
            raise ValueError('Ring: dimension must be at least 2')

        self.d = d
        self.sigma = sigma
        
        # Ensure radia is always a numpy array, even if it's a scalar
        if np.isscalar(radia):
            self.radia = np.array([radia])
        else:
            self.radia = np.array(radia)

        self.name = 'ring'

    def sample(self, N):
        """
        Samples N points from the distribution.

        :param N: Number of samples
        :return: A sample array of size (N, d)
        """
        # Draw angles and noise
        angles = np.random.rand(N) * 2 * np.pi
        noise = np.random.randn(N) * self.sigma

        # Compute weights
        weights = 2 * np.pi * self.radia
        weights = weights / np.sum(weights)

        # Sample from possible radia vector
        radia_idx = np.random.choice(len(self.radia), N, p=weights)
        radius_samples = self.radia[radia_idx] + noise

        # Compute points in first two dimensions (X and Y)
        xs = radius_samples * np.sin(angles)
        ys = radius_samples * np.cos(angles)
        X = np.column_stack([xs, ys])

        # Append independent Gaussian directions for d > 2
        if self.d > 2:
            gauss_noise = self.sigma * np.random.randn(N, self.d - 2)
            X = np.column_stack([X, gauss_noise])

        return X
    def log_pdf(self, X):
        """
        Evaluates the log-probability density function for the given samples X.

        :param X: Input samples of shape (N, d)
        :return: Log-probabilities of shape (N,)
        """
        if X.shape[1] != self.d:
            raise ValueError('Ring: dimension mismatch for input samples')

        # Compute weights
        weights = 2 * np.pi * self.radia
        weights = weights / np.sum(weights)

        # Compute norms in the first two directions
        norms = np.sqrt(np.sum(X[:, :2] ** 2, axis=1))

        # Compute log_pi for the first two dimensions
        log_pi = np.zeros(X.shape[0])
        for i in range(X.shape[0]):
            log_pdf_ring_comps = -0.5 * ((norms[i] - self.radia) ** 2) / (self.sigma ** 2) - \
                                 0.5 * np.log(2 * np.pi * self.sigma ** 2) - \
                                 np.log(2 * np.pi * self.radia)
            log_pi[i] = self.log_sum_exp(log_pdf_ring_comps + np.log(weights))

        # Add the log_pdf for the remaining Gaussian directions (if d > 2)
        if self.d > 2:
            for i in range(2, self.d):
                log_pi += self.norm_log_pdf(X[:, i], 0, self.sigma)

        return log_pi

    def log_sum_exp(self, X):
        """
        Evaluates log-sum-exp for a vector X.

        :param X: Input vector
        :return: log(sum(exp(X)))
        """
        X_min = np.min(X)
        X_without_min = np.delete(X, np.argmin(X))

        # Compute log-sum-exp
        return X_min + np.log(1 + np.sum(np.exp(X_without_min - X_min)))

    def norm_log_pdf(self, X, mean, sigma):
        """
        Evaluates the log of a normal probability density function.

        :param X: Input data
        :param mean: Mean of the normal distribution
        :param sigma: Standard deviation of the normal distribution
        :return: Log-probability of the input data
        """
        return -0.5 * np.log(2 * np.pi * sigma ** 2) - 0.5 / (sigma ** 2) * (X - mean) ** 2
