from __future__ import division
import numpy as np
from scipy.stats import norm


def h2_weights(ldScores, N, M, hsq):
	'''
	Computes appropriate regression weights to correct for heteroskedasticity in the LD 
	Score regression under an infinitesimal model. These regression weights are 
	approximately equal to the reciprocal of the conditional variance function
	1 / var(chi^2 | LD Score)
	
	Parameters
	----------
	M : int > 0
		Number of SNPs used for estimating LD Score (need not equal number of SNPs included in
		the regression).
	N :  np.ndarray of ints > 0 with shape (M, )
		Number of individuals sampled for each SNP.
	ldScores : np.array 
		LD Scores. 
	hsq : float in [0,1]
		Heritability estimate.
	
	Returns
	-------
	weights : np.array
		Regression weights. Approx equal to reciprocal of conditional variance function.
	
	'''
	ldScores = np.fmax(ldScores, 1)
	c = hsq * N / M
	weights = 1 / (1+c*ldScores)**2
	return weights
	
	
def gencov_weights(M, ldScores, N1, N2, No, h1, h2, rho_g, rho):
	'''
	Computes appropriate regression weights to correct for heteroskedasticity in the 
	bivariate LDScore regression under and infinitesimal model. These regression weights are 
	approximately equal to the reciprocal of the conditional variance function
	1 / var(betahat1*betahat2 | LD Score)
	
	Parameters
	----------
	M : int > 0
		Number of SNPs used for estimating LD Score (need not equal number of SNPs included in
		the regression).
	ldScores : np.array 
		LD Scores. 
	rhog : float in [0,1]
		Genetic covariance estimate.
	rho : float in [0,1]
		Phenotypic correlation estimate.
	N1 : np.ndarray of ints > 0 with shape (M, )
		Number of individuals sampled for each SNP in study 1.
	N2 : np.ndarray of ints > 0 with shape (M, )
		Number of individuals sampled for each SNP in study 2.
	No : np.ndarray of ints > 0 with shape (M, )
		Number of overlapping individuals sampled for each SNP.
	
	Returns
	-------
	weights : np.array
		Regression weights. Approx equal to reciprocal of conditional variance function.
	
	'''
	ldScores = np.fmax(ldScores, 1) # set negative LD Score estimates to one
	a = h1*ldScores / M + (1-h1) / N1 
	b = h2*ldScores / M + (1-h2) / N2
	c = rho_g*ldScores/M + No*rho/(N1*N1)
	weights = a*b + 2*c**2
	return 1 / weights


def h2g(chisq, ldScores, N, M, block_size=1000):
	'''
	Current best-practice function for estimating h2g / partitioned h2g via LD Score.
	
	If the GWAS in question is an ascertained case-control sample, N = total sample size
	and the returned value of h2g (and the corresponding standard errors!) will be on the 
	observed scale. Use obs_to_liab along with an estimate of prevalence in order to 
	transform to the liability scale. 
	
	Parameters
	----------
	chisq : np.ndarray of floats with shape (M, )
		Chi-square statistics from GWAS.
	ldScores : np.matrix of floats with shape (M, # of annotations)
		LD Scores.
	N : np.ndarray of ints > 0
		Sample size. Use (# cases) + (# controls) for case/control GWAS. The case/control
		ratio is taken into account when converting to the liability scale.
	M : np.matrix of ints > 0
		Number of SNPs per annotation in the reference panel.
	
	Returns
	-------
	x : LstsqJackknife
		Jackknife object with h2g estimate on the observed scale.
	
	'''
	# aggregate estimate of h2, used only in regression weights
	agg_h2 = (np.mean(chisq) - 1) / np.mean(N*ldScores / M)
	weights = h2_weights(ldscores, N, M, agg_h2)
	x = ldscore_reg(chisq, ldScores, weights=weights, block_size=block_size)
	return x


def gencov(betahat1, betahat2, ldScores, N1, N2, M, N_overlap=None, rho=None, 
	block_size=1000):
	'''
	Current best-practice function for estimating genetic covariance / partitioned genetic
	covariance via LD Score.
	
	Genetic covariance is not on a particularly interpretable scale, especially when one or
	both studies are ascertained case/control sample. However, genetic correlation is always
	on an interpretable scale (e.g., gencov_obs / sqrt(h2_obs1 * h2_obs2) is on the 
	un-ascertained liability scale, because the observed-to-liability factors in the 
	numerator and denominator cancel. 
	
	Genetic covariance is more useful for significance testing, because the estimate of 
	genetic correlation is a ratio estimator, and obtaining reliable confidence intervals 
	for ratio estimators is tricky.
	
	Parameters
	----------
	betahat1 : np.ndarray of floats with shape (M, )
		Effect-size estimates from study 1.
	betahat2 : np.ndarray of floats with shape (M, )
		Effect-size estimates from study 2.
	ldScores : np.matrix of floats with shape (M, # of annotations)
		LD Scores.
	N1 : np.ndarray of ints > 0
		Sample size of study 1. Use (# cases) + (# controls) for case/control GWAS. The case/
		control ratio is taken into account when converting to the liability scale.
	N2 : np.ndarray of ints > 0
		Sample size of study 2.
	M : np.matrix of ints > 0
		Number of SNPs per annotation in the reference panel.
	N_overlap : np.ndarray of ints > 0, default = 0.
		Number of individuals included in both studies. If N_overlap is not zero, then the LD
		Score regression estimates of genetic covariance will be unbiased
	
	Returns
	-------
	x : LstsqJackknife
		Jackknife object with genetic covariance estimate on observed scale.

	'''	
	betaprod = betahat1*betahat2
	agg_gencov = np.mean(betaprod)/ np.mean(ldScores / M)
	weights = gencov_weights(M, ldScores, N1, N2, No, h1, h2, agg_gencov, rho)
	x = ldscore_reg(betaprod, ldScores, weights=weights, block_size=block_size)
	return x

	
def gencor(betahat1, betahat2, ldScores, N1, N2, M, N_overlap=None, rho=None, 
	block_size=1000):
	'''
	Current best-practice function for estimating genetic correlation via LD Score.
	
	Parameters
	----------
	Same as gencov()
	
	Returns
	-------
	TBD
	
	'''
	hsq1 = h2g(N1*betahat1**2, ldScores, N1, M)
	hsq2 = h2g(N1*betahat2**2, ldScores, N2, M)
	cov = gencov(betahat1, betahat2, N1, N2, M, N_overlap, rho, block_size)
	biased_cor = cov.est / np.sqrt(hsq1.est * hsq2.est)
	numer_delete_vals = cov.delete_cals
	denom_delete_vals = np.sqrt(hsq1.delete_vals*hsq2.delete_vals)
	cor = RatioJackknife(biased_cov, cov.delete_vals, numer_delete_vals, denom_delete_vals)
	return (hsq1, hsq2, cov, cor)
	

def ldscore_reg(y, ldScores, weights=None, block_size=1000):
	'''
	Function to estimate heritability / partitioned heritability / genetic covariance/ 
	partitioned genetic covariance from summary statistics. 
	
	NOTE: Weights should be 1 / variance.

	Parameters
	----------
	y : np.matrix
		Response variable (additive chi-square statistics if estimating h2, dominance
		deviation chi-square statistics if estimating H2[DOMDEV], betahat1*betahat2 if
		estimating genetic covariance).
	ldScores : np.matrix
		LD Scores or partitioned LD Scores.
	weights : np.matrix
		Regression weights.
	block_size : int > 0
		Size of blocks for block jackknife standard error. 
	
	Returns
	-------
	output : LinearJackknife
		LD Score regression parameter + standard error estimates.
		NOTE: these are the LD Score regression parameter estimates, *NOT* hsq or genetic 
		covariance estimates. Have to multiply by M/N or M, respectively. 
	
	'''
	if len(ldScores.shape) <= 1:
		ldScores = np.atleast_2d(ldScores).T
	if len(y.shape) > 1 or len(y.shape) == 0:
		raise ValueError('y must have shape (M, )')
	else:
		y = np.atleast_2d(y).T
	
	num_snps = y.shape[0]
	num_annots = ldScores.shape[1]
	if weights is None:
		weights = np.ones(num_snps)
		
	sqrtWeights = np.atleast_2d(np.sqrt(weights)).T
	y = y * sqrtWeights
	x = np.zeros((len(ldScores), num_annots + 1))
	x[:,0:num_annots] = ldScores * sqrtWeights
	x[:,num_annots] = np.squeeze(sqrtWeights) # intercept 
	x = np.matrix(x); y = np.matrix(y)
	output = LstsqJackknife(x, y, block_size)
	return output


def obs_to_liab(h2_obs, P, K):
	'''
	Converts heritability on the observed scale in an ascertained sample to heritability 
	on the liability scale in the population.

	Parameters
	----------
	h2_obs : float	
		Heritability on the observed scale in an ascertained sample.
	P : float in [0,1]
		Prevalence of the phenotype in the sample.
	K : float in [0,1]
		Prevalence of the phenotype in the population.
		
	Returns
	-------
	h2_liab : float
		Heritability of liability in the population and standard error.
		
	'''
	if K <= 0 or K >= 1:
		raise ValueError('K must be in the range (0,1)')
	if P <= 0 or P >= 1:
		raise ValueError('P must be in the range (0,1)')
	
	thresh = norm.isf(K)
	z = norm.pdf(thresh)
	conversion_factor = K**2 * (1-K)**2 / (P * (1-P) * z**2)
	h2_liab = h2_obs * conversion_factor
	return h2_liab
	
	
class LstsqJackknife(object):
	'''
	Least-squares block jackknife.
	
	Terminology 
	-----------
	For a delete-k block jackknife with nb blocks, define
	
	full estimate : 
		The value of the estimator applied to all the data.
	i_th block-value : 
		The value of estimator applied to the ith block of size k.
	i_th delete-k value : 
		The value of the estimator applied to the data with the ith block of size k removed.
	i_th pseudovalue :
		nb*(full estimate) - (nb - 1) * (i_th delete-k value)
	jackknife estimate:
		The mean psuedovalue (with the mean taken over all nb values of i).
	i_th psuedoerror : 
		(i_th pseudovalue) - (jackknife estimate)
	jackknife standard error / variance :
		Standard error / variance of the pseudoerrors		
	output_dim : int > 0
		Number of output parameters. Default is x.shape[1].
	
	Parameters
	----------
	x : np.matrix with shape ()
		Predictors.
	y : np.matrix with shape ()
		Response variable.
	block_size: int, > 0
		Size of jackknife blocks.

	Attributes
	----------
	num_blocks : int 
		Number of jackknife blocks
	est : np.matrix
		Value of estimator applied to full data set.
	pseudovalues : np.matrix
		Jackknife pseudovalues.
	psuedoerrors : np.matrix
		Jackknife pseudoerrors.
	jknife_est : np.matrix
		Mean pseudovalue.,
	jknife_var : np.matrix
		Estimated variance of jackknife estimator.
	jknife_se : np.matrix
		Estimated standard error of jackknife estimator.

	Methods
	-------
	autocov(lag) : 
		Returns lag-[lag] autocovariance in jackknife pseudoerrors.
	autocor(lag) :
		Returns lag-[lag] autocorrelation (autocovariance divided by standard deviation) in 
		jackknife pseudoerrors.
	__block_vals_to_psuedovals__(block_vals, est) :
		Converts block values and full estimate to pseudovalues.
	__block_vals_to_est__() :
		Converts block values to full estimate.
		
	Attributes
	----------
	block_size : int
		Size of all blocks except possibly the last block.
	block_vals : np.matrix
		Block jackknife block values.
	delete_vals : np.matrix
		Block jackknife delete-(block_size) values.
	pseudovalues : np.matrix
		Block jackknife pseudovalues.
	est : np.matrix
		Non-jackknifed estimate.
	jknife_est : np.matrix
		Jackknife estimate.
	jknife_var : np.matrix
	  Jackknife estimate of variance the jackknife estimate.
	jknife_se : np.matrix
		Jackknife estimate of the standard error of the jackknife estimate.
	jknife_cov : np.matrix
		Jackknife estimate of the variance-covariance matrix of the jackknife estimate. 
	
	
	Possible TODO: impute FFT de-correlation (NP)
	'''

	def __init__(self, x, y, block_size):
		if len(x.shape) <= 1:
			x = np.atleast_2d(x).T
		if len(y.shape) <= 1:
			y = np.atleast_2d(y).T
	
		self.N = y.shape[0]
		if self.N != x.shape[0]:
			raise ValueError('Number of data points in y != number of data points in x')
		
		self.output_dim = x.shape[1] 
		if block_size > self.N / 2:
			raise ValueError('Block size must be < N/2')

		self.block_size = block_size		
		self.num_blocks = int(np.ceil(self.N / self.block_size))
		if self.num_blocks > self.N:
			raise ValueError('Number of blocks > number of data points')

		self.block_vals = self.__compute_block_vals__(x, y, block_size)
		self.est = self.__block_vals_to_est__(self.block_vals)
		(self.pseudovalues, self.delete_values) = self.__block_vals_to_pseudovals__(self.block_vals, self.est)
		(self.jknife_est, self.jknife_var, self.jknife_se, self.jknife_cov) = self.__jknife__(self.pseudovalues, self.num_blocks)
		
	def __jknife__(self, pseudovalues, num_blocks):
		jknife_est = np.mean(pseudovalues, axis=0) 
		jknife_var = np.var(pseudovalues, axis=0) / (num_blocks - 1) 
		jknife_se = np.std(pseudovalues, axis=0) / np.sqrt(num_blocks - 1)
		jknife_cov = np.cov(pseudovalues.T) / (num_blocks )
		return (jknife_est, jknife_var, jknife_se, jknife_cov)

	def __compute_block_vals__(self, x, y, block_size):
		xty_block_vals = []; xtx_block_vals = []
		for s in xrange(0, self.N, block_size):	
			# s = block start SNP index; e = block end SNP index
			e = min(self.N, s+block_size)
			xty = np.dot( x[s:e,...].T, y[s:e,...] )
			xtx = np.dot( x[s:e,...].T, x[s:e,...] )
			xty_block_vals.append(xty)
			xtx_block_vals.append(xtx)
			
		block_vals = (xty_block_vals, xtx_block_vals)
		return block_vals

	def __block_vals_to_pseudovals__(self, block_vals, est):
		pseudovalues = np.matrix(np.zeros((self.num_blocks, self.output_dim)))
		delete_values = np.matrix(np.zeros((self.num_blocks, self.output_dim)))
		xty_blocks = block_vals[0]
		xtx_blocks = block_vals[1]
		xty_tot = np.sum(xty_blocks, axis=0)
		xtx_tot = np.sum(xtx_blocks, axis=0)
		for j in xrange(0,self.num_blocks):
			delete_xty_j = xty_tot - xty_blocks[j]
			delete_xtx_inv_j = np.linalg.inv(xtx_tot - xtx_blocks[j])
			delete_value_j = np.dot(delete_xtx_inv_j, delete_xty_j).T
			pseudovalues[j,...] = self.num_blocks*est - (self.num_blocks-1)*delete_value_j
			delete_values[j,...] = delete_value_j

		return (pseudovalues, delete_values)
		
	def __block_vals_to_est__(self, block_vals):
		xty_blocks = block_vals[0]
		xtx_blocks = block_vals[1]
		xty = np.sum(xty_blocks, axis=0)
		xtx_inv = np.linalg.inv(np.sum(xtx_blocks, axis=0))
		return np.matrix(np.dot(xtx_inv, xty).T)
		
	def autocov(self, lag):
		pseudoerrors = self.pseudovalues - np.mean(self.pseudovalues)
		if lag <= 0 or lag >= self.num_blocks:
			error_msg = 'Lag >= number of blocks ({L} vs {N})'
			raise ValueError(error_msg.format(L=lag, N=self.num_blocks))

		v = pseudoerrors[lag:len(pseudoerrors),...]
		w = pseudoerrors[0:len(pseudoerrors) - lag,...]
		autocov = np.diag(np.dot(v.T, w)) / (self.num_blocks - lag)
		return autocov

	def autocor(self, lag):
		return self.autocov(lag) / self.jknife_se
		

class RatioJackknife(LstsqJackknife):
	'''
	Block jackknife class for genetic correlation estimation.
	
	Inherits from LstsqJackknife. 1-D only
	
	Parameters
	----------
	est : float or np.array with shape (# of ratios, )
		(Biased) ratio estimate (e.g., if we are estimate a = b / c, est should be \
		\hat{a} = \hat{b} / \hat{c}.
	numer_delete_vals : np.matrix with shape (# of blocks, # of ratios) or (# of blocks, )
		Delete-k values for the numerator.
	denom_delete_vals: np.matrix with shape (# of blocks, # of ratios) or (# of blocks, )
		Delete-k values for the denominator.
		
		
	'''

	def __init__(self, est, numer_delete_vals, denom_delete_vals):
		if len(numer_delete_vals.shape) <= 1:
			numer_delete_vals = np.atleast_2d(numer_delete_vals).T
		if len(denom_delete_vals.shape) <= 1:
			denom_delete_vals = np.atleast_2d(denom_delete_vals).T
		if numer_delete_vals.shape != denom_delete_vals.shape:
			raise ValueError('numer_delete_vals.shape != denom_delete_vals.shape.')
	
		self.est = np.atleast_1d(np.array(est))
		self.numer_delete_vals = numer_delete_vals 
		self.denom_delete_vals = denom_delete_vals 
		self.num_blocks = numer_delete_vals.shape[0]
		self.output_dim = numer_delete_vals.shape[1]
		self.pseudovalues = self.__delete_vals_to_pseudovals__(self.est, self.denom_delete_vals, self.numer_delete_vals)
		(self.jknife_est, self.jknife_var, self.jknife_se, self.jknife_cov) = self.__jknife__(self.pseudovalues, self.num_blocks)

	def __delete_vals_to_pseudovals__(self, est, denom, numer):
		pseudovalues = np.matrix(np.zeros((self.num_blocks, self.output_dim)))
		for j in xrange(0,self.num_blocks):
			pseudovalues[j,...] = self.num_blocks*est - (self.num_blocks-1)*numer[j,...]/denom[j,...]

		return pseudovalues