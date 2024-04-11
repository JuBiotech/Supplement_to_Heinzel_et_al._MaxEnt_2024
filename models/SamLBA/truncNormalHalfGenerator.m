function data = truncNormalHalfGenerator(g)
	data = drawTruncatedNormal([1,1], g.mean, g.sigma / 2, g.mean-g.sigma, g.mean+g.sigma);
end
