function data = truncNormalGenerator(g)
	data = drawTruncatedNormal([1,1], g.mean, g.sigma, g.mean-2*g.sigma, g.mean+2*g.sigma);
end
