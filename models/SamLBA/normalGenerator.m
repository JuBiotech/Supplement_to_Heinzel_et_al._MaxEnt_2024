function data = normalGenerator(g)
	data = randn(1, 1) .* g.sigma + g.mean;
end
