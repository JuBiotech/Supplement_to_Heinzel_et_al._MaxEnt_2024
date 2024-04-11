function data = absoluteUniformGenerator(g)
	data = rand(1, 1) * (g.bounds(2) - g.bounds(1)) + g.bounds(1);
end
