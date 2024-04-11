function data = absoluteUniformIntegerGenerator(g)
	data = randi(g.bounds(2) - g.bounds(1) + 1) + g.bounds(1) - 1;
end
