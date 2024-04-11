function data = absoluteUniformNonNullIntegerGenerator(g)
	data = randi(g.bounds(2) - g.bounds(1) + 1) + g.bounds(1) - 1;
	while (data == 0)
		data = randi(g.bounds(2) - g.bounds(1) + 1) + g.bounds(1) - 1;
	end
end