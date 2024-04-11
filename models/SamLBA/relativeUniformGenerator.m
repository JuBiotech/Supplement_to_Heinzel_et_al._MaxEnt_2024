function data = relativeUniformGenerator(g, multiplicator)
	data = multiplicator .* 2 .* (rand(1, 1)-0.5) .* g.sigma + g.mean;
end
