project('expalgo-csr', 'cpp')

omp = dependency('openmp')

executable('transpose',
	'main.cpp',
	override_options: ['cpp_std=c++17'],
	dependencies: omp,
	install: true)